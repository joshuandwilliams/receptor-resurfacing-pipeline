"""A numpy-backed stand-in for the torch operations bin/mpnn_design_region_score.py uses.

torch is not in this repository's `test` extra and is not installed in CI, so
`mpnn_design_region_score.py`'s top-level `import torch` fails there, HAS_MPNN
goes False, and `main()` exits 1 before any of the scoring arithmetic runs.
On a laptop where torch happens to be installed it fails differently: torch
2.2's `Tensor.numpy()` raises "Numpy is not available" against numpy 2.x,
because that wheel was compiled for numpy 1.x. Either way the module's
arithmetic is unreachable, and it is the only part of the module worth
testing.

This shim makes it reachable and makes it reachable identically everywhere. It
is loaded only by `fixtures_mpnn_stage.mpnn_module_with_stub`, which puts this
directory at the front of sys.path for the duration of one `exec_module` and
takes it back off afterwards, so nothing else can pick it up.

What is and is not under test
-----------------------------
The module performs exactly these operations: build a 2-D index tensor from a
sequence, pad or truncate it along dim 1, read `.shape[1]`, squeeze dim 0 and
hand the result to numpy, and index a (L, 21) log-probability array. All of
them are elementwise or shape-only, and numpy's semantics for them are
torch's. What this shim cannot vouch for is anything about autograd, devices
or dtype promotion, and the module does none of that: `no_grad` wraps a
forward pass that this suite stubs anyway, and the device is only ever printed.

So a test built on this exercises the module's own arithmetic — the alphabet
mapping, the pad value, the free-position mask, the two means — and not
torch's. That is the part that has a right answer stated at the call site.

Deliberate narrowness
---------------------
`squeeze(dim)` raises on a dimension that is not size 1, where torch would
return the tensor unchanged. Every squeeze in the module is `squeeze(0)` on a
leading batch dimension of 1, so the stricter behaviour turns an unexpected
shape into a failure instead of a silently different answer.

`load` is `pickle.load`, not torch's zip archive reader, so a checkpoint
fixture is a pickle. `load_model` only ever indexes the mapping that comes
back, so the serialisation format is not under test.
"""
from __future__ import annotations

import contextlib
import pickle

import numpy as np

# dtype markers. The module only names `torch.long`; float32 is here because
# the ProteinMPNN stub beside this file builds float tensors.
long = np.int64
float32 = np.float32


class device:  # noqa: N801 - lowercase to match torch
    """`torch.device("cpu")`. Only its string form is used, by main()'s log."""

    def __init__(self, name: str):
        self.type = name

    def __repr__(self) -> str:
        return self.type

    __str__ = __repr__


class cuda:  # noqa: N801 - lowercase to match torch
    """`torch.cuda.is_available()`, always False: the logic tier has no GPU,
    and the scoring arithmetic does not depend on which device it ran on."""

    @staticmethod
    def is_available() -> bool:
        return False


class Tensor:
    """A thin wrapper over an ndarray, with the methods the module calls."""

    def __init__(self, array):
        self._array = np.asarray(array)

    @property
    def shape(self):
        return self._array.shape

    def squeeze(self, dim: int) -> "Tensor":
        if self._array.shape[dim] != 1:
            raise ValueError(
                f"squeeze({dim}) on a dimension of size "
                f"{self._array.shape[dim]}; the module only ever squeezes a "
                "leading batch dimension of 1")
        return Tensor(np.squeeze(self._array, axis=dim))

    def unsqueeze(self, dim: int) -> "Tensor":
        return Tensor(np.expand_dims(self._array, axis=dim))

    def cpu(self) -> "Tensor":
        return self

    def to(self, _device) -> "Tensor":
        return self

    def numpy(self):
        return self._array

    def __getitem__(self, item) -> "Tensor":
        return Tensor(self._array[item])

    def __len__(self) -> int:
        return len(self._array)

    def __repr__(self) -> str:
        return f"stub.Tensor({self._array!r})"


def tensor(data, dtype=None, device=None) -> Tensor:
    return Tensor(np.array(data, dtype=dtype))


def zeros(*shape, dtype=None, device=None) -> Tensor:
    if len(shape) == 1 and isinstance(shape[0], (tuple, list)):
        shape = tuple(shape[0])
    return Tensor(np.zeros(shape, dtype=dtype or np.float32))


def ones(*shape, dtype=None, device=None) -> Tensor:
    if len(shape) == 1 and isinstance(shape[0], (tuple, list)):
        shape = tuple(shape[0])
    return Tensor(np.ones(shape, dtype=dtype or np.float32))


def full(shape, fill_value, dtype=None, device=None) -> Tensor:
    return Tensor(np.full(tuple(shape), fill_value, dtype=dtype))


def arange(stop, dtype=None, device=None) -> Tensor:
    return Tensor(np.arange(stop, dtype=dtype))


def cat(tensors, dim: int = 0) -> Tensor:
    return Tensor(np.concatenate([t.numpy() for t in tensors], axis=dim))


@contextlib.contextmanager
def no_grad():
    yield


def load(path, map_location=None, weights_only=False):
    with open(path, "rb") as handle:
        return pickle.load(handle)


def save(obj, path) -> None:
    with open(path, "wb") as handle:
        pickle.dump(obj, handle)
