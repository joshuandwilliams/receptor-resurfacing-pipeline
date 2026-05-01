from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ComparisonResult:
    passed: bool
    reference: Path
    actual: Path
    strategy: str
    message: str = ""
    differences: list[str] = field(default_factory=list)

    def __bool__(self) -> bool:
        return self.passed

    def assert_passed(self) -> None:
        if self.passed:
            return
        lines = [
            f"[{self.strategy}] characterization comparison failed",
            f"  reference: {self.reference}",
            f"  actual:    {self.actual}",
        ]
        if self.message:
            lines.append(f"  message:   {self.message}")
        if self.differences:
            lines.append("  differences:")
            for diff in self.differences:
                lines.append(f"    - {diff}")
        raise AssertionError("\n".join(lines))
