#!/usr/bin/env python3
"""AF3 no-MSA inputs for the unsteered sequences of the never-steered designs.

The pipeline ran AlphaFold3 only on each survivor's *steered* sequence, so
there is no orthogonal prediction of the original ProteinMPNN sequences. That
leaves the question of whether AlphaFold3 and Boltz-2 agree on the designs
Boltz-2 posed correctly without any steering unanswered.

This writes one AF3 JSON per design the pipeline never steered
(cold_start_all_clean in its plan.json), carrying the original ProteinMPNN
receptor sequence rather than the steered one. Everything else matches the
pipeline's AF3_NOMSA_ON_SURVIVORS process: the same three seeds, no MSA, and
the effector template forced on chain B.

Usage:
    generate_inputs.py
"""

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(HERE, "..", "..", "runs", "crystal_full_test_contig",
                   "results")
RUNS = os.path.join(RUN, "negative_steering", "runs")
OUT = os.path.join(HERE, "inputs")
SEEDS = [42, 123, 456]


def never_steered():
    out = []
    for unit in sorted(os.listdir(RUNS)):
        plan = os.path.join(RUNS, unit, "cycle_0", "plan.json")
        if not os.path.isfile(plan) or "control" in unit:
            continue
        with open(plan) as fh:
            if json.load(fh).get("cold_start_all_clean"):
                out.append(unit)
    return out


def main():
    import csv
    os.makedirs(OUT, exist_ok=True)

    seqs = {}
    with open(os.path.join(RUN, "sequences", "scored_metadata.csv")) as fh:
        for r in csv.DictReader(fh):
            seqs[f"design_{r['design']}_seq_{r['seq']}"] = r["corrected_receptor"]

    # Effector sequence and its template come from any existing AF3 input, so
    # this run is identical to the pipeline's except for the receptor sequence.
    ref = os.path.join(RUN, "orthogonal_metrics", "af3_nomsa",
                       "design_0_seq_0", "input.json")
    with open(ref) as fh:
        template = json.load(fh)
    effector = next(c for c in template["sequences"]
                    if c["protein"]["id"] == "B")["protein"]

    units = never_steered()
    for unit in units:
        payload = {
            "name": f"{unit}_unsteered",
            "modelSeeds": SEEDS,
            "dialect": template["dialect"],
            "version": template["version"],
            "sequences": [
                # templates must be present on every chain under
                # --norun_data_pipeline, empty for the untemplated receptor.
                {"protein": {"id": "A", "sequence": seqs[unit],
                             "unpairedMsa": "", "pairedMsa": "",
                             "templates": []}},
                {"protein": {"id": "B", "sequence": effector["sequence"],
                             "unpairedMsa": "", "pairedMsa": "",
                             "templates": effector["templates"]}},
            ],
        }
        with open(os.path.join(OUT, f"{unit}.json"), "w") as fh:
            json.dump(payload, fh, indent=2)

    with open(os.path.join(HERE, "units.txt"), "w") as fh:
        fh.write("\n".join(units) + "\n")
    print(f"wrote {len(units)} AF3 inputs to {OUT}")
    print("units:", ", ".join(units))


if __name__ == "__main__":
    main()
