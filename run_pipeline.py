"""
scripts/run_pipeline.py
────────────────────────
Run the complete scRNA-seq Scanpy pipeline end-to-end.

Usage:
    python scripts/run_pipeline.py
    python scripts/run_pipeline.py --steps 1 2 3       # run only steps 1-3
    python scripts/run_pipeline.py --start 3            # start from step 3
"""

import sys
import os
import argparse
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.step01_preprocessing   import run_preprocessing
from scripts.step02_normalization   import run_normalization
from scripts.step03_dim_reduction   import run_dim_reduction
from scripts.step04_clustering      import run_clustering
from scripts.step05_diff_expression import run_diff_expression

# Re-import using the script names we actually created
import importlib

STEPS = {
    1: ("Preprocessing & QC",                "scripts.01_preprocessing",   "run_preprocessing"),
    2: ("Normalization & Feature Selection", "scripts.02_normalization",   "run_normalization"),
    3: ("Dimensionality Reduction",          "scripts.03_dim_reduction",   "run_dim_reduction"),
    4: ("Clustering",                        "scripts.04_clustering",      "run_clustering"),
    5: ("Differential Expression",           "scripts.05_diff_expression", "run_diff_expression"),
}


def run_step(step_num: int) -> None:
    name, module_path, func_name = STEPS[step_num]
    print(f"\n{'='*60}")
    print(f"  STEP {step_num}: {name}")
    print(f"{'='*60}")
    t0 = time.time()

    # Dynamic import
    mod = importlib.import_module(module_path)
    func = getattr(mod, func_name)
    func()

    elapsed = time.time() - t0
    print(f"\n  ✓ Step {step_num} finished in {elapsed:.1f}s")


def main():
    parser = argparse.ArgumentParser(
        description="Run scRNA-seq Scanpy pipeline"
    )
    parser.add_argument(
        "--steps", nargs="+", type=int,
        help="Specific steps to run (e.g. --steps 1 2 3)"
    )
    parser.add_argument(
        "--start", type=int, default=1,
        help="Start from this step (default: 1)"
    )
    args = parser.parse_args()

    if args.steps:
        steps_to_run = sorted(set(args.steps))
    else:
        steps_to_run = list(range(args.start, 6))

    print("\n" + "═"*60)
    print("   Single-Cell RNA-seq Pipeline  ·  Scanpy")
    print("═"*60)
    print(f"Running steps: {steps_to_run}\n")

    total_start = time.time()

    for step in steps_to_run:
        if step not in STEPS:
            print(f"[warn] Step {step} does not exist — skipping.")
            continue
        try:
            run_step(step)
        except Exception as e:
            print(f"\n[ERROR] Step {step} failed: {e}")
            print("Pipeline halted. Fix the error and rerun from this step with:")
            print(f"  python scripts/run_pipeline.py --start {step}\n")
            sys.exit(1)

    total_elapsed = time.time() - total_start
    print(f"\n{'═'*60}")
    print(f"  ✅ Pipeline complete in {total_elapsed:.1f}s")
    print(f"{'═'*60}\n")


if __name__ == "__main__":
    main()
