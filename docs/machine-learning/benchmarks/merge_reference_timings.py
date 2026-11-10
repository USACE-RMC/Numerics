"""Merges the scikit-learn and DiceKriging timings into gpr-reference-timings.json.

The merged file is read by the Benchmark test in Test_GaussianProcessRegression.cs, which times the C#
fit and prediction on the same data with the given number of repetitions per size.

Usage:
    python merge_reference_timings.py <sklearn.json> <dicekriging.json> <processor> <out.json> [csharpReps]
"""
import json
import sys


def main():
    sk = json.load(open(sys.argv[1], encoding="utf-8"))
    dk = json.load(open(sys.argv[2], encoding="utf-8"))
    reps = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    if sk["d"] != dk["d"] or sk["m"] != dk["m"]:
        raise ValueError("The two timing files use different dimensions.")
    dk_sizes = {s["n"]: s for s in dk["sizes"]}
    sizes = []
    for s in sk["sizes"]:
        other = dk_sizes[s["n"]]
        sizes.append({
            "n": s["n"], "reps": reps,
            "scikit-learn": {"fitSeconds": s["fitSeconds"], "predictSeconds": s["predictSeconds"], "reps": s["reps"]},
            "DiceKriging": {"fitSeconds": other["fitSeconds"], "predictSeconds": other["predictSeconds"], "reps": other["reps"]},
        })
    out = {
        "processor": sys.argv[3], "d": sk["d"], "m": sk["m"],
        "scikit-learn": {k: sk[k] for k in ("sklearn", "numpy", "scipy", "python", "blas")},
        "DiceKriging": {k: dk[k] for k in ("DiceKriging", "R", "blas")},
        "sizes": sizes,
    }
    with open(sys.argv[4], "w", encoding="utf-8", newline="\n") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
