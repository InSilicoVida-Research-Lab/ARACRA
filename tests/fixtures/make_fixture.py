#!/usr/bin/env python3
"""
make_fixture.py — build a tiny, deterministic dose-response dataset.

Why synthetic: the real regression risk in ARACRA is the statistics
(DESeq2 -> DRomics -> BMD -> tPOD), not the alignment. Alignment needs 40 GB of
reference and hours of compute; the statistics need a count matrix. So the
fixture enters through Direct Mode and exercises everything downstream of
quantification in a couple of minutes.

Design notes that matter for what the test can assert:

  * Fixed seed, so the count matrix is byte-identical on every machine. A
    change in output therefore means a change in the CODE, not sampling noise.
  * Counts are negative binomial, which is what DESeq2 assumes, so the fit is
    not fighting the generative model.
  * Three gene classes are planted deliberately:
      - responders   : Hill curves with BMDs INSIDE the tested range
      - low_movers   : responses so steep they land BELOW the lowest dose,
                       which is what triggers the NTP extrapolation flag
      - null         : no dose dependence
    The low_movers exist specifically so --bmd_extrap_filter can be tested in
    both states and shown to change the tPOD.
  * Dose spacing is half-log, 8 levels plus vehicle, which is more than
    DRomics' 3-4 level minimum and typical of a TempO-Seq concentration series.
"""

import argparse
import numpy as np
import pandas as pd

SEED = 20260914
DOSES = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]   # µM
N_REPS = 3

N_NULL = 420
N_RESPONDER = 60
N_LOW_MOVER = 20
BASE_MEAN = 300.0        # TempO-Seq-like per-probe depth
DISPERSION = 0.15


def hill(dose, ec50, hill_n, amplitude):
    """Monotone Hill response; dose 0 gives exactly the baseline."""
    d = np.asarray(dose, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(d > 0, d**hill_n / (d**hill_n + ec50**hill_n), 0.0)
    return 1.0 + amplitude * frac


def nb_counts(rng, mean, dispersion):
    """Negative binomial via gamma-Poisson, matching the DESeq2 model."""
    mean = np.maximum(mean, 1e-6)
    shape = 1.0 / dispersion
    scale = mean * dispersion
    return rng.poisson(rng.gamma(shape, scale))


def build(outdir, chemical="TestChem"):
    rng = np.random.default_rng(SEED)

    samples, dose_col, type_col, rep_col = [], [], [], []
    for d in DOSES:
        for r in range(1, N_REPS + 1):
            tag = "VEH" if d == 0 else f"D{str(d).replace('.', 'p')}"
            samples.append(f"{tag}_R{r}")
            dose_col.append(d)
            type_col.append("control" if d == 0 else "treatment")
            rep_col.append(r)

    genes, curves, truth = [], [], []

    # 1. Null genes — flat.
    for i in range(N_NULL):
        g = f"ENSG_NULL_{i:04d}"
        genes.append(g)
        curves.append(np.ones(len(DOSES)))
        truth.append((g, "null", np.nan))

    # 2. Responders — EC50 inside the tested range (0.1-10 µM).
    for i in range(N_RESPONDER):
        g = f"ENSG_RESP_{i:04d}"
        ec50 = float(10 ** rng.uniform(-1.0, 1.0))
        n = float(rng.uniform(1.2, 3.0))
        amp = float(rng.choice([-1, 1]) * rng.uniform(0.8, 2.5))
        genes.append(g)
        curves.append(hill(DOSES, ec50, n, amp))
        truth.append((g, "responder", ec50))

    # 3. Low movers — EC50 far below the lowest tested dose (0.01 µM). These
    #    produce BMDs the NTP extrapolation rule is meant to catch.
    for i in range(N_LOW_MOVER):
        g = f"ENSG_LOWMOVE_{i:04d}"
        ec50 = float(10 ** rng.uniform(-4.5, -3.2))
        n = float(rng.uniform(1.5, 3.5))
        amp = float(rng.choice([-1, 1]) * rng.uniform(1.0, 2.0))
        genes.append(g)
        curves.append(hill(DOSES, ec50, n, amp))
        truth.append((g, "low_mover", ec50))

    curves = np.vstack(curves)                      # genes x doses
    per_gene_base = rng.lognormal(np.log(BASE_MEAN), 0.5, size=len(genes))

    counts = np.zeros((len(genes), len(samples)), dtype=int)
    col = 0
    for di in range(len(DOSES)):
        for _ in range(N_REPS):
            mu = per_gene_base * curves[:, di]
            counts[:, col] = nb_counts(rng, mu, DISPERSION)
            col += 1

    counts_df = pd.DataFrame(counts, index=genes, columns=samples)
    counts_df.index.name = "gene_id"

    meta_df = pd.DataFrame({
        "Sample_Name": samples,
        "Treatment": [chemical if t == "treatment" else "Control" for t in type_col],
        "Dose": dose_col,
        "Type": type_col,
        "Batch": [1] * len(samples),
        "Replicate": rep_col,
    })

    truth_df = pd.DataFrame(truth, columns=["gene_id", "class", "true_ec50_uM"])

    counts_df.to_csv(f"{outdir}/counts_matrix.csv")
    meta_df.to_csv(f"{outdir}/metadata.csv", index=False)
    truth_df.to_csv(f"{outdir}/ground_truth.csv", index=False)

    print(f"  genes            : {len(genes)} "
          f"({N_NULL} null, {N_RESPONDER} responder, {N_LOW_MOVER} low-mover)")
    print(f"  samples          : {len(samples)} "
          f"({len(DOSES)} dose levels x {N_REPS} reps)")
    print(f"  doses (uM)       : {', '.join(str(d) for d in DOSES)}")
    print(f"  tested range     : {DOSES[1]} - {DOSES[-1]}")
    print(f"  extrapolation cut: {DOSES[1] / 10:g} (lowest dose / 10)")
    print(f"  expect ~{N_LOW_MOVER} genes flagged as extrapolated")
    print(f"  median depth     : {int(np.median(counts.sum(axis=0)))} counts/sample")
    return counts_df, meta_df, truth_df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--chemical", default="TestChem")
    a = ap.parse_args()
    build(a.outdir, a.chemical)
