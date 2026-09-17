#!/usr/bin/env python3
"""
make_pathway_fixture.py — build a dose-response count matrix using REAL
Ensembl gene IDs (from dump_real_pathway.R), to actually exercise
run_dromics.R's pathway-level tPOD code path for the first time.

04_fixture.sh's main fixture uses fake gene IDs that never match GO/KEGG/
MSigDB annotations, so that code path (gene -> pathway mapping, per-pathway
median BMD, "lowest median pathway wins" tPOD) has never run in the test
suite at all. This fixture supplies real IDs so it can.

Two scenarios, same background genes and dosing grid as 04_fixture.sh's
main fixture (so results are comparable):

  null  — every gene (background AND pathway genes) is flat. No true
          dose-response signal anywhere. Characterizes the "noise floor":
          how low a tPOD the NTP "lowest of many tested gene sets" method
          reports from pure chance, given hundreds of real gene sets get
          tested against ~400 background genes with real annotations.

  spike — every gene in the chosen real pathway gets an identical Hill-curve
          response with a known EC50; background genes stay flat. Tests
          whether the algorithm both (a) identifies that exact real pathway
          as the top hit and (b) reports a tPOD close to the true EC50.
"""

import argparse
import json
import numpy as np
import pandas as pd

SEED = 20260916
DOSES = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]   # µM, matches 04_fixture.sh
N_REPS = 3
BASE_MEAN = 300.0
DISPERSION = 0.15
SPIKE_EC50 = 1.0          # µM — comfortably inside the tested range
SPIKE_HILL_N = 2.0
SPIKE_AMPLITUDE = 2.0


def hill(dose, ec50, hill_n, amplitude):
    d = np.asarray(dose, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(d > 0, d**hill_n / (d**hill_n + ec50**hill_n), 0.0)
    return 1.0 + amplitude * frac


def nb_counts(rng, mean, dispersion):
    mean = np.maximum(mean, 1e-6)
    shape = 1.0 / dispersion
    scale = mean * dispersion
    return rng.poisson(rng.gamma(shape, scale))


def build(outdir, pathway_genes, background_genes, scenario, chemical="PathwayTestChem"):
    rng = np.random.default_rng(SEED)

    samples, dose_col, type_col = [], [], []
    for d in DOSES:
        for r in range(1, N_REPS + 1):
            tag = "VEH" if d == 0 else f"D{str(d).replace('.', 'p')}"
            samples.append(f"{tag}_R{r}")
            dose_col.append(d)
            type_col.append("control" if d == 0 else "treatment")

    genes = list(background_genes) + list(pathway_genes)
    is_pathway = np.array([False] * len(background_genes) + [True] * len(pathway_genes))

    if scenario == "null":
        curves = np.ones((len(genes), len(DOSES)))
    elif scenario == "spike":
        flat = np.ones(len(DOSES))
        spiked = hill(DOSES, SPIKE_EC50, SPIKE_HILL_N, SPIKE_AMPLITUDE)
        curves = np.vstack([spiked if p else flat for p in is_pathway])
    else:
        raise ValueError(scenario)

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
    })

    counts_df.to_csv(f"{outdir}/counts_matrix.csv")
    meta_df.to_csv(f"{outdir}/metadata.csv", index=False)

    meta_out = {
        "scenario": scenario,
        "n_genes": len(genes),
        "n_background": len(background_genes),
        "n_pathway": len(pathway_genes),
        "spike_ec50": SPIKE_EC50 if scenario == "spike" else None,
    }
    with open(f"{outdir}/scenario_meta.json", "w") as f:
        json.dump(meta_out, f, indent=2)

    print(f"  scenario         : {scenario}")
    print(f"  genes            : {len(genes)} ({len(background_genes)} background, "
          f"{len(pathway_genes)} pathway)")
    if scenario == "spike":
        print(f"  spike EC50 (uM)  : {SPIKE_EC50}")
    return counts_df, meta_df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--pathway-genes-csv", required=True)
    ap.add_argument("--background-genes-csv", required=True)
    ap.add_argument("--scenario", choices=["null", "spike"], required=True)
    a = ap.parse_args()

    pathway_genes = pd.read_csv(a.pathway_genes_csv)["ensembl_id"].tolist()
    background_genes = pd.read_csv(a.background_genes_csv)["ensembl_id"].tolist()
    build(a.outdir, pathway_genes, background_genes, a.scenario)
