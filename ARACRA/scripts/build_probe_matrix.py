#!/usr/bin/env python3
# =============================================================================
#  build_probe_matrix.py — probe-level TempO-Seq count matrix
#
#  WHY THIS EXISTS
#    temposeq_counts.py SUMS probes to gene level before modelling. The EPA
#    HTTr workflow (Harrill et al. 2024, Toxicology 501:153694, §2.5) does NOT
#    do that: it models at the PROBE level with DESeq2, and only afterwards
#    aggregates to gene level by taking "the highest magnitude fold change for
#    any associated probe in either direction".
#
#    Summing probe counts and taking max|L2FC| of probes are different
#    estimators. Running the EPA pipeline on a summed matrix silently
#    substitutes the aggregation rule. So: probe-level in, gene-level out at
#    the L2FC step (done inside run_tcplfit2_epa.R).
#
#  Usage:
#    python3 build_probe_matrix.py \
#        --counts_dir path/to/idxstats \
#        --outfile probe_counts.csv
# =============================================================================

import csv, argparse, sys, json
from pathlib import Path


def parse_idxstats(filepath):
    """samtools idxstats -> {probe_name: mapped_count}"""
    counts = {}
    with open(filepath) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4 or parts[0] == '*':
                continue
            counts[parts[0]] = int(parts[2])
    return counts


def find_count_files(counts_dir):
    counts_dir = Path(counts_dir)
    for pattern in ["*_idxstats.tsv", "*.idxstats", "*_idxstats.*",
                    "*_probe_counts.tsv", "**/*_idxstats.tsv"]:
        files = sorted(counts_dir.glob(pattern))
        if files:
            return files
    print(f"ERROR: No idxstats files found in {counts_dir}", file=sys.stderr)
    return []


def main(counts_dir, outfile):
    files = find_count_files(counts_dir)
    if not files:
        sys.exit(1)

    samples, per_sample = [], {}
    all_probes = set()
    for f in files:
        s = f.stem
        for sfx in ['_idxstats', '_probe_counts', '.idxstats']:
            s = s.replace(sfx, '')
        samples.append(s)
        c = parse_idxstats(f)
        per_sample[s] = c
        all_probes.update(c.keys())

    probes = sorted(all_probes)
    print(f"Probe-level matrix: {len(probes)} probes x {len(samples)} samples")

    with open(outfile, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow([""] + samples)
        for p in probes:
            w.writerow([p] + [per_sample[s].get(p, 0) for s in samples])

    qc = {s: {'total_reads': sum(per_sample[s].values()),
              'probes_detected': sum(1 for v in per_sample[s].values() if v > 0)}
          for s in samples}
    Path(outfile).parent.joinpath("probe_matrix_summary.json").write_text(
        json.dumps({'n_probes': len(probes), 'n_samples': len(samples),
                    'sample_qc': qc}, indent=2))
    print(f"Output: {outfile}")
    for s in samples:
        print(f"  {s}: {qc[s]['total_reads']:,} reads | "
              f"{qc[s]['probes_detected']} probes detected")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Build probe-level count matrix (no gene aggregation)")
    ap.add_argument("--counts_dir", required=True)
    ap.add_argument("--outfile", default="probe_counts.csv")
    a = ap.parse_args()
    main(a.counts_dir, a.outfile)
