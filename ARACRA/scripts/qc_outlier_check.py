#!/usr/bin/env python3
# =============================================================================
#  Layer 1: Alignment QC outlier detection
#  Reads STAR _Log.final.out or HISAT2 _align_summary.txt files,
#  runs iterative Grubbs + IQR on mapping rates.
#  Both tests must agree before a sample is flagged.
#  Writes outlier_flag.json for the Streamlit popup warning.
#
#  Usage: python3 qc_outlier_check.py --star_dir PATH --outdir PATH
#         [--aligner star|hisat2]
#
#  Changes from v1:
#   - Bidirectional Grubbs test (checks both min AND max extremes)
#   - Iterative Grubbs (removes flagged sample, re-tests until clean)
#   - IQR now checks both lower AND upper fences
#   - Proper quartile calculation (method-7 linear interpolation)
#   - HISAT2 alignment summary parsing support
#   - Warnings for malformed/truncated log files
#   - Robust error handling per-sample
# =============================================================================

import re, json, math, argparse, sys
from pathlib import Path

# Grubbs critical values, two-sided alpha=0.05
GRUBBS = {
    3: 1.155, 4: 1.481, 5: 1.715, 6: 1.887, 7: 2.020, 8: 2.126,
    9: 2.215, 10: 2.290, 11: 2.355, 12: 2.412, 13: 2.462, 14: 2.507,
    15: 2.549, 16: 2.585, 17: 2.620, 18: 2.651, 19: 2.681, 20: 2.709,
    21: 2.733, 22: 2.758, 23: 2.781, 24: 2.802, 25: 2.822,
    30: 2.908, 35: 2.979, 40: 3.036, 50: 3.128, 60: 3.199,
    80: 3.305, 100: 3.383, 120: 3.428, 150: 3.520, 200: 3.607
}


def gcrit(n):
    """Interpolate Grubbs critical value for sample size n."""
    if n < 3:
        return float("inf")
    ks = sorted(GRUBBS)
    if n in GRUBBS:
        return GRUBBS[n]
    if n > ks[-1]:
        return GRUBBS[ks[-1]]
    lo = max(k for k in ks if k <= n)
    hi = min(k for k in ks if k >= n)
    t = (n - lo) / (hi - lo) if hi != lo else 0
    return GRUBBS[lo] + t * (GRUBBS[hi] - GRUBBS[lo])


def quantile_linear(sorted_vals, p):
    """Method-7 linear interpolation (numpy/R default)."""
    n = len(sorted_vals)
    if n == 0:
        return 0
    if n == 1:
        return sorted_vals[0]
    h = (n - 1) * p
    lo = int(math.floor(h))
    hi = min(lo + 1, n - 1)
    frac = h - lo
    return sorted_vals[lo] + frac * (sorted_vals[hi] - sorted_vals[lo])


def parse_star_log(path):
    """Parse STAR _Log.final.out file."""
    pats = {
        "pct_unique":    r"Uniquely mapped reads %\s*\|\s*([\d.]+)%",
        "pct_multi":     r"% of reads mapped to multiple loci\s*\|\s*([\d.]+)%",
        "pct_too_short": r"% of reads unmapped: too short\s*\|\s*([\d.]+)%",
        "total_reads":   r"Number of input reads\s*\|\s*(\d+)",
    }
    try:
        txt = open(path).read()
    except Exception as e:
        print(f"  WARNING: Cannot read {path}: {e}", file=sys.stderr)
        return None

    result = {}
    for k, p in pats.items():
        m = re.search(p, txt)
        if m:
            result[k] = float(m.group(1))

    if "pct_unique" not in result:
        print(f"  WARNING: Could not parse mapping rate from {path.name} "
              f"(file may be truncated or malformed)", file=sys.stderr)
        return None

    return result


def parse_hisat2_log(path):
    """Parse HISAT2 _align_summary.txt file."""
    try:
        txt = open(path).read()
    except Exception as e:
        print(f"  WARNING: Cannot read {path}: {e}", file=sys.stderr)
        return None

    result = {}

    # Total reads
    m = re.search(r"(\d+) reads; of these:", txt)
    if m:
        result["total_reads"] = float(m.group(1))

    # Aligned exactly 1 time (uniquely)
    m = re.search(r"(\d+) \(([\d.]+)%\) aligned exactly 1 time", txt)
    if m:
        result["pct_unique"] = float(m.group(2))

    # Multi-mapped
    m2 = re.search(r"(\d+) \(([\d.]+)%\) aligned >1 times", txt)
    if m2:
        result["pct_multi"] = float(m2.group(2))

    # Fallback: overall alignment rate
    if "pct_unique" not in result:
        m = re.search(r"([\d.]+)% overall alignment rate", txt)
        if m:
            result["pct_unique"] = float(m.group(1))

    if "pct_unique" not in result:
        print(f"  WARNING: Could not parse alignment rate from {path.name} "
              f"(file may be truncated or malformed)", file=sys.stderr)
        return None

    return result


def iqr_test(vals, names):
    """IQR outlier test — checks BOTH lower and upper fences."""
    n = len(vals)
    if n < 3:
        return []

    sv = sorted(vals)
    q1 = quantile_linear(sv, 0.25)
    q3 = quantile_linear(sv, 0.75)
    iqr = q3 - q1
    lower_fence = q1 - 1.5 * iqr
    upper_fence = q3 + 1.5 * iqr

    flagged = []
    for nm, v in zip(names, vals):
        if v < lower_fence or v > upper_fence:
            flagged.append({
                "sample": nm, "value": round(v, 2),
                "lower_fence": round(lower_fence, 2),
                "upper_fence": round(upper_fence, 2),
                "fence": round(lower_fence, 2),  # backward compat
                "q1": round(q1, 2), "q3": round(q3, 2),
                "direction": "low" if v < lower_fence else "high"
            })
    return flagged


def grubbs_test_iterative(vals, names, max_iter=5):
    """
    Iterative bidirectional Grubbs test.
    Checks the most extreme value (min OR max) at each step.
    Removes flagged sample and re-tests until clean or max_iter reached.
    """
    remaining_vals = list(vals)
    remaining_names = list(names)
    flagged = []

    for _ in range(max_iter):
        n = len(remaining_vals)
        if n < 3:
            break

        mu = sum(remaining_vals) / n
        std = math.sqrt(sum((v - mu) ** 2 for v in remaining_vals) / (n - 1))
        if std == 0:
            break

        # Find the most extreme value (bidirectional)
        deviations = [(abs(v - mu), i, v) for i, v in enumerate(remaining_vals)]
        max_dev, max_idx, max_val = max(deviations, key=lambda x: x[0])

        G = max_dev / std
        Gc = gcrit(n)

        if G > Gc:
            flagged.append({
                "sample": remaining_names[max_idx],
                "value": round(max_val, 2),
                "G": round(G, 3),
                "G_crit": round(Gc, 3),
                "p": "<0.05",
                "direction": "low" if max_val < mu else "high"
            })
            remaining_vals.pop(max_idx)
            remaining_names.pop(max_idx)
        else:
            break

    return flagged


def main(star_dir, outdir, aligner="star"):
    sd = Path(star_dir)
    od = Path(outdir)
    od.mkdir(parents=True, exist_ok=True)

    # Discover log files based on aligner
    if aligner == "hisat2":
        logs = sorted(sd.glob("*_align_summary.txt"))
        if not logs:
            logs = sorted(sd.glob("**/*_align_summary.txt"))
        suffix = "_align_summary.txt"
        parse_fn = parse_hisat2_log
    else:
        logs = sorted(sd.glob("*_Log.final.out"))
        if not logs:
            logs = sorted(sd.glob("**/*_Log.final.out"))
        suffix = "_Log.final.out"
        parse_fn = parse_star_log

    if not logs:
        (od / "outlier_flag.json").write_text(json.dumps({
            "layer": "alignment", "status": "no_logs",
            "aligner": aligner, "samples": [], "flagged": []
        }, indent=2))
        print(f"No {aligner.upper()} logs found in {star_dir}")
        return

    # Parse all logs — warn on failures instead of silently skipping
    samples = []
    skipped = []
    for lf in logs:
        m = parse_fn(lf)
        if m is not None and "pct_unique" in m:
            samples.append({"sample": lf.name.replace(suffix, ""), **m})
        else:
            skipped.append(lf.name.replace(suffix, ""))

    if skipped:
        print(f"  WARNING: Skipped {len(skipped)} sample(s) with unparseable logs: "
              f"{', '.join(skipped)}", file=sys.stderr)

    if len(samples) < 3:
        (od / "outlier_flag.json").write_text(json.dumps({
            "layer": "alignment", "status": "too_few_samples",
            "aligner": aligner, "n_samples": len(samples),
            "skipped": skipped, "samples": samples, "flagged": []
        }, indent=2))
        print(f"Only {len(samples)} parseable samples — need >=3 for outlier detection")
        return

    names = [s["sample"] for s in samples]
    rates = [s["pct_unique"] for s in samples]
    libs = [s.get("total_reads", 0) / 1e6 for s in samples]
    mu = sum(rates) / len(rates)
    std_r = (math.sqrt(sum((r - mu) ** 2 for r in rates) / (len(rates) - 1))
             if len(rates) > 1 else 0)

    # Run both tests
    iq = iqr_test(rates, names)
    gr = grubbs_test_iterative(rates, names)
    iqs = {f["sample"] for f in iq}
    grs = {f["sample"] for f in gr}
    confirmed = iqs & grs  # BOTH tests must agree

    # Build per-sample report
    report = []
    for s, r, lib in zip(names, rates, libs):
        sd2 = next(x for x in samples if x["sample"] == s)
        report.append({
            "sample":            s,
            "pct_unique":        round(r, 2),
            "total_reads_M":     round(lib, 2),
            "pct_multi":         round(sd2.get("pct_multi", 0), 2),
            "pct_too_short":     round(sd2.get("pct_too_short", 0), 2),
            "iqr_flag":          s in iqs,
            "grubbs_flag":       s in grs,
            "confirmed_outlier": s in confirmed,
            "iqr_details":       next((f for f in iq if f["sample"] == s), None),
            "grubbs_details":    next((f for f in gr if f["sample"] == s), None),
        })

    flagged = [s["sample"] for s in report if s["confirmed_outlier"]]
    result = {
        "layer":        "alignment",
        "aligner":      aligner,
        "status":       "outliers_found" if flagged else "all_pass",
        "n_samples":    len(samples),
        "skipped":      skipped,
        "mean_mapping": round(mu, 2),
        "std_mapping":  round(std_r, 2),
        "samples":      report,
        "flagged":      flagged,
    }
    (od / "outlier_flag.json").write_text(json.dumps(result, indent=2))
    print(f"Alignment QC: {len(samples)} samples, mean={mu:.1f}%, flagged={len(flagged)}")
    if skipped:
        print(f"  Skipped (unparseable): {', '.join(skipped)}")
    for s in report:
        if s["confirmed_outlier"]:
            print(f"  OUTLIER: {s['sample']}  {s['pct_unique']}%  "
                  f"(IQR + Grubbs both agree)")
    print(f"Report: {od}/outlier_flag.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Alignment QC outlier detection (Layer 1)")
    ap.add_argument("--star_dir", required=True,
                    help="Directory containing alignment log files")
    ap.add_argument("--outdir", required=True,
                    help="Output directory for outlier_flag.json")
    ap.add_argument("--aligner", default="star", choices=["star", "hisat2"],
                    help="Aligner type (default: star)")
    a = ap.parse_args()
    main(a.star_dir, a.outdir, a.aligner)
