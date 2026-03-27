#!/usr/bin/env python3
# =============================================================================
#  temposeq_counts.py — TempO-Seq probe → gene count aggregation
#
#  Takes samtools idxstats output + probe-to-gene mapping and produces
#  a gene-level count matrix compatible with DRomics/DESeq2.
#
#  Key behaviors:
#    - Multiple probes per gene → counts are SUMMED
#    - Gene IDs: ENSEMBL preferred, symbol fallback, consistent per gene
#    - Probes not in mapping are excluded
#    - Output format matches featureCounts: first column = gene ID, rest = samples
#
#  Usage:
#    python3 temposeq_counts.py \
#        --mapping probe_to_gene.tsv \
#        --counts_dir . \
#        --outfile gene_counts.csv
# =============================================================================

import csv, argparse, sys, json
from pathlib import Path
from collections import defaultdict


def load_probe_mapping(mapping_path):
    """Load probe-to-gene mapping with consistent gene ID resolution.
    Returns:
       probe2gene: {probe_name: canonical_gene_id}
       gene_info:  {canonical_gene_id: {symbol, ensembl, probes: [...]}}
    """
    raw_probes = []
    with open(mapping_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            raw_probes.append({
                'probe': row['probe_name'].strip(),
                'symbol': row.get('gene_symbol', '').strip(),
                'ensembl': row.get('ensembl_id', '').strip(),
                'entrez': row.get('entrez_id', '').strip(),
            })

    # Build symbol↔ensembl cross-reference from all probes
    symbol_to_ensembl = {}
    ensembl_to_symbol = {}
    for p in raw_probes:
        if p['ensembl'] and p['symbol']:
            symbol_to_ensembl.setdefault(p['symbol'], p['ensembl'])
            ensembl_to_symbol.setdefault(p['ensembl'], p['symbol'])

    # Assign canonical gene ID per probe
    probe2gene = {}
    gene_info = {}

    for p in raw_probes:
        if not p['probe']:
            continue

        ensembl = p['ensembl']
        symbol = p['symbol']

        # Canonical ID: ENSEMBL preferred
        if ensembl:
            gene_id = ensembl
        elif symbol and symbol in symbol_to_ensembl:
            gene_id = symbol_to_ensembl[symbol]
        elif symbol:
            gene_id = symbol
        else:
            continue  # no annotation → skip

        probe2gene[p['probe']] = gene_id

        if gene_id not in gene_info:
            gene_info[gene_id] = {
                'symbol': ensembl_to_symbol.get(gene_id, symbol),
                'ensembl': ensembl or symbol_to_ensembl.get(symbol, ''),
                'probes': [],
            }
        gene_info[gene_id]['probes'].append(p['probe'])

    n_multi = sum(1 for g in gene_info.values() if len(g['probes']) > 1)
    n_ensg = sum(1 for gid in gene_info if gid.startswith('ENSG'))

    print(f"Probe-to-gene mapping:")
    print(f"  Probes mapped:         {len(probe2gene)}")
    print(f"  Unique genes:          {len(gene_info)}")
    print(f"  With ENSEMBL ID:       {n_ensg}")
    print(f"  Symbol-only:           {len(gene_info) - n_ensg}")
    print(f"  Multi-probe genes:     {n_multi}")

    if n_multi > 0:
        examples = [(gid, info) for gid, info in gene_info.items()
                     if len(info['probes']) > 1][:5]
        for gid, info in examples:
            print(f"    {gid} ({info['symbol']}): "
                  f"{len(info['probes'])} probes → {', '.join(info['probes'][:4])}")

    return probe2gene, gene_info


def parse_idxstats(filepath):
    """Parse samtools idxstats. Returns {probe_name: mapped_count}"""
    counts = {}
    with open(filepath) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            if parts[0] == '*':
                continue
            counts[parts[0]] = int(parts[2])
    return counts


def find_count_files(counts_dir):
    """Find idxstats files in directory."""
    counts_dir = Path(counts_dir)
    for pattern in ["*_idxstats.tsv", "*.idxstats", "*_idxstats.*",
                     "*_probe_counts.tsv", "**/*_idxstats.tsv"]:
        files = sorted(counts_dir.glob(pattern))
        if files:
            return files

    print(f"ERROR: No count files found in {counts_dir}", file=sys.stderr)
    print(f"  Contents:", file=sys.stderr)
    for f in sorted(counts_dir.iterdir()):
        print(f"    {f.name}", file=sys.stderr)
    return []


def aggregate_counts(mapping_path, counts_dir, outfile):
    probe2gene, gene_info = load_probe_mapping(mapping_path)

    stats_files = find_count_files(counts_dir)
    if not stats_files:
        sys.exit(1)

    print(f"\nFound {len(stats_files)} sample file(s)")

    sample_names = []
    gene_counts = defaultdict(lambda: defaultdict(int))
    sample_qc = {}

    for sf in stats_files:
        sample = sf.stem
        for sfx in ['_idxstats', '_probe_counts', '.idxstats']:
            sample = sample.replace(sfx, '')
        sample_names.append(sample)

        probe_counts = parse_idxstats(sf)
        total_reads = sum(probe_counts.values())
        assigned_reads = 0
        probes_hit = sum(1 for c in probe_counts.values() if c > 0)

        for probe_name, count in probe_counts.items():
            if probe_name in probe2gene:
                gene_id = probe2gene[probe_name]
                gene_counts[gene_id][sample] += count
                assigned_reads += count

        pct = round(100 * assigned_reads / max(1, total_reads), 1)
        sample_qc[sample] = {
            'total_reads': total_reads,
            'assigned_reads': assigned_reads,
            'assignment_pct': pct,
            'probes_hit': probes_hit,
        }
        print(f"  {sample}: {total_reads:,} reads → {assigned_reads:,} assigned "
              f"({pct}%) | {probes_hit} probes detected")

    # Write count matrix
    all_genes = sorted(gene_counts.keys())

    with open(outfile, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + sample_names)
        for gene_id in all_genes:
            row = [gene_id] + [gene_counts[gene_id].get(s, 0) for s in sample_names]
            writer.writerow(row)

    # Validation
    n_with_reads = sum(1 for g in all_genes
                       if any(gene_counts[g].get(s, 0) > 0 for s in sample_names))
    n_ensg = sum(1 for g in all_genes if g.startswith('ENSG'))

    print(f"\n{'='*50}")
    print(f"Gene count matrix: {len(all_genes)} genes × {len(sample_names)} samples")
    print(f"  Expressed genes:   {n_with_reads}")
    print(f"  ENSEMBL IDs:       {n_ensg} ({round(100*n_ensg/max(1,len(all_genes)))}%)")
    print(f"  Symbol-only:       {len(all_genes) - n_ensg}")

    if len(all_genes) - n_ensg > 0:
        syms = [g for g in all_genes if not g.startswith('ENSG')][:10]
        print(f"  Symbol examples:   {', '.join(syms)}")

    print(f"Output: {outfile}")

    # Summary JSON
    Path(outfile).parent.joinpath("temposeq_count_summary.json").write_text(
        json.dumps({
            'n_genes': len(all_genes),
            'n_expressed': n_with_reads,
            'n_samples': len(sample_names),
            'n_ensembl': n_ensg,
            'n_symbol_only': len(all_genes) - n_ensg,
            'multi_probe_genes': sum(1 for g in gene_info.values() if len(g['probes']) > 1),
            'sample_qc': sample_qc,
            'sample_names': sample_names,
        }, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Aggregate TempO-Seq probe counts to gene-level matrix")
    ap.add_argument("--mapping", required=True,
                    help="probe_to_gene.tsv from parse_temposeq_manifest.py")
    ap.add_argument("--counts_dir", required=True,
                    help="Directory containing *_idxstats.tsv files")
    ap.add_argument("--outfile", default="gene_counts.csv",
                    help="Output gene-level count matrix")
    args = ap.parse_args()
    aggregate_counts(args.mapping, args.counts_dir, args.outfile)
