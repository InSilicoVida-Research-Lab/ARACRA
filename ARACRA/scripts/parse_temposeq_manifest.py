#!/usr/bin/env python3
# =============================================================================
#  parse_temposeq_manifest.py — TempO-Seq manifest parser
#  Converts BioSpyder manifest CSV to:
#    1. probes.fa       — FASTA file of probe sequences (for STAR index)
#    2. probe_to_gene.tsv — mapping of probe_name → gene_symbol, ensembl_id
#
#  Usage:
#    python3 parse_temposeq_manifest.py \
#        --manifest manifest.csv \
#        --outdir .
# =============================================================================

import csv, argparse, sys
from pathlib import Path


def parse_manifest(manifest_path, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fa_path   = outdir / "probes.fa"
    map_path  = outdir / "probe_to_gene.tsv"

    probes = []
    seen_names = set()

    with open(manifest_path, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)

        # Normalize column names (strip whitespace, uppercase)
        if reader.fieldnames:
            col_map = {c.strip().upper(): c for c in reader.fieldnames}
        else:
            print("ERROR: Empty manifest file", file=sys.stderr)
            sys.exit(1)

        # Find required columns
        def find_col(candidates):
            for c in candidates:
                if c in col_map:
                    return col_map[c]
            return None

        col_probe_name = find_col(["PROBE_NAME", "PROBE_ID"])
        col_gene       = find_col(["GENE_SYMBOL", "GENE_NAME", "GENE"])
        col_sequence   = find_col(["PROBE_SEQUENCE", "SEQUENCE", "SEQ"])
        col_ensembl    = find_col(["ENSEMBL_GENE_ID", "ENSEMBL_ID", "ENSG"])
        col_entrez     = find_col(["ENTREZ_ID", "ENTREZ", "GENE_ID"])

        if not col_probe_name or not col_sequence:
            print(f"ERROR: Cannot find PROBE_NAME and PROBE_SEQUENCE columns",
                  file=sys.stderr)
            print(f"  Available columns: {list(reader.fieldnames)}", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            name = row.get(col_probe_name, "").strip()
            seq  = row.get(col_sequence, "").strip().upper()
            gene = row.get(col_gene, "").strip() if col_gene else ""
            ensg = row.get(col_ensembl, "").strip() if col_ensembl else ""
            entrez = row.get(col_entrez, "").strip() if col_entrez else ""

            if not name or not seq:
                continue
            # Skip if sequence contains non-ATCGN characters
            if not all(c in "ATCGN" for c in seq):
                print(f"  WARNING: Skipping {name} — invalid sequence characters",
                      file=sys.stderr)
                continue
            # Handle duplicate probe names
            if name in seen_names:
                continue
            seen_names.add(name)

            probes.append({
                "probe_name": name,
                "gene_symbol": gene,
                "ensembl_id": ensg,
                "entrez_id": entrez,
                "sequence": seq,
            })

    if not probes:
        print("ERROR: No valid probes found in manifest", file=sys.stderr)
        sys.exit(1)

    # Write FASTA
    with open(fa_path, 'w') as fa:
        for p in probes:
            fa.write(f">{p['probe_name']}\n{p['sequence']}\n")

    # Write probe-to-gene mapping
    with open(map_path, 'w') as mp:
        mp.write("probe_name\tgene_symbol\tensembl_id\tentrez_id\n")
        for p in probes:
            mp.write(f"{p['probe_name']}\t{p['gene_symbol']}\t"
                     f"{p['ensembl_id']}\t{p['entrez_id']}\n")

    # Summary stats
    genes = set(p['gene_symbol'] for p in probes if p['gene_symbol'])
    multi_probe_genes = {}
    for p in probes:
        g = p['gene_symbol']
        if g:
            multi_probe_genes[g] = multi_probe_genes.get(g, 0) + 1
    n_multi = sum(1 for v in multi_probe_genes.values() if v > 1)

    print(f"TempO-Seq manifest parsed:")
    print(f"  Total probes:       {len(probes)}")
    print(f"  Unique genes:       {len(genes)}")
    print(f"  Multi-probe genes:  {n_multi}")
    print(f"  Probe length range: {min(len(p['sequence']) for p in probes)}"
          f"-{max(len(p['sequence']) for p in probes)} bp")
    print(f"  Output FASTA:       {fa_path}")
    print(f"  Output mapping:     {map_path}")

    # Also write a summary JSON
    import json
    summary = {
        "total_probes": len(probes),
        "unique_genes": len(genes),
        "multi_probe_genes": n_multi,
        "min_probe_length": min(len(p['sequence']) for p in probes),
        "max_probe_length": max(len(p['sequence']) for p in probes),
    }
    (outdir / "manifest_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Parse TempO-Seq manifest to FASTA + mapping")
    ap.add_argument("--manifest", required=True, help="BioSpyder manifest CSV file")
    ap.add_argument("--outdir", default=".", help="Output directory")
    args = ap.parse_args()
    parse_manifest(args.manifest, args.outdir)
