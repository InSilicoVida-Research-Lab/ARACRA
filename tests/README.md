# ARACRA — tests

Five tiers, ordered by cost. Nothing here installs anything, touches
`~/databases`, or writes outside `tests/.out/`.

```bash
bash tests/run_tests.sh              # tiers 1-2   ~5 seconds
bash tests/run_tests.sh --solve      # + tier 3    ~2-5 minutes
bash tests/run_tests.sh --fixture    # + tier 4    ~5 minutes
bash tests/run_tests.sh --pathway    # + tier 5    ~1-2 minutes
bash tests/run_tests.sh --all        # everything
bash tests/run_tests.sh --complexity # report only
```

**A `skip` is never a `pass`.** If a check could not run on this machine — no
Rscript, no conda, no network — it is reported as skipped and counted
separately. The suite is meant to be honest about what was actually verified,
not to look green.

---

## Tier 1 — `01_static.sh` (seconds, no dependencies)

Syntax, cross-file agreement, and regression guards against each defect
previously found in `main`: the stray `=1.10` redirect file, the doubled quote
in `nextflow.config`, the literal `~` that Groovy never expands, the hardcoded
`rtracklayer_1.66.0` tarball, the 8501/8502 port split, and the `.env` key that
`run_app.sh` wrote but the app needed.

The `.env` parity check is the useful one to keep: it parses the keys
`aracra_star_app.py` reads and the keys `lib/aracra_common.sh` writes, and fails
if they diverge. That class of bug is invisible until a user hits it.

**This is the tier to put in CI.** It needs nothing but bash and Python.

## Tier 2 — `02_preflight.sh` (seconds, no downloads)

Run this *before* `setup.sh`. Reports cores, RAM, disk, architecture and which
aligner will be chosen, then estimates what the install will cost in GB and
minutes. Finds in five seconds what would otherwise surface ninety minutes into
a failed install — most often insufficient disk, or arm64 hardware that needs
`dgx_install_branch` rather than `main`.

## Tier 3 — `03_solve.sh` (minutes, network, no install)

The most valuable check, because it tests the one thing about `setup.sh` that
cannot be verified by reading it: whether the version-constrained solve actually
resolves. Uses `--dry-run`, so it downloads no packages.

It solves the pinned spec *and* the unpinned fallback. If the pinned one fails
you immediately know whether the version floors are the cause or the channels
are simply broken — and the floors live in one file,
`ARACRA/lib/aracra_common.sh`, so the fix is one line.

Requires conda or mamba on PATH; skips cleanly otherwise.

## Tier 4 — `04_fixture.sh` (minutes, needs the ARACRA env)

End-to-end statistics on a synthetic dataset, entering through **Direct Mode**
so no genome, index or aligner is needed. Runs the real `run_deseq2.R` and
`run_dromics.R`, not mocks.

```bash
conda activate ~/miniforge3/envs/test_ARACRA
bash tests/04_fixture.sh
```

`fixtures/make_fixture.py` builds 500 genes × 27 samples (8 half-log doses plus
vehicle, 3 replicates) from a **fixed seed**, so the count matrix is identical on
every machine. A change in output therefore means a change in the code, not
sampling noise.

Three gene classes are planted on purpose:

| class | n | purpose |
| ----- | - | ------- |
| `null` | 420 | no dose dependence — measures false positives |
| `responder` | 60 | Hill curves with BMDs inside the tested range |
| `low_mover` | 20 | responses steep enough to land **below** the lowest dose |

The `low_mover` genes exist specifically to trigger the NTP extrapolation flag.
The script runs DRomics twice, with `--bmd_extrap_filter FALSE` and `TRUE`, and
asserts that:

- with the filter **off**, `n_removed_extrap == 0` — the default must not have
  changed, because that is what reproduces published ARACRA runs;
- with the filter **on**, every flagged gene is actually removed;
- turning it on **raises** the rank25 tPOD rather than lowering it.

That third assertion is the regression test for the defect originally found on
BPA, where the extrapolation filter flagged but never removed.

## Tier 5 — `05_pathway_fixture.sh` (minutes, needs org.Hs.eg.db/GO.db)

The main fixture's gene IDs (`ENSG_RESP_0000` etc.) never match real GO/KEGG/
MSigDB annotations, so `run_dromics.R`'s pathway-level tPOD code — the part
that maps genes to gene sets, takes each set's median BMD, and reports the
single lowest one as the tPOD (NTP 2018's own stated method — verified
against [NTP Research Report 5](https://www.ncbi.nlm.nih.gov/books/NBK531562/))
— had never actually run in this suite. This tier uses real Ensembl IDs
(`fixtures/dump_real_pathway.R` picks one real, well-sized GO:BP term plus a
background sample) so it finally does, on two scenarios:

- **null** — no true signal anywhere. Informational: reports whatever tPOD
  (if any) the method turns up from pure chance, given hundreds of real gene
  sets get tested with no multiple-comparisons correction — the noise floor
  a real result should be judged against.
- **spike** — every gene in the chosen pathway gets an identical known-EC50
  response. Asserts the pathway is detected (present in the ranked output,
  high coverage, BMD near the true EC50) — but *not* that it wins the
  single-minimum contest, because a first real run already showed a 3-gene
  overlapping subset of the same spiked genes out-noise it and win instead.
  That is reproduced NTP-2018-faithful behavior (no correction for gene-set
  overlap), not a code defect — see the printed ranked table for a concrete
  example.

## Complexity — `complexity_report.py`

```bash
python3 tests/complexity_report.py          # human readable
python3 tests/complexity_report.py --json   # machine readable
```

Not a quality score. It answers: which files would a new maintainer have to
read, how many branch points does a run cross, how many knobs exist and are they
wired end to end, and what does the install actually cost.

Two findings worth acting on, as of this writing:

- **70% of the code sits in three files**, and `aracra_star_app.py` alone holds
  ~2,600 lines with ~500 branch points — far more than `main.nf` or
  `run_dromics.R`. The GUI, not the science, is the largest maintenance surface,
  and it is the least covered by these tests.
- **A full install is ~71 GB**; the leanest viable configuration
  (`--skip-index --skip-screen`) is ~26 GB. Worth stating in the manual, since
  the README's stated 100 GB minimum is the right number but the reason for it
  isn't obvious.

## Adding a check

Source `lib/assert.sh` and use `pass` / `fail` / `skip`, or the helpers
`assert_ok`, `assert_contains`, `assert_absent`, `need_tool`. Call `summary` at
the end; it sets the exit status. No framework, no install.
