#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def _totalCores    = Runtime.runtime.availableProcessors()

// ── Direct analysis mode: skip straight to DESeq2/DRomics ──
params.count_matrix     = null
params.run_name         = null
params.treatment        = "Bisphenol A"
params.control          = "Control"
params.run_deg          = true
params.run_dromics      = true
params.fdr_strict       = 0.01
params.fdr_relaxed      = 0.05
params.log2fc_threshold = 1.0
params.read_threshold   = 1000000
params.dr_fdr           = 0.05
params.dr_criterion     = "AIC"
params.perform_bmd      = true
params.bmd_bootstrap    = false

// ── Aligner selection: "star" or "hisat2" ──
params.aligner          = "star"
params.hisat2_index     = "${HOME}/databases/hg38_reference/hisat2_index_hg38/genome_tran"
params.run_hisat2       = false

def direct_mode = (params.count_matrix != null)
def run_base = params.run_name ? "${params.outdir}/runs/${params.run_name}" : params.outdir
def use_star   = (params.aligner == "star")
def use_hisat2 = (params.aligner == "hisat2")

log.info """
╔══════════════════════════════════════════════════════════╗
║         RNA-seq / TempO-Seq Analysis Pipeline           ║
╠══════════════════════════════════════════════════════════╣
║  Mode            : ${direct_mode ? '⚡ DIRECT ANALYSIS' : '🔬 FULL PIPELINE'}
║  Aligner         : ${params.aligner.toUpperCase()}
║  Layout          : ${params.layout}
║  Metadata        : ${params.metadata}
║  Output          : ${params.outdir}
${params.run_name ? '║  Run name        : ' + params.run_name : ''}
${direct_mode ? '║  Count matrix    : ' + params.count_matrix : ''}
╠══════════════════════════════════════════════════════════╣
║  Core strategy    : ${_totalCores} total cores
║    Download       : 1 core × 3 parallel (network-bound)
║    Fastp          : 1 core × ${Math.max(1, _totalCores - 1)} parallel
║    FastQ Screen   : 2 cores (independent)
║    ${params.aligner.toUpperCase().padRight(5)}         : ${Math.max(2, (int)(_totalCores / 2))} cores × 2 parallel${use_star ? ' (shared memory)' : ''}
║    Post-align QC  : 1 core × ${Math.max(1, _totalCores - 2)} parallel
║    Final steps    : ${_totalCores} cores
╠══════════════════════════════════════════════════════════╣
║  Steps enabled:
║   1.  Download SRA       : ${params.run_download}
║   2.  fastp (trim+QC)    : ${params.run_fastp}
║   3.  FastQ Screen       : ${params.run_fastq_screen}
║   4a. STAR alignment     : ${use_star}
║   4b. HISAT2 alignment   : ${use_hisat2}
║   5.  samtools stats     : ${params.run_samtools_stats}
║   6.  Strandedness       : ${params.run_strandedness}
║   7a. Coverage (HK)      : ${params.run_coverage_hk}
║   7b. Coverage (full)    : ${params.run_coverage_full}
║   8.  Read distribution  : ${params.run_read_distribution}
║   9.  Picard RNA metrics : ${params.run_picard}
║  10.  featureCounts      : ${params.run_featurecounts}
║  11.  Salmon             : ${params.run_salmon}
║  12.  Merge counts       : ${params.run_merge_counts}
║  13.  MultiQC            : ${params.run_multiqc}
║  14.  DESeq2             : ${params.run_deg}
║  15.  DRomics/BMD        : ${params.run_dromics}
╚══════════════════════════════════════════════════════════╝
"""

// ============================================================
//  PROCESSES
// ============================================================

process PARSE_METADATA {
    label 'minimal'
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:  path metadata
    output:
    path "sample_ids.txt", emit: sample_ids
    path "metadata.csv",   emit: metadata_csv

    script:
    """
    #!/usr/bin/env python3
    import openpyxl, csv
    wb = openpyxl.load_workbook("${metadata}", read_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    header = [str(h).strip() for h in rows[0]]
    sample_col = 0
    for i, h in enumerate(header):
        if 'sample' in h.lower() or 'srr' in h.lower():
            sample_col = i; break
    with open("sample_ids.txt","w") as fi, open("metadata.csv","w",newline="") as fc:
        w = csv.writer(fc); w.writerow(header)
        for row in rows[1:]:
            if row[sample_col] is not None:
                sid = str(row[sample_col]).strip()
                if sid:
                    fi.write(sid+"\\n")
                    w.writerow([str(v).strip() if v else "" for v in row])
    """
}

// ---- LIGHT: 1 core, many parallel ----

process DOWNLOAD_SRA {
    label 'download'
    tag "${sample_id}"
    publishDir "${params.outdir}/data", mode: 'copy'

    input:  val sample_id
    output: tuple val(sample_id), path("${sample_id}*.fastq.gz"), emit: fastq

    script:
    if (params.layout == "PE")
        """
        prefetch ${sample_id} --max-size 50G
        fasterq-dump --split-3 ${sample_id} --threads ${task.cpus} --progress
        pigz -p ${task.cpus} ${sample_id}_1.fastq ${sample_id}_2.fastq 2>/dev/null || \
            gzip ${sample_id}_1.fastq ${sample_id}_2.fastq
        rm -rf ${sample_id}/ ${sample_id}.fastq
        """
    else
        """
        prefetch ${sample_id} --max-size 50G
        fasterq-dump --split-3 ${sample_id} --threads ${task.cpus} --progress
        [ -f ${sample_id}_1.fastq ] && [ ! -f ${sample_id}_2.fastq ] && \
            mv ${sample_id}_1.fastq ${sample_id}.fastq
        pigz -p ${task.cpus} ${sample_id}.fastq 2>/dev/null || gzip ${sample_id}.fastq
        rm -rf ${sample_id}/
        """
}

process FASTP {
    label 'pre_qc'
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: "*.{html,json}"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*_trimmed*"

    input:  tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("*_trimmed*.fastq.gz"), emit: trimmed_reads
    path("${sample_id}_fastp.html"),                    emit: html
    path("${sample_id}_fastp.json"),                    emit: json

    script:
    if (params.layout == "PE")
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${sample_id}_trimmed_1.fastq.gz -O ${sample_id}_trimmed_2.fastq.gz \
            -h ${sample_id}_fastp.html -j ${sample_id}_fastp.json --thread ${task.cpus}
        """
    else
        """
        fastp -i ${reads[0]} -o ${sample_id}_trimmed.fastq.gz \
            -h ${sample_id}_fastp.html -j ${sample_id}_fastp.json --thread ${task.cpus}
        """
}

// ---- SCREEN: 2 cores, independent ----

process FASTQ_SCREEN {
    label 'screen'
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/fastq_screen", mode: 'copy'

    input:  tuple val(sample_id), path(reads)
    output:
    path("*.txt"),  emit: report
    path("*.html"), emit: html
    path("*.png"),  optional: true, emit: png

    script:
    """
    fastq_screen --conf ${params.fastq_screen_conf} --aligner bowtie2 \
        --threads ${task.cpus} --outdir . ${reads[0]}
    """
}

// ---- HEAVY: All cores, 1 at a time ----

// ---- STAR: Load genome into shared memory (once, before all alignments) ----
process STAR_LOAD_GENOME {
    label 'final'
    
    output: val(true), emit: genome_loaded

    script:
    """
    # Safety: remove any orphaned shared memory from a previous crashed run
    STAR --genomeDir ${params.star_index} --genomeLoad Remove 2>/dev/null || true
    # Now load fresh
    STAR --genomeDir ${params.star_index} --genomeLoad LoadAndExit
    """
}

// ---- STAR: Remove genome from shared memory (after all alignments) ----
process STAR_REMOVE_GENOME {
    label 'minimal'

    input: val(all_done)

    script:
    """
    STAR --genomeDir ${params.star_index} --genomeLoad Remove
    """
}

process STAR_ALIGN {
    label 'heavy'
    tag "${sample_id}"
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val genome_ready    // ensures genome is loaded first

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), emit: bam
    path("${sample_id}_Log.final.out"),    emit: log_final
    path("${sample_id}_Log.out"),          emit: log_out
    path("${sample_id}_ReadsPerGene.out.tab"), emit: reads_per_gene
    path("${sample_id}_SJ.out.tab"),       emit: splice_junctions

    script:
    def read_files = (params.layout == "PE") ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    def read_cmd   = reads[0].name.endsWith('.gz') ? "--readFilesCommand zcat" : ""
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${params.star_index} \
        --genomeLoad LoadAndKeep \
        --readFilesIn ${read_files} ${read_cmd} \
        --outFileNamePrefix ${sample_id}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --quantMode GeneCounts \
        --limitBAMsortRAM ${params.star_ram}
    """
}

// ---- LIGHT: Post-alignment QC (1 core, many parallel) ----

// ---- HISAT2 Alignment (alternative to STAR, uses ~8 GB RAM) ----
process HISAT2_ALIGN {
    label 'heavy'
    tag "${sample_id}"
    publishDir "${params.outdir}/hisat2", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam
    path("${sample_id}_align_summary.txt"),                 emit: log_summary

    script:
    def read_input = (params.layout == "PE") ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"
    def read_cmd   = reads[0].name.endsWith('.gz') ? "" : ""  // hisat2 handles .gz natively
    def sort_mb    = Math.max(200, (int)((task.memory.toMega() * 0.3) / task.cpus))
    """
    hisat2 -x ${params.hisat2_index} \
        ${read_input} \
        -p ${task.cpus} --dta --no-unal \
        --summary-file ${sample_id}_align_summary.txt \
    | samtools sort -@ ${task.cpus} -m ${sort_mb}M \
        -T ${sample_id}_tmp -o ${sample_id}_sorted.bam
    """
}

process SAMTOOLS_INDEX {
    label 'post_align'
    tag "${sample_id}"
    publishDir "${params.outdir}/${params.aligner}", mode: 'copy'

    input:  tuple val(sample_id), path(bam)
    output: tuple val(sample_id), path(bam), path("${bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index ${bam}
    """
}

process SAMTOOLS_FLAGSTAT {
    label 'post_align'
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/samtools", mode: 'copy'

    input:  tuple val(sample_id), path(bam), path(bai)
    output: path("${sample_id}_flagstat.txt"), emit: flagstat

    script:
    """
    samtools flagstat ${bam} > ${sample_id}_flagstat.txt
    """
}

process SAMTOOLS_IDXSTATS {
    label 'post_align'
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/samtools", mode: 'copy'

    input:  tuple val(sample_id), path(bam), path(bai)
    output: path("${sample_id}_idxstats.txt"), emit: idxstats

    script:
    """
    samtools idxstats ${bam} > ${sample_id}_idxstats.txt
    """
}

process INFER_STRANDEDNESS {
    label 'post_align'
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/strandedness", mode: 'copy'

    input:  tuple val(sample_id), path(bam), path(bai)
    output:
    tuple val(sample_id), path("${sample_id}_strandedness.txt"), emit: strand_report
    tuple val(sample_id), env(STRAND_FEATURECOUNTS),             emit: strand_fc
    tuple val(sample_id), env(STRAND_SALMON),                    emit: strand_salmon

    script:
    """
    infer_experiment.py -r ${params.housekeeping_bed} -i ${bam} > ${sample_id}_strandedness.txt

    FRAC_REVERSE=\$(grep '+-,-+' ${sample_id}_strandedness.txt | awk '{print \$NF}')
    FRAC_FORWARD=\$(grep '++,--' ${sample_id}_strandedness.txt | awk '{print \$NF}')

    if (( \$(echo "\$FRAC_REVERSE > 0.75" | bc -l) )); then
        STRAND_FEATURECOUNTS="2"
        STRAND_SALMON="${params.layout == 'PE' ? 'ISR' : 'SR'}"
    elif (( \$(echo "\$FRAC_FORWARD > 0.75" | bc -l) )); then
        STRAND_FEATURECOUNTS="1"
        STRAND_SALMON="${params.layout == 'PE' ? 'ISF' : 'SF'}"
    else
        STRAND_FEATURECOUNTS="0"
        STRAND_SALMON="IU"
    fi

    echo "featureCounts_strand=\$STRAND_FEATURECOUNTS" >> ${sample_id}_strandedness.txt
    echo "salmon_libtype=\$STRAND_SALMON" >> ${sample_id}_strandedness.txt
    """
}

process GENEBODY_COVERAGE_HK {
    label 'post_align'
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/genebody_coverage/housekeeping", mode: 'copy'

    input:  tuple val(sample_id), path(bam), path(bai)
    output: path("${sample_id}_hk_coverage*"), emit: coverage

    script:
    """
    geneBody_coverage.py -r ${params.housekeeping_bed} -i ${bam} -o ${sample_id}_hk_coverage
    """
}

process GENEBODY_COVERAGE_FULL {
    label 'post_align'
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/genebody_coverage/full", mode: 'copy'

    input:  tuple val(sample_id), path(bam), path(bai)
    output: path("${sample_id}_full_coverage*"), emit: coverage

    script:
    """
    geneBody_coverage.py -r ${params.bed} -i ${bam} -o ${sample_id}_full_coverage
    """
}

process READ_DISTRIBUTION {
    label 'post_align'
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/read_distribution", mode: 'copy'

    input:  tuple val(sample_id), path(bam), path(bai)
    output: path("${sample_id}_read_dist.txt"), emit: read_dist

    script:
    """
    read_distribution.py -r ${params.bed} -i ${bam} > ${sample_id}_read_dist.txt
    """
}

process PICARD_RNA_METRICS {
    label 'post_align'
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/picard", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val strand

    output:
    path("${sample_id}_rnaseq_metrics.txt"), emit: metrics
    path("${sample_id}_coverage.pdf"),       emit: coverage_pdf

    script:
    def ps = (strand == "2") ? "SECOND_READ_TRANSCRIPTION_STRAND" :
             (strand == "1") ? "FIRST_READ_TRANSCRIPTION_STRAND" : "NONE"
    """
    picard CollectRnaSeqMetrics I=${bam} O=${sample_id}_rnaseq_metrics.txt \
        REF_FLAT=${params.refflat} STRAND=${ps} CHART_OUTPUT=${sample_id}_coverage.pdf
    """
}

// ---- FINAL: All cores ----

process FEATURECOUNTS {
    label 'final'
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    path bam_files
    path bai_files
    val strand

    output:
    path("gene_counts.txt"),         emit: counts
    path("gene_counts.txt.summary"), emit: summary

    script:
    def pe_flag = (params.layout == "PE") ? "-p --countReadPairs" : ""
    """
    featureCounts -a ${params.gtf} -o gene_counts.txt \
        -T ${task.cpus} -s ${strand} ${pe_flag} *.bam
    """
}

process SALMON_QUANT {
    label 'final'
    tag "${sample_id}"
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val libtype

    output: path("${sample_id}_quant"), emit: quant_dir

    script:
    def ri = (params.layout == "PE") ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${reads[0]}"
    """
    salmon quant -i ${params.salmon_index} -l ${libtype} ${ri} \
        -o ${sample_id}_quant --threads ${task.cpus} --validateMappings
    """
}

process MAKE_COUNT_MATRIX {
    label 'minimal'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path counts_file
    path metadata_csv

    output:
    path("gene_count_matrix.csv"), emit: count_matrix
    path("sample_metadata.csv"),   emit: sample_info

    script:
    """
    #!/usr/bin/env python3
    import csv, re, shutil
    with open("${counts_file}") as fh, open("gene_count_matrix.csv","w",newline="") as out:
        w = csv.writer(out)
        for line in fh:
            if line.startswith("#"): continue
            parts = line.strip().split("\\t")
            if parts[0] == "Geneid":
                samples = [re.sub(r'_Aligned.*','',c.split("/")[-1]) for c in parts[6:]]
                w.writerow(["gene_id"] + samples)
            else:
                w.writerow([parts[0]] + parts[6:])
    shutil.copy("${metadata_csv}","sample_metadata.csv")
    print(f"Count matrix: {len(samples)} samples")
    """
}

process MERGE_SALMON {
    label 'minimal'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:  path quant_dirs
    output:
    path("salmon_tpm_matrix.csv"),    emit: tpm_matrix
    path("salmon_counts_matrix.csv"), emit: counts_matrix

    script:
    """
    #!/usr/bin/env python3
    import os, csv, glob
    qfiles = sorted(glob.glob("*_quant/quant.sf"))
    samples = [os.path.basename(os.path.dirname(f)).replace("_quant","") for f in qfiles]
    tpm, cts = {}, {}
    for qf, s in zip(qfiles, samples):
        with open(qf) as fh:
            for row in csv.DictReader(fh, delimiter="\\t"):
                tx = row["Name"]
                if tx not in tpm: tpm[tx], cts[tx] = {}, {}
                tpm[tx][s] = row["TPM"]; cts[tx][s] = row["NumReads"]
    for fname, data in [("salmon_tpm_matrix.csv",tpm),("salmon_counts_matrix.csv",cts)]:
        with open(fname,"w",newline="") as f:
            w = csv.writer(f); w.writerow(["transcript_id"]+samples)
            for tx in sorted(data): w.writerow([tx]+[data[tx].get(s,0) for s in samples])
    print(f"Salmon: {len(tpm)} transcripts × {len(samples)} samples")
    """
}

process MULTIQC {
    label 'minimal'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:  path('qc_files/*')
    output:
    path("multiqc_report.html"),  emit: report
    path("multiqc_report_data"),  emit: data

    script:
    """
    multiqc . --force -n multiqc_report.html -o .
    """
}

// ---- DESeq2 Analysis ----
process DESEQ2_ANALYSIS {
    label 'final'
    publishDir "${run_base}/deg_results", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path count_matrix
    path metadata

    output:
    path("*.csv"),  optional: true, emit: csvs
    path("*.png"),  optional: true, emit: plots
    path("*.json"), optional: true, emit: summary

    script:
    """
    Rscript ${projectDir}/scripts/run_deseq2.R \
        --counts      ${count_matrix} \
        --metadata    ${metadata} \
        --treatment   "${params.treatment}" \
        --control     "${params.control}" \
        --outdir      . \
        --fdr_strict  ${params.fdr_strict} \
        --fdr_relaxed ${params.fdr_relaxed} \
        --log2fc      ${params.log2fc_threshold} \
        --read_thresh ${params.read_threshold}
    """
}

// ---- DRomics / BMD Analysis ----
process DROMICS_ANALYSIS {
    label 'final'
    publishDir "${run_base}/dromics_results", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path count_matrix
    path metadata

    output:
    path("*.csv"),  optional: true, emit: csvs
    path("*.png"),  optional: true, emit: plots
    path("*.json"), optional: true, emit: summary

    script:
    """
    Rscript ${projectDir}/scripts/run_dromics.R \
        --counts    ${count_matrix} \
        --metadata  ${metadata} \
        --treatment "${params.treatment}" \
        --control   "${params.control}" \
        --outdir    . \
        --fdr       ${params.dr_fdr} \
        --criterion "${params.dr_criterion}" \
        --bmd       ${params.perform_bmd} \
        --bootstrap ${params.bmd_bootstrap}
    """
}

// ============================================================
//  WORKFLOW
// ============================================================
workflow {

    if (direct_mode) {
        // ══════════════════════════════════════════════════════════
        // DIRECT MODE: Skip everything, go straight to DESeq2/DRomics
        // ══════════════════════════════════════════════════════════
        log.info "⚡ Direct analysis mode — skipping SRA/FASTQ/alignment/QC"

        count_ch    = Channel.fromPath(params.count_matrix)
        metadata_ch = Channel.fromPath(params.metadata)

        if (params.run_deg == true || params.run_deg == "true") {
            DESEQ2_ANALYSIS(count_ch, metadata_ch)
        }
        if (params.run_dromics == true || params.run_dromics == "true") {
            DROMICS_ANALYSIS(
                Channel.fromPath(params.count_matrix),
                Channel.fromPath(params.metadata)
            )
        }

    } else {
        // ══════════════════════════════════════════════════════════
        // FULL PIPELINE MODE
        // ══════════════════════════════════════════════════════════

        // 0. Parse metadata
        PARSE_METADATA(Channel.fromPath(params.metadata))
        sample_ids_ch = PARSE_METADATA.out.sample_ids.splitText().map{it.trim()}.filter{it.length()>0}

        // 1. Download
        if (params.run_download) {
            DOWNLOAD_SRA(sample_ids_ch)
            raw_reads_ch = DOWNLOAD_SRA.out.fastq
        }

        // 2. Fastp
        if (params.run_fastp) {
            FASTP(raw_reads_ch)
            reads_ch = FASTP.out.trimmed_reads
        } else {
            reads_ch = raw_reads_ch
        }

        // 3. FastQ Screen (independent)
        if (params.run_fastq_screen) {
            FASTQ_SCREEN(reads_ch)
        }

        // 4. Alignment — STAR or HISAT2
        if (use_star) {
            STAR_LOAD_GENOME()
            STAR_ALIGN(reads_ch, STAR_LOAD_GENOME.out.genome_loaded)
            SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
            indexed_bam_ch = SAMTOOLS_INDEX.out.indexed_bam
            STAR_REMOVE_GENOME(STAR_ALIGN.out.bam.collect().map { true })
        } else if (use_hisat2) {
            HISAT2_ALIGN(reads_ch)
            SAMTOOLS_INDEX(HISAT2_ALIGN.out.bam)
            indexed_bam_ch = SAMTOOLS_INDEX.out.indexed_bam
        }

        // 5-9. Post-alignment QC
        if (params.run_samtools_stats) {
            SAMTOOLS_FLAGSTAT(indexed_bam_ch)
            SAMTOOLS_IDXSTATS(indexed_bam_ch)
        }

        if (params.run_strandedness) {
            INFER_STRANDEDNESS(indexed_bam_ch)
            strand_fc_ch     = INFER_STRANDEDNESS.out.strand_fc
            strand_salmon_ch = INFER_STRANDEDNESS.out.strand_salmon
        }

        if (params.run_coverage_hk)       { GENEBODY_COVERAGE_HK(indexed_bam_ch) }
        if (params.run_coverage_full)      { GENEBODY_COVERAGE_FULL(indexed_bam_ch) }
        if (params.run_read_distribution)  { READ_DISTRIBUTION(indexed_bam_ch) }

        if (params.run_picard) {
            bam_strand_ch = indexed_bam_ch.join(strand_fc_ch)
            PICARD_RNA_METRICS(
                bam_strand_ch.map { sid, bam, bai, strand -> tuple(sid, bam, bai) },
                bam_strand_ch.map { sid, bam, bai, strand -> strand }
            )
        }

        // 10. featureCounts — ALL BAMs in single call
        if (params.run_featurecounts) {
            all_bams_ch = indexed_bam_ch.map { sid, bam, bai -> bam }.collect()
            all_bais_ch = indexed_bam_ch.map { sid, bam, bai -> bai }.collect()
            strand_val  = strand_fc_ch.map { sid, strand -> strand }.first()
            FEATURECOUNTS(all_bams_ch, all_bais_ch, strand_val)
        }

        // 11. Salmon
        if (params.run_salmon) {
            salmon_in = reads_ch.join(strand_salmon_ch)
            SALMON_QUANT(
                salmon_in.map { sid, reads, lib -> tuple(sid, reads) },
                salmon_in.map { sid, reads, lib -> lib }
            )
        }

        // 12. Count matrix
        if (params.run_merge_counts && params.run_featurecounts) {
            MAKE_COUNT_MATRIX(FEATURECOUNTS.out.counts, PARSE_METADATA.out.metadata_csv)
        }
        if (params.run_merge_counts && params.run_salmon) {
            MERGE_SALMON(SALMON_QUANT.out.quant_dir.collect())
        }

        // 13. MultiQC
        if (params.run_multiqc) {
            mqc_ch = Channel.empty()
            if (params.run_fastp)             { mqc_ch = mqc_ch.mix(FASTP.out.json) }
            if (params.run_fastq_screen)      { mqc_ch = mqc_ch.mix(FASTQ_SCREEN.out.report) }
            if (use_star)                     { mqc_ch = mqc_ch.mix(STAR_ALIGN.out.log_final) }
            if (use_hisat2)                   { mqc_ch = mqc_ch.mix(HISAT2_ALIGN.out.log_summary) }
            if (params.run_samtools_stats)    { mqc_ch = mqc_ch.mix(SAMTOOLS_FLAGSTAT.out.flagstat)
                                                mqc_ch = mqc_ch.mix(SAMTOOLS_IDXSTATS.out.idxstats) }
            if (params.run_strandedness)      { mqc_ch = mqc_ch.mix(INFER_STRANDEDNESS.out.strand_report.map{it[1]}) }
            if (params.run_coverage_hk)       { mqc_ch = mqc_ch.mix(GENEBODY_COVERAGE_HK.out.coverage) }
            if (params.run_coverage_full)     { mqc_ch = mqc_ch.mix(GENEBODY_COVERAGE_FULL.out.coverage) }
            if (params.run_read_distribution) { mqc_ch = mqc_ch.mix(READ_DISTRIBUTION.out.read_dist) }
            if (params.run_picard)            { mqc_ch = mqc_ch.mix(PICARD_RNA_METRICS.out.metrics) }
            if (params.run_featurecounts)     { mqc_ch = mqc_ch.mix(FEATURECOUNTS.out.summary) }
            if (params.run_salmon)            { mqc_ch = mqc_ch.mix(SALMON_QUANT.out.quant_dir) }
            MULTIQC(mqc_ch.collect())
        }

        // 14-15. DESeq2 / DRomics (after counts are ready)
        if (params.run_deg == true || params.run_deg == "true") {
            if (params.run_featurecounts) {
                DESEQ2_ANALYSIS(FEATURECOUNTS.out.counts, file(params.metadata))
            }
        }
        if (params.run_dromics == true || params.run_dromics == "true") {
            if (params.run_featurecounts) {
                DROMICS_ANALYSIS(FEATURECOUNTS.out.counts, file(params.metadata))
            }
        }
    }
}

workflow.onComplete {
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║  Pipeline Complete!                                      ║
    ║  Status   : ${workflow.success ? '✅ SUCCESS' : '❌ FAILED'}
    ║  Duration : ${workflow.duration}
    ║  Output   : ${params.outdir}
    ╚══════════════════════════════════════════════════════════╝
    """
}
