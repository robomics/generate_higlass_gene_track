#!/usr/bin/env nextflow

// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {

    filter_genbank_data(Channel.of(file(params.gene2refseq),
                                   file(params.gene_info),
                                   file(params.gene2pubmed)),
                        params.taxid)
        .branch {
                  gene2refseq: it.getSimpleName().startsWith(file(params.gene2refseq).getSimpleName())
                  geneinfo: it.getSimpleName().startsWith(file(params.gene_info).getSimpleName())
                  gene2pubmed: it.getSimpleName().startsWith(file(params.gene2pubmed).getSimpleName())
                }
        .set { genbank_filtered }

    refgene_sorted = process_refgene(file(params.refgene))

    chrom_sizes = file(params.chrom_sizes)

    gene2pubmed_count = count_citations(genbank_filtered.gene2pubmed)
    geneid_refseqid = process_gene2refseq(genbank_filtered.gene2refseq)

    geneid_refseqid_count = join_gene_with_refseq(geneid_refseqid,
                                                  gene2pubmed_count)
    geneid_refgene_count = join_refseq_with_geneid(geneid_refseqid_count,
                                                   refgene_sorted)
    gene_subinfo_citation_count = join_cit_with_genes(genbank_filtered.geneinfo,
                                                      gene2pubmed_count)

    annotation = generate_annotation_bed(gene_subinfo_citation_count,
                                         geneid_refgene_count)

    run_clodius(chrom_sizes,
                run_exonu(annotation),
                params.assembly)
}

process filter_genbank_data {
    input:
        path file
        val taxid

    output:
        path "*.zst"

    shell:
        '''
        set -o pipefail

        gzip -dc '!{file}'   |
            grep '^!{taxid}' |
            zstd --adapt -T!{task.cpus} -o '!{file.simpleName}.!{taxid}.zst'
        '''
}

process process_refgene {
    input:
        path refgene

    output:
        path "*.zst"

    shell:
        '''
        set -o pipefail

        gzip -dcf '!{refgene}'                        |
            awk -F $'\\t' '{if (!($3 ~ /_/)) print;}' |
            sort -k 2,2                               |
            zstd --adapt -T!{task.cpus} -o 'refgene_sorted.zst'
        '''
}

process count_citations {
    input:
        path gene2pubmed

    output:
        path "*.zst"

    shell:
        '''
        set -o pipefail

        zstd -dcf '!{gene2pubmed}'    |
            awk '{print $2}'          |
            sort                      |
            uniq -c                   |
            awk '{print $2 "\\t" $1}' |
            sort -k 1,1               |
            zstd --adapt -T!{task.cpus} -o 'gen2pubmed_count.zst'
        '''
}

process process_gene2refseq {
    input:
        path gene2refseq

    output:
        path "*.zst"

    shell:
        '''
        set -o pipefail

        zstd -dcf '!{gene2refseq}'           |
            awk -F $'\\t' '{ split($4,a,"."); if (a[1] != "-") print $2 "\\t" a[1];}' |
            sort -k 1,1 -u |
            zstd --adapt -T!{task.cpus} -o 'geneid_refseqid.zst'
        '''
}

process join_gene_with_refseq {
    input:
        path geneid_refseqid
        path gene2pubmed_count

    output:
        path "*.zst"

    shell:
        '''
        set -o pipefail

        join -j 1 \
             <(zstd -dcf '!{geneid_refseqid}')   \
             <(zstd -dcf '!{gene2pubmed_count}') |
             sort -k 2,2                         |
             zstd --adapt -T!{task.cpus} -o geneid_refseqid_count.zst
        '''
}

process join_refseq_with_geneid {
    input:
        path geneid_refseqid_count
        path refgene_sorted

    output:
        path "*.zst"

    shell:
        '''
        set -o pipefail

        join -j 2                                     \
             <(zstd -dcf '!{geneid_refseqid_count}')  \
             <(zstd -dcf '!{refgene_sorted}')         |
             awk 'BEGIN{ OFS="\t" } { print $2,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$3; }' |
             sort -k 1,1                              |
             zstd --adapt -T!{task.cpus} -o geneid_refgene_count.zst
        '''
}

process join_cit_with_genes {
    input:
        path geneinfo
        path gene2pubmed_count

    output:
        path "*.zst"

    shell:
        '''
        set -o pipefail

        join -1 2 -2 1                 \
             <(zstd -dcf '!{geneinfo}' |
               sort -k 2,2 ) \
             <(zstd -dcf '!{gene2pubmed_count}') |
             awk 'BEGIN{ OFS="\\t" } { print $1,$3,$10,$12,$16 }' |
             sort -k 1,1                         |
             zstd --adapt -T!{task.cpus} -o gene_subinfo_citation_count.zst
        '''
}

process generate_annotation_bed {
    input:
        path gene_subinfo_citation_count
        path geneid_refgene_count

    output:
        path "*.bed.zst"

    shell:
        '''
        set -o pipefail

        join <(zstd -dcf '!{gene_subinfo_citation_count}') \
             <(zstd -dcf '!{geneid_refgene_count}')        |
             awk 'BEGIN{ OFS="\\t" } { print $7,$9,$10,$2,$16,$8,$6,$1,$3,$4,$11,$12,$14,$15; }' |
             zstd --adapt -T!{task.cpus} -o gene_annotation.bed.zst
        '''
}

process run_exonu {
    input:
        path gene_annotation

    output:
        path "*.bed"

    shell:
        '''
        set -o pipefail

        exonU.py <(zstd -dcf '!{gene_annotation}') > gene_annotation_exon_unions.bed
        '''
}

process run_clodius {
    publishDir "${params.outdir}", mode: 'copy',
                                   saveAs: { "gene-annotations-${params.assembly}.db" }
    input:
        path chrom_sizes
        path gene_annotation_exon_unions
	    val assembly_name

    output:
        path "annotation.db"

    shell:
        '''
        clodius aggregate bedfile \
            --max-per-tile 20 \
            --importance-column 5 \
	        --assembly '!{assembly_name}' \
	        --chromsizes-filename '!{chrom_sizes}' \
            --output-file annotation.db \
            --delimiter $'\\t' \
            '!{gene_annotation_exon_unions}'
        '''
}
