#!/usr/bin/env nextflow

input_fastqs = Channel.fromPath("data/*fq")
file("./tmp").mkdirs()

process read_input_fastq {
  input:
    file in_fq from input_fastqs
  output:
    set val(in_fq), file("catted") into z
  shell:
    """
    cat ${in_fq} > catted
    """
}

process collase_catted {
  publishDir  "tmp", mode: 'copy'
  input:
    set name, catted from z
  output: 
    file "${name}" into finalout
  shell:
    """
    cat ${catted} > ${name}
    """
}
