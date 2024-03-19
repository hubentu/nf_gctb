nextflow.enable.dsl=2

params.ldmDir="/path/to/ldm/"
params.stat="/path/to/ma"
params.fscript="/format/script"
params.thread=16
params.outdir="outdir"
params.gctb_mem="256.GB"
params.ext="/dev/null"

process prepare_ma {
    publishDir params.outdir, mode: "copy"

    input:
    path fscript
    path stat
    path info
    path ext

    output:
    path "*.ma"

    script:
    """
    chmod +x $fscript
    ./$fscript $stat $info $ext
    rm $stat
    rm $ext
    """
}

process gctb_sbayesS {
    cpus params.thread
    memory params.gctb_mem

    publishDir params.outdir, mode: "copy"

    input:
    path ldmDir
    path ma

    output:
    path "${ma.baseName}.*"

    script:
    """
    ls $ldmDir/*.bin | sed 's/.bin\$//' > ldm.list
    gctb --sbayes S --thread ${params.thread} --mldm ldm.list --gwas-summary $ma --out ${ma.baseName} --impute-n
    """
}
    
workflow {
    stat = Channel.fromPath(params.stat)
    ldmDir = Channel.fromPath(params.ldmDir)
    info = Channel.fromPath("${params.ldmDir}/*hr*.info").collectFile(name: "all.info")
    fscript = Channel.fromPath(params.fscript)
    ext_file = Channel.fromPath(params.ext).ifEmpty(file('/dev/null'))

    ma = prepare_ma(fscript, stat, info, ext_file)
    gctb_sbayesS(ldmDir, ma)
}
