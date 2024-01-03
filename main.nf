nextflow.enable.dsl=2

params.ldmDir="/path/to/ldm/"
params.stat="/path/to/ma"
params.fscript="format/script"
params.thread=16
params.outdir="outdir"

process prepare_ma {
    disk '500 GB'
    publishDir params.outdir, mode: "copy"

    input:
    path stat
    path ldmDir

    output:
    path "*.ma"

    script:
    """
    ${params.fscript} $stat $ldmDir
    """
}

process generate_LDMList {
    input:
    path ldmDir

    output:
    path "ldm.mlist"
    
    script:
    """
    for i in \$(seq 1 22); do
        echo $ldmDir/*_chr\${i}[_.]*.bin | sed 's/.bin\$//'
    done > ldm.mlist
    """
}

process gctb_sbayesS {
    disk '500 GB'
    publishDir params.outdir, mode: "copy"

    input:
    path ldmDir
    path ldmList
    path ma

    output:
    path "${ma.baseName}.*"

    script:
    """
    gctb --sbayes S --thread ${params.thread} --mldm $ldmList --gwas-summary $ma --out ${ma.baseName} --impute-n
    """
}
    
workflow {
    stat = Channel.fromPath(params.stat)
    fscript = Channel.fromPath(params.fscript)
    ldmDir = Channel.fromPath(params.ldmDir)

    ma = prepare_ma(stat, ldmDir)
    mlist = generate_LDMList(ldmDir)
    gctb_sbayesS(ldmDir, mlist, ma)
}
