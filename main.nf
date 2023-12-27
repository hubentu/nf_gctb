nextflow.enable.dsl=2

params.ldmDir="/path/to/ldm/"
params.ma="/path/to/ma"
params.thread=16

process generateLDMList {
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

process run_gctb {
    publishDir ".", mode: "copy"

    input:
    path ldmDir
    path ldmList
    path ma
    var thread

    output:
    path "test.log"

    script:
    """
    gctb --sbayes S --thread $thread --mldm $ldmList --gwas-summary $ma --out $ma --impute-n > test.log
    """
}
    
workflow {
    inputMA = Channel.fromPath(params.ma)
    ldmDir = Channel.fromPath(params.ldmDir)
    mlist = generateLDMList(ldmDir)
    run_gctb(ldmDir, mlist, inputMA)
}
