#!/usr/bin/env python

import sys
import re
import os
import shutil
import subprocess
from argparse import ArgumentParser
sys.path.append('/home/ubuntu/arjun/toil-scripts/src/')
from toil_scripts.lib.programs import docker_call

def fai_chunk(faipath, blocksize):
    """
    yield chr, start, end of sequence
    """
    seq_map = {}
    with open( faipath ) as handle:
        for line in handle:
            tmp = line.split("\t")
            seq_map[tmp[0]] = long(tmp[1])

    for seq in seq_map:
        l = seq_map[seq]
        for i in xrange(1, l, blocksize):
            #yield seq, i, min(i+blocksize-1, l)
            yield '{}:{}-{}'.format(seq, i, min(i+blocksize-1, l))


def run_muse(workdir, genome, normal, tumor, block_size, dbsnp, outvcf):
    """
    Runs muse in parallel on chunks of the genome
    If genome and bam files are not indexed, runs samtools to create them.
    and dbsnp VCF for running muse sump
    """
    
    # check that the index files exist or make them
    genome_idx = genome + '.fai'
    if not os.path.exists(genome_idx):
        subprocess.check_call( ['samtools', 'faidx', genome_idx] )

    normal_idx = normal + '.bai'
    if not os.path.exists(normal_idx):
       subprocess.check_call( ['samtools', 'index', normal] )

    tumor_idx = tumor + '.bai'
    if not os.path.exists(tumor_idx):
       subprocess.check_call( ['samtools', 'index', tumor] )

    # create the commands
    outfiles = []
    output_base = 'output.file'
    for block_num, block in enumerate(fai_chunk( genome_idx, block_size ) ):
        parameters = ['call', 
                   '-f', genome, 
                   '-r', block, 
                   '-O', '{}.{}'.format(output_base, block_num),
                   tumor, normal]

        outfiles.append('{}.{}.MuSE.txt'.format(output_base, block_num))

# NOTE TO ARJUN: docker_calls are not paralellized in this test script
        docker_call(tool='aarjunrao/muse',
                work_dir=workdir, parameters=parameters )

    # merge the resulting files in the same order as run
#    outfiles = outfiles[:1]
    first = True
    merge = 'merge.output'
    with open(merge, "w") as ohandle:
        for f in outfiles:
            with open(f) as handle:
                for line in handle:
                    if first or not line.startswith("#"):
                        ohandle.write(line)
            first = False
            #os.unlink(f)

    # run MuSE sump on the merged file, with SNP file if available
    # The -E parameter means 'whole exome' (as opposed to -G, for genome)
    parameters = ['sump',
              '-I', merge,
#              '-D', dbsnp,
              '-O', outvcf,
              '-E']
    docker_call(tool='aarjunrao/muse',
                work_dir=workdir, parameters=parameters)
     



if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-g', '--genome', help='faidx indexed reference sequence file', required=True)
    parser.add_argument('-b', '--blocksize', type=long, help='Parallel Block Size', default=50000000)
    parser.add_argument('-o', '--outvcf', help='output file name (VCF)', default='out.vcf')
    parser.add_argument('-s', '--snpvcf', help='dbSNP vcf file')

    parser.add_argument('-t', '--tumor_bam', required=True)
    parser.add_argument('-n', '--normal_bam', required=True)
    args = parser.parse_args()

    workdir = '/mnt/data/muse'
    run_muse(workdir, args.genome, args.normal_bam, args.tumor_bam, args.blocksize, args.snpvcf, args.outvcf)
