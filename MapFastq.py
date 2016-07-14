import glob, os
import zlib
import pandas as pd
import gzip
from Bio.Seq import Seq
import numpy as np
from Bio.Alphabet import generic_dna
import subprocess
from subprocess import PIPE, Popen
CutAdaptCMD = '/usr/local/lib/python2.7/dist-packages/cutadapt/bin/cutadapt'
TnPrimerAndTail = 'GTATTTTACCGACCGTTACCGACCGTTTTCATCCCTA'
BowtiePath = '/specific/fprivate/bnet/yaelsilb/delib/bowtie2-2.2.6/'
CInd = '/home/bnet/yaelsilb/delib/data/5314/5314_s05-m05-r12'


def GetReads(Val,Text):
    Sub = Text[Text.find(Val)+len(Val)+1:].lstrip(' ')
    Sub = Sub[:Sub.find('\n')]
    Vals = Sub.split(' ')
    return int(Vals[0].replace(',',''))
    
def ProcessOutput(output):
    Summary = output[output.find('=== Summary ===')+17:]
    TotalReads = GetReads('Total read pairs processed',Summary)
    Adapter = GetReads('Reads with adapters',Summary)
    Total = GetReads('Reads written (passing filters)',Summary)
    return TotalReads, Adapter, Total

def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()

def RemoveTn(InputFileName, TnSeq, overlap, OutputFName):
    #-g is from the beggining and -a is from the end or cive versa (if we use -a we get only a string length 6, which is probably the 6 Ns between the primer and Tn)
    #-o out put, afterwards the input. --overlap determines the minimal number of bases found from the Tn
    cmdl = CutAdaptCMD + ' -g ' + TnSeq + ' -o '+ OutputFName + ' ' + InputFileName + ' --discard-untrimmed --overlap ' + str(overlap)
    Output,err = cmdline(cmdl)
    return ProcessOutput(Output)
    
    

usage = """USAGE: MapFastq.py  
   -o  --OutDir          [str]   Output fir, if not specified creating all file in curr dir.
   -i  --InputFileName   [str]   Input fastq file name 
   -h, --help                    print this help 
"""

if __name__ == '__main__':
    import sys
    import getopt
    
    OutputDir = ""; FastqFName = ""; 
    
    try:                                
        opts, args = getopt.getopt(sys.argv[1:], "o:i:h", ["OutDir=","InputFileName=","help"])
    except getopt.GetoptError:
        print usage
        sys.exit(2)                     
    for opt, arg in opts:                
        if opt in ("-h", "--help"):      
            print usage
            sys.exit()                  
        elif opt in ("-o", "--OutDir"): 
            OutputDir = arg
        elif opt in ("-i", "--InputFileName"): 
            FastqFName = arg

    if  len(FastqFName)==0 or len(CInd) ==0:
        print "Input file or index file not specified, existing."
        print usage
        sys.exit(2)

    CleanfName = os.path.join(OutputDir, os.path.basename(FastqFName) + '.clean.fq')
    CurrRes = RemoveTn(FastqFName,TnPrimerAndTail,37, CleanfName)
    #aligning fastq file to the reference genome
    #-x reference_genome; -U fq file of unpaired reads; -S output SAM alignment file; X is the max fragment size to consider (relevant only in paired end)
    SamfName = os.path.join(OutputDir, os.path.splitext(os.path.basename(FastqFName))[0] + '.sam')
    BowtieCmd = BowtiePath + 'bowtie2 -x '+ CInd + ' -U ' + CleanfName + ' -S ' + SamfName
    Output,err = cmdline(BowtieCmd)
    # converting to bam file 
    cmdline('samtools view -bS '+ SamfName + ' > '+ os.path.splitext(SamfName)[0]+".bam")
    #sorting and indexing 
    cmdline('samtools sort '+os.path.splitext(SamfName)[0]+'.bam ' + os.path.splitext(SamfName)[0]+'.sorted')
    cmdline('samtools index '+os.path.splitext(SamfName)[0]+'.sorted.bam')
