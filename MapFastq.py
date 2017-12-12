import os
from subprocess import PIPE, Popen
import Shared

usage = '''MapFastq.py  
   -o  --out-dir            [str]   Output directory. Defaults to current directory if left unspecified.
   -i  --input-fastq-file   [str]   REQUIRED. Input fastq file (Need to include path and '.fastq.gz' at end of filename)
   -t  --transposon-seq     [str]   Transposon sequence, default is Berman Tn with primer
   -a  --clean-adapters     [bool]  Clean Illumina universl adapters. Default is False
   -d  --delete-originals   [bool]  Delete input FASTQ files. Default is False
   -k  --keep-clean-fqs     [bool]  Keep the cleaned FASTQ files. Default is False
   -r  --run-reverse        [bool]  Run program for reverse read (R2) sequence data. Default is False
   -h  --help                       Show this help message and exit 
'''

TnPrimerAndTail = 'GTATTTTACCGACCGTTACCGACCGTTTTCATCCCTA'
AdapterSeq = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

RevTn = 'TAGGGATGAAAACGGTCGGTAACGGTCGGTAAAATAC'
RevAdapter = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'

# NB: bowtie2 requires spaces to be escapes with a backslash for the -x parameter.
CInd = Shared.get_dependency("albicans", "reference genome", "C_albicans_SC5314_version_A22-s07-m01-r08_chromosomes_HapA").replace(' ', '\ ')
CORES = 4 # How many cores on the machine = how many threads should the external tools utilize

CutadaptShellPath = os.popen('command -v cutadapt').read()
CutAdaptPath = CutadaptShellPath.strip('\n')

BowtieShellPath = os.popen('command -v bowtie2').read()
BowtiePath = BowtieShellPath.strip('\n')

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
    #-g is from the beggining and -a is from the end or vice versa (if we use -a we get only a string length 6, which is probably the 6 Ns between the primer and Tn)
    #-o out put, afterwards the input. --overlap determines the minimal number of bases found from the Tn
    CutAdaptCmd = CutAdaptPath + ' -m 2 -g ' + TnSeq + ' -o "'+ OutputFName + '" "' + InputFileName + '" --discard-untrimmed --overlap ' + str(overlap)
    Output,err = cmdline(CutAdaptCmd)
    return ProcessOutput(Output)
    

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--out-dir", default='.')
    parser.add_argument("-i", "--input-file-name", required=True)
    parser.add_argument("-a", "--clean-adapters", default=True, action='store_true')
    parser.add_argument("-d", "--delete-originals", default=False, action='store_true')
    parser.add_argument("-k", "--keep-clean-fastqs", default=False, action='store_true')
    parser.add_argument("-r", "--run-reverse", default=False, action='store_true')
    
    args = parser.parse_args()
    
    OutputDir = args.out_dir
    FastqFName = args.input_file_name
    clean_adapters = args.clean_adapters
    delete_originals = args.delete_originals
    keep_clean_fqs = args.keep_clean_fastqs
    run_reverse = args.run_reverse
    
    # setting sequennces to run on reverse read sequence file
    if run_reverse:
        TnPrimerAndTail = RevTn
        AdapterSeq = RevAdapter
    
    # NB: we're first removing the transposon head from the beginning and then 
    # removing the adapter. This is because the sequencing tech, for whatever reason,
    # removes the adapater from the 5' Tn end, but keeps the tail adapters (sometimes),
    # presumably if the reads are too short.
    # After we have all of the reads that have a transposon, we then trim the adapter from the tail.
    CleanfName = os.path.join(OutputDir, os.path.basename(FastqFName) + '.clean.fq')
    CurrRes = RemoveTn(FastqFName,TnPrimerAndTail,37, CleanfName)
    
    if clean_adapters:
        temp_fq = CleanfName + ".no_adapters.fq"
        cmdline(r'''%s -m 2 -a %s -o "%s" "%s"''' %
                (CutAdaptPath, AdapterSeq, temp_fq, CleanfName))
        os.unlink(CleanfName)
        os.rename(temp_fq, CleanfName)
    
    #aligning fastq file to the reference genome
    #-x reference_genome; -U fq file of unpaired reads; -S output SAM alignment file; X is the max fragment size to consider (relevant only in paired end)
    SamfName = os.path.join(OutputDir, os.path.splitext(os.path.basename(FastqFName))[0] + '.sam')
    BowtieCmd = BowtiePath + ' -p ' + str(CORES) + ' -x "'+ CInd + \
                '" -U "' + CleanfName + '" -S "' + SamfName + '"'
    Output,err = cmdline(BowtieCmd)
    # converting to bam file 
    cmdline('samtools view -@ ' + str(CORES) + ' -bS "'+ SamfName + '" > "' +
            os.path.splitext(SamfName)[0] + '.bam"')
    #sorting and indexing 
    cmdline('samtools sort -@ ' + str(CORES) + ' "' +
            os.path.splitext(SamfName)[0] + '.bam" > "' +
            os.path.splitext(SamfName)[0] + '.sorted.bam"')
    cmdline('samtools index "' + os.path.splitext(SamfName)[0] + '.sorted.bam"')

    # Remove intermediats:
    if delete_originals:
        os.remove(FastqFName)
    if not keep_clean_fqs:
        os.remove(CleanfName)
    os.remove(SamfName)
    os.remove(os.path.splitext(SamfName)[0] + '.bam')
    
