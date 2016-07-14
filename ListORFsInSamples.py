import glob, os
import pandas as pd
import numpy as np
from scipy.stats import hypergeom, ranksums
ChrFeaturefName = '/home/bnet/yaelsilb/delib/data/5314/C_albicans_SC5314_version_A22-s05-m05-r12_chromosomal_feature.tab'


#This script goes over all ORFs and list for each ORF the number of hits and reads in each sample

ChrFeatCols = ['FeatureName', 'GeneName','Aliases','FeatureType','Chromosome','StartCoord','StopCoord','Strand','PrimaryCGDID','SecondaryCGDID',\
        'Description','DateCreated','SeqCoordVerDate','Blank1','Blank2','GeneNameReserDate','ReservedIsstandardName','SC_ortholog']
ChrFeature = pd.read_table(ChrFeaturefName, skiprows =8, names=ChrFeatCols)
ChrFeature = ChrFeature[ChrFeature.FeatureType.str.contains('ORF')] #now interesting only in ORFs
ChrFeature = ChrFeature[ChrFeature.Chromosome.str.contains('A_C_albicans_SC5314')] #only features in chromosome A to which we mapped


usage = """USAGE: DrawORFNei.py  
   -i  --InDir              [str]   Input directory containing relevant hit files (defualt current directory)
   -o  --OutFileName        [str]   Output file name, mandatory
   -h, --help                  print this help 
"""

if __name__ == '__main__':
    import sys
    import getopt
    
    InDir = ""; OutFileName = ""; 
    
    try:                                
        opts, args = getopt.getopt(sys.argv[1:], "o:i:h", ["OutFileName=","InDir=","help"])
    except getopt.GetoptError:
        print usage
        sys.exit(2)                     
    for opt, arg in opts:                
        if opt in ("-h", "--help"):      
            print usage
            sys.exit()                  
        elif opt in ("-o", "--OutFileName"): 
            OutFileName = arg
        elif opt in ("-i", "--InDir"): 
            InDir = arg

    if  len(OutFileName)==0:
        print "Output file names is not specified, existing."
        print usage
        sys.exit(2)


    #reading all hit files from the input directory:
    dfs = {}
    fNames = glob.glob(os.path.join(InDir, '*Hits.txt'))
    for Name in fNames:
        BaseName = os.path.basename(Name)
        BaseName = BaseName[:BaseName.find('.')]
        dfs[BaseName] = (pd.read_table(Name))
        #dfs[BaseName]['Hit count'] = dfs[BaseName]['Hit count'].astype(int)
        
        
    #creating the ORF file - for each feature in the feature map iterate over all samples and list number of hits and reads
    f = open(OutFileName,'w')
    f.write('Feature name\tGene name\tStart coord\tStop coord\t')
    for Sample in dfs.keys():
        f.write('Num Of Hits '+Sample+'\tTotal reads ' + Sample +'\tMax reads per hit ' + Sample+'\t') 
    f.write('\n') 
    for l in ChrFeature.iterrows(): # going over each ORF in chromosome feature file
        f.write('\t'.join([l[1].FeatureName, str(l[1].GeneName), str(l[1].StartCoord), str(l[1].StopCoord)]) + '\t')
        for SampleName, SampleFile in dfs.items():
            ORFsRec = SampleFile[(SampleFile['Up feature name']==l[1].FeatureName) & (SampleFile['Down feature name']==l[1].FeatureName)]
            if ORFsRec.shape[0] == 0:
                f.write('0\t0\t0\t')
            else:
                f.write(str(ORFsRec.shape[0]) + '\t' + str(int(ORFsRec['Hit count'].sum())) + '\t' + str(int(ORFsRec['Hit count'].max())) + '\t')
        f.write('\n')
    f.close()

    #now reading the file as a matrix and creating the statistic across all samples
    ORFsSumm = pd.read_csv(OutFileName, sep='\t') #not 80P!!
    ORFsSummV = ORFsSumm > 0
    HitCols = [col for col in ORFsSummV.columns if col.startswith('Num Of Hits')]
    ORFsSumm['SampleWithHits'] = ORFsSummV[HitCols].sum(axis=1) #which samples has at least one hits
    ORFsSumm['AvgHitsPerORFAllSamples'] = ORFsSumm[HitCols].mean(axis=1)
    ORFsSumm['GeneLen'] = (abs(ORFsSumm['Start coord']-ORFsSumm['Stop coord']))
    for Sample in dfs.keys():
        ORFsSumm['HitsPerORFPer100bp' + Sample] = ORFsSumm['Num Of Hits ' + Sample] / ORFsSumm.GeneLen * 100 ## not 80P!!
    ORFsSumm['AvgHitsPerORFPer100bpAllSamples'] = ORFsSumm[[col for col in ORFsSumm.columns if col.startswith('HitsPerORFPer100bp')]].mean(axis=1)
    #ORFsSumm.to_csv(GeneListFileDir + Seq + '_ORFs80PWithSumAllSamples.txt', sep='\t')##80P!!
    ORFsSumm.to_csv(OutFileName, sep='\t')## not 80P!!
