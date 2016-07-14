import glob, os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Ellipse, Rectangle
from subprocess import PIPE, Popen
import datetime
ChrFeaturefName = '/home/bnet/yaelsilb/delib/data/5314/C_albicans_SC5314_version_A22-s05-m05-r12_chromosomal_feature.tab'             


ChrFeatCols = ['FeatureName', 'GeneName','Aliases','FeatureType','Chromosome','StartCoord','StopCoord','Strand','PrimaryCGDID','SecondaryCGDID',\
        'Description','DateCreated','SeqCoordVerDate','Blank1','Blank2','GeneNameReserDate','ReservedIsstandardName','SC_ortholog']
ChrFeature = pd.read_table(ChrFeaturefName, skiprows =8, names=ChrFeatCols)
ChrFeature = ChrFeature[ChrFeature.FeatureType.str.contains('ORF')] #now interesting only in ORFs
ChrFeature = ChrFeature[ChrFeature.Chromosome.str.contains('A_C_albicans_SC5314')] #only features in chromosome A to which we mapped
        
        
def DrawORFZoom(Hits,SampleName, bp2zoom,ylim, FileName, ORFName):
    plt.figure(figsize=(30,10))
    #searching for ORF in chromosome file
    ORF2Zoom = ChrFeature[ChrFeature.GeneName == ORFName]
    if ORF2Zoom.shape[0] == 0:
        print 'ORF name does not exists in feature file. aborting'
        return
    ChrName = ORF2Zoom.Chromosome.values[0]
    ORFPos = min(ORF2Zoom.StartCoord.values[0],ORF2Zoom.StopCoord.values[0])
    ORFLen = abs(ORF2Zoom.StartCoord.values[0] - ORF2Zoom.StopCoord.values[0])
    CurrChrFeat = ChrFeature[ChrFeature.Chromosome == ChrName][['FeatureName','GeneName','StartCoord','StopCoord']]

    ChrHits = Hits[Hits.Chromosome == ChrName]
    ax = plt.subplot(111)
    ax.bar(ChrHits['Hit position'].values, ChrHits['Hit count'].values)
    ax.set_ylabel('Reads', fontsize = 20, fontweight = 'bold')
    #ax.set_xlabel('nt', fontsize = 20, fontweight = 'bold')
    ax.set_ylim([0,ylim])
    plt.title("{} - {}  (Max reads = {:,}; Distance from ORF point = {:,}".format(SampleName, ORFName, ylim, bp2zoom), fontsize = 26)
    ax.set_xlim([ORFPos - bp2zoom,ORFPos + bp2zoom])
    ax.tick_params(axis='both', which='major', labelsize=14)
    #adding genes :
    ORFs = CurrChrFeat[(abs(CurrChrFeat.StartCoord - ORFPos)<bp2zoom) | (abs(CurrChrFeat.StopCoord - ORFPos)<bp2zoom)] #should also add ORF len for both..
    print ORFs.shape[0], 'ORFs found in that range'
    for l in ORFs.iterrows():
        ORFLoc = min(l[1].StartCoord, l[1].StopCoord)
        ORFLen = abs(l[1].StartCoord - l[1].StopCoord)        
        ORFShape = Rectangle(xy=(ORFLoc, 0), width=ORFLen, height=ylim/100,  edgecolor='g', fc='g', lw=2, clip_on=False) 
        if l[1].GeneName == ORFName:
            ORFShape = Rectangle(xy=(ORFLoc, 0), width=ORFLen, height=ylim/100,  edgecolor='r', fc='r', lw=2, clip_on=False) 
        print l[1].FeatureName, l[1].GeneName,l[1].StartCoord, l[1].StopCoord, ORFLen, ORFLoc
        ax.add_patch(ORFShape)
        ax.text(ORFLoc, -ylim*0.07 ,l[1].GeneName, size=16, fontweight = 'bold')

    plt.savefig(FileName)
        
        
usage = """USAGE: DrawORFNei.py  
   -i  --InFileName         [str]   Input hit file name, mandatory
   -o  --OutFileName        [str]   Output figure file name, mandatory
   -y  --yLim               [int]   Maximum hits per pos (truncate the hit height on y axis), default is 5,000
   -z  --xZoom              [int]   nt to zoom from both sides of the ORF, default is 5,000
   -g  --GeneName           [str]   The ORF to zoom 
   -h, --help                  print this help 
"""

if __name__ == '__main__':
    import sys
    import getopt
    
    InFile = ""; OutFileName = ""; yLim = 5000; Zoom = 5000; ORFName=""
    
    try:                                
        opts, args = getopt.getopt(sys.argv[1:], "o:i:y:z:g:h", ["OutFileName=","InFileName=","yLim=","xZoom=","GeneName=","help"])
    except getopt.GetoptError:
        print usage
        sys.exit(2)                     
    for opt, arg in opts:                
        if opt in ("-h", "--help"):      
            print usage
            sys.exit()                  
        elif opt in ("-o", "--OutFileName"): 
            OutFileName = arg
        elif opt in ("-i", "--InFileName"): 
            InFile = arg
        elif opt in ("-y", "--yLim"): 
            yLim = int(arg)
        elif opt in ("-z", "--xZoom"): 
            Zoom = int(arg)
        elif opt in ("-g", "--GeneName"): 
            ORFName = arg

    if  len(InFile)==0 or not os.path.isfile(InFile) or len(OutFileName)==0 or len(ORFName) ==0:
        print "Input or output files names or ORF bane not specified, existing."
        print usage
        sys.exit(2)
        
    Hits = pd.read_table(InFile)
    DrawORFZoom(Hits,os.path.basename(InFile), Zoom,yLim, OutFileName, ORFName)
    
