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
ChrfName = '/home/bnet/yaelsilb/delib/data/5314/SC5314_A22_Haploid_A.fa' #move to configuration file
CentromersfName = '/home/bnet/yaelsilb/delib/data/5314/CentormersA.txt'

#do we use MRS and centromeres? sould be taken from the feature file..
MRS = {'Ca22chrRA_C_albicans_SC5314':[(1,1407),(429937,428531),(1360337,1346420)],'Ca22chr2A_C_albicans_SC5314':\
        [(1,978),(1639716,1638739),(2141592,2155585)],'Ca22chr1A_C_albicans_SC5314':[(1,1527),(259642,261168),(1680765,1695385)],\
        'Ca22chr4A_C_albicans_SC5314':[(1,330),(345406,345735),(1194500,1208753)],'Ca22chr3A_C_albicans_SC5314':[(1,2013),\
        (1641576,1639564)],'Ca22chr5A_C_albicans_SC5314':[(723050,730604)],'Ca22chr6A_C_albicans_SC5314':[(875304,860661)],\
        'Ca22chr7A_C_albicans_SC5314':[(242081,228342),(575236,589526)]}
Centormers = pd.read_csv(CentromersfName, sep = '\t') # left "Left and Right" are for Chr A and right "left and right" are for Chr B

#read chromosome Len from fasta file
ChrLen = {}

with open(ChrfName,'r') as f:
    for l in f:
        if l[0] == '>' : #>Ca22chr1A_C_albicans_SC5314 (3188363 nucleotides)
            ChrName = l[1:l.find(' ')]
            ChrLen[ChrName]=int(l[l.find('(')+1:l.find('nuc')])

RelativeLen = {x:float(ChrLen[x])*100/max(ChrLen.values()) for x in ChrLen.keys()}

ChrsOrder = ['Ca22chrRA_C_albicans_SC5314', 'Ca22chr1A_C_albicans_SC5314','Ca22chr2A_C_albicans_SC5314','Ca22chr3A_C_albicans_SC5314',\
             'Ca22chr4A_C_albicans_SC5314', 'Ca22chr5A_C_albicans_SC5314', 'Ca22chr6A_C_albicans_SC5314','Ca22chr7A_C_albicans_SC5314']
             
        
        
def DrwaChrMapRelative(Hits,SampleName, FileName ='', ylim = 0):
    plt.figure(figsize=(30,4*len(ChrsOrder)))
    plt.gcf().suptitle(SampleName, fontsize=22,fontweight = 'bold', y = 0.92, x = 0.2)
    for i in range(8):#(len(ChrsOrder)):
        ChrHits = Hits[Hits.Chromosome == ChrsOrder[i]]
        ax = plt.subplot2grid((len(ChrsOrder),100), (i,0), colspan=int(RelativeLen[ChrsOrder[i]]))
        #adding the last position as 0/ 1 so the plot won't be shorter    
        ax.bar(ChrHits['Hit position'].values, ChrHits['Hit count'].values)
        ax.set_ylabel(ChrsOrder[i][:ChrsOrder[i].find('_')], fontsize = 20, fontweight = 'bold')
        ax.set_ylim([0,ylim])
        ax.set_xlim([0,ChrLen[ChrsOrder[i]]])
        ax.tick_params(axis='both', which='major', labelsize=14)
    
        #add cenromere
        #xPos = (Centormers[Centormers.Chr == ChrsOrder[i]].Left.values[0] + Centormers[Centormers.Chr == ChrsOrder[i]].Right.values[0] )/2 
        xPos = min(Centormers[Centormers.Chr == ChrsOrder[i]].Left.values[0], Centormers[Centormers.Chr == ChrsOrder[i]].Right.values[0] ) #the location is left and not center!
        xWidth = abs(Centormers[Centormers.Chr == ChrsOrder[i]].Left.values[0]- Centormers[Centormers.Chr == ChrsOrder[i]].Right.values[0] ) 
        Cen = Ellipse(xy=(xPos, 0), width=20000, height=ylim/10,  edgecolor='r', fc='r', lw=2, clip_on=False) #used to be height= 500
        ax.add_patch(Cen) # the axis are not scaled so the x radius is larger than the u axis...
        
        #add MRS
        for startPos,stopPos in MRS[ChrsOrder[i]]: 
            #xPos = ((startPos+stopPos)/2)#*RelativeLen[ChrsOrder[i]]/100
            xPos = min(startPos, stopPos) 
            xWidth = abs(startPos-stopPos) 
            ax.add_patch(Rectangle((xPos,0), width=max(xWidth,5000), height=ylim/10, ec='g', fc='g', clip_on=False))
            
    if len(FileName):
        plt.savefig(FileName)

usage = """USAGE: DrawChrMap.py  
   -i  --InDir         [str]   Input directory with .sam files to parse (deafualt is curr dir)
   -o  --OutDir        [str]   Output directory to which the hit file will be writen (deafualt is curr dir)
   -y  --yLim          [int]   Maximum hits per pos (truncate the hit height on y axis), default is 5,000
   -h, --help                  print this help 
"""

if __name__ == '__main__':
    import sys
    import getopt
    
    MapDir = ""; GeneListFileDir = ""; yLim = 5000;
    
    try:                                
        opts, args = getopt.getopt(sys.argv[1:], "o:i:y:h", ["OutDir=","InDir=","yLim=","help"])
    except getopt.GetoptError:
        print usage
        sys.exit(2)                     
    for opt, arg in opts:                
        if opt in ("-h", "--help"):      
            print usage
            sys.exit()                  
        elif opt in ("-o", "--OutDir"): 
            MapDir = arg
        elif opt in ("-i", "--InDir"): 
            GeneListFileDir = arg
        elif opt in ("-y", "--yLim"): 
            yLim = int(arg)


    if len(MapDir)>0 and not os.path.isdir(MapDir):
        os.makedirs(MapDir)
        
    fNames = glob.glob(os.path.join(GeneListFileDir, '*Hits.txt'))
    for fName in fNames:
        print fName
        BaseName = os.path.basename(fName)
        BaseName = BaseName[:-4]
        OutName = os.path.join(MapDir, BaseName + '_ChrMap_' +str(yLim)+ '.png')
        if os.path.isfile(OutName):
            continue
        Hits = pd.read_table(fName)
        DrwaChrMapRelative(Hits,BaseName, OutName, yLim)
