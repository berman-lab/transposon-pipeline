import glob, os
import pandas as pd
import numpy as np
ChrFile = '/home/bnet/yaelsilb/delib/data/5314/SC5314_A22_Haploid_A.fa' #move to config file. upload to nas
FeatureFName = '/home/bnet/yaelsilb/delib/data/5314/C_albicans_SC5314_version_A22-s05-m05-r12_chromosomal_feature.tab' #move to config file. upload to nas

SamCols = ['ReadName','FlagSum', 'Chr','Pos','MapQuality','CIGAR','MateChr','MatePos','InsertSize','ReadSeq','ReadQualities',\
           'OptField1','OptField2','OptField3','OptField4','OptField5','OptField6','OptField7','OptField8','OptField9','OptField10']
SamColsTypes = {'ReadName':str,'FlagSum':int, 'Chr':str,'Pos':float,'MapQuality':float,'CIGAR':str,'MateChr':str,'MatePos':int,'InsertSize':int,'ReadSeq':str,'ReadQualities':str,\
           'OptField1':str,'OptField2':str,'OptField3':str,'OptField4':str,'OptField5':str,'OptField6':str,'OptField7':str,'OptField8':str,'OptField9':str,'OptField10':str}
ChrFeatCols = ['FeatureName', 'GeneName','Aliases','FeatureType','Chromosome','StartCoord','StopCoord','Strand','PrimaryCGDID','SecondaryCGDID',\
        'Description','DateCreated','SeqCoordVerDate','Blank1','Blank2','GeneNameReserDate','ReservedIsstandardName','SC_ortholog']

		
#this procedure traverse the sam file, check for high confidence alingment, unite unique positions, and if they should / could be aligned with adjucent positions
def FindHitsPerSample(SamAlign,ChrFeatMap, Sep0N = 2,MapQ=10):
    Chromosomes = SamAlign.Chr.unique()
    SamAlign = SamAlign[(SamAlign.MapQuality>=MapQ)] #consider only reads with a decent map quality 
    ChrHitList = {}
    TotalHits =0;TotalUniqueHits=0
    for Chr in Chromosomes:
        if Chr not in ChrFeatMap:
            print '{} not found in Chromosome feature file'.format(Chr)
            continue
            
        PosCount,bins = np.histogram(SamAlign[SamAlign.Chr==Chr].Pos.values, bins = range(ChrLen[Chr]+1)) #returns number of hits in each position along the chromosome
        HitsPos = np.where(PosCount>0)[0] #returns only the positions in the chromosome which have at least one read
        TotalHits += len(HitsPos)
        #now we want to count how many hits are "unique" after merging close hits with at most N seperating non-reads between them.
        ChrHitList[Chr] = []
        i=0
        
        while i < len(HitsPos): #checking all hits positions
            StartI=i #the index of the ith position with read in the chromosome
            while (i < len(HitsPos)-1) and (HitsPos[i+1] - HitsPos[i] <= Sep0N): #checking all following reads, if they are close to the current - unite them into one merged hit.
                i+=1
            #adding a merged hit with following properties: position (on chromosome) of the first hit found, position of the last hit which was merged, length of the merged section, number of reads found in that merges hit section, chromosome
            ChrHitList[Chr].append((HitsPos[StartI],HitsPos[i], HitsPos[i]-HitsPos[StartI], sum(PosCount[HitsPos[StartI]:HitsPos[i]+1]), Chr))
            i+=1
        
        TotalUniqueHits += len(ChrHitList[Chr])
    return ChrHitList,TotalHits, TotalUniqueHits, SamAlign.shape[0]
       
#this class gets the position of a hit and finds its nearest ORF in the specified strand (W/C) up or down stream of that index 
class NearestORF():
    def __init__(self,Ind, Dist, Strand, UpDown):
        self.Ind = Ind #the index of the ORF in the chromosome feature list
        self.Dist = Dist #the distance between the ORF position and the hit position on the chromowome
        self.Strand = Strand
        self.UpDown = UpDown #is it up or down stream ORF 
    def __gt__(self, ORF2): #greater than function
        return self.Dist > ORF2.Dist
    def __lt__(self, ORF2): #less than function
        return self.Dist < ORF2.Dist
    def GetFeatureName(self):
        if self.Ind >= 0:
            return ChrFeature.loc[self.Ind].FeatureName
        else:
            return 'None' #e.g. if it has no neighbor ORFS downstream - that is, the index is prior to any ORF
    def GetGeneName(self):
        if self.Ind >= 0:
            return ChrFeature.loc[self.Ind].GeneName
        else:
            return 'None'
    def GetFeatureType(self):
        if self.Ind >= 0:
            return ChrFeature.loc[self.Ind].FeatureType
        else:
            return 'None'
    def GetStrand(self):
        if self.Dist == 0:
            return "ORF(" + self.Strand +')'#ORF
        elif self.Ind < 0: # no neighbor found
            return " - "
        else :
            return self.Strand
            

#this procedure find the neares ORF for a position given 
#StartI, StopI: the position on the chormosome of the merged hit; ChrFeatMap: a vector of chromosome length with feature index if exists or 0 in any chromosome position; strand: (W/C)
def  FindNearestORFInStrand(StartI,StopI, ChrFeatMap,Strand):
    if ChrFeatMap[StartI]>=0: #ORF exists at hit position (first merged hit)
        return NearestORF(ChrFeatMap[StartI],0, Strand,'Up'), NearestORF(ChrFeatMap[StartI],0, Strand,'Down') 
    elif  ChrFeatMap[StopI]>=0:  #ORF exists at hit position (last merged hit)
        return NearestORF(ChrFeatMap[StopI],0, Strand,'Up'), NearestORF(ChrFeatMap[StopI],0, Strand,'Down')
    else:
        upi=StartI-1 #upstream is down the index (upi will hold the chromosome position of the nearest ORF
        upF =-1      #upF will hold the feature index of the nearest ORF we found, if exists
        while upi >= 0: #going back in chromosome till we reach its first nucleotide            
            if ChrFeatMap[upi]>=0: #if we find a feature in current position on the chromoeoms
                upF = ChrFeatMap[upi] #we save index of the feature we found 
                break
            upi -= 1
        #again with downstream ORFS
        downi=StopI+1
        DownF = -1
        while downi < len(ChrFeatMap):            
            if ChrFeatMap[downi]>=0:
                DownF = ChrFeatMap[downi]
                break
            downi+=1 
        #returns the closest ORF found up and down stream
        return NearestORF(upF,StartI - upi, Strand,'Up'), NearestORF(DownF, downi-StopI, Strand,'Down')
    
#this procedure takes the mapped hits and for each hit check its nearst ORF and output all the data into the hit file.
def ListHitProp(ChrHitList, FileName,ChrFeatC, ChrFeatW):
    Featuref = open(FileName,'w')
    Featuref.write('Chromosome\tUp feature type\tUp feature name\tUp gene name\tUp feature dist\tDown feature type\tDown feature name\tDown gene name\tDown feature dist\tIntergenicType\tHit position\tHit count\n')
    for Chr in ChrHitList.keys():
        for StartI,StopI,Len,Count,chr_ in ChrHitList[Chr]:
            #find out if we are in ORF, other feature or intergenic
            CORFUp,CORFDown = FindNearestORFInStrand(StartI,StopI, ChrFeatC[Chr], 'C')
            WORFUp,WORFDown = FindNearestORFInStrand(StartI,StopI, ChrFeatW[Chr], 'W')
            #ORF,IntergenicType =  FindClosestGeneInChr(CORF,WORF)
            UpORF = CORFUp if CORFUp< WORFUp else WORFUp
            DownORF = CORFDown if CORFDown< WORFDown else WORFDown
            #now returning a record of this type: StartI, Count, Chr, UpORF, UpDist,DownORF,DownDist, Type (WW,WC,CW,CC)
            Featuref.write('\t'.join([Chr,UpORF.GetFeatureType(), UpORF.GetFeatureName(), str(UpORF.GetGeneName()), str(UpORF.Dist),\
                                      DownORF.GetFeatureType(), DownORF.GetFeatureName(), str(DownORF.GetGeneName()), str(DownORF.Dist),\
                                      (UpORF.GetStrand()+'-'+DownORF.GetStrand()),str(StartI), str(Count)])+'\n')
    Featuref.close()
                        
					
					
usage = """USAGE: MapFastq.py  
   -i  --InDir         [str]   Input directory with .sam files to parse (deafualt is curr dir)
   -o  --OutDir        [str]   Output directory to which the hit file will be writen (deafualt is curr dir)
   -q  --MapQ          [int]   Map Quality - hits to parse from the sam file (default is 10)
   -m  --MergeDist     [int]   Hits to merge with at most x nt distance between two hits (deafult is 2, that is that hit in pos 1 & 3 (3-1=2) will be united into a single hit
   -h, --help                    print this help 
"""

if __name__ == '__main__':
    import sys
    import getopt
    
    SamFileDir = ""; GeneListFileDir = ""; MapQ=10; MergeDist=2;
    
    try:                                
        opts, args = getopt.getopt(sys.argv[1:], "o:i:q:m:h", ["OutDir=","InDir=","MapQ=","MergeDist=","help"])
    except getopt.GetoptError:
        print usage
        sys.exit(2)                     
    for opt, arg in opts:                
        if opt in ("-h", "--help"):      
            print usage
            sys.exit()                  
        elif opt in ("-o", "--OutDir"): 
            GeneListFileDir = arg
        elif opt in ("-i", "--InDir"): 
            SamFileDir = arg
        elif opt in ("-q", "--MapQ"): 
            MapQ = int(arg)
        elif opt in ("-m", "--MergeDist"): 
            MergeDist = int(arg)

	
    #read chromosome Len from fasta file
    ChrLen = {}
    with open(ChrFile,'r') as f:
        for l in f:
            if l[0] == '>' : #>Ca22chr1A_C_albicans_SC5314 (3188363 nucleotides)
                ChrName = l[1:l.find(' ')]
                ChrLen[ChrName]=int(l[l.find('(')+1:l.find('nuc')])
                
    #read all features (ORFs location)
    ChrFeature = pd.read_table(FeatureFName, skiprows =8, names=ChrFeatCols)
    Chromosomes = ChrFeature.Chromosome.unique()


    #Creating the a dictionary of chromosome and its features for the watson strand
    ChrFeatW={} 
    for Chr in Chromosomes:
        #each chromosome have a vector in the chromosome length
        if Chr not in ChrLen.keys(): #e.g. only chromosome B which I analyze
            continue
        ChrFeatW[Chr] = -1*np.ones(ChrLen[Chr])
        Feat = ChrFeature[(ChrFeature.Chromosome==Chr) & (ChrFeature.Strand == 'W') ]
        for row in Feat.iterrows():
            #each position in the chromosome vector is assigned with the feature index relevant to it        
                ChrFeatW[Chr][row[1].StartCoord: row[1].StopCoord+1] = row[0]

    #again, creating the a dictionary of chromosome and its feature but for the crick strand
    ChrFeatC={}
    for Chr in Chromosomes:
        #each chromosome have a vector in the chromosome length
        if Chr not in ChrLen.keys(): #e.g. only chromosome B which I analyze
            continue
        ChrFeatC[Chr] = -1*np.ones(ChrLen[Chr])
        Feat = ChrFeature[(ChrFeature.Chromosome==Chr) & (ChrFeature.Strand == 'C')]
        for row in Feat.iterrows():
            #each position in the chromosome vector is assigned with the feature index relevant to it        
            ChrFeatC[Chr][row[1].StopCoord: row[1].StartCoord+1] = row[0]

    

    if len(GeneListFileDir) > 0 and not os.path.isdir(GeneListFileDir): 
        os.makedirs(GeneListFileDir)


    fNames = glob.glob(os.path.join(SamFileDir, '*.sam'))
    #going over all sam files in a given directory and for each creating a hit file, by uniting close hits and check thier nearest ORFS
    for Name in fNames:
        BaseName = os.path.basename(Name)
        OutFileName = os.path.join(GeneListFileDir, BaseName[:-4] + '_Hits.txt')
        if os.path.isfile(OutFileName): #check first if the file already exists as it is very time consuming. if we want to re create the hit file, just delete or rename it
            continue;
        Sami = pd.read_table(Name, skiprows =19, names=SamCols, dtype = SamColsTypes)
        print "parsing file {}, with {} hits".format(BaseName, Sami.shape[0])
        MyList,TotalHits, TotalUniqueHits, TotalReads = FindHitsPerSample(Sami,ChrFeatW, MergeDist,MapQ) 
        ListHitProp(MyList, OutFileName, ChrFeatC, ChrFeatW ) 
        print "total hits positions (Map quality > {}): {:,}".format(MapQ,TotalHits)
        print "Unique hits position (minimal dist = {} bp) : {:,}".format(MergeDist,TotalUniqueHits)
        print 'total reads found: {:,} (Map quality = {})'.format(TotalReads,MapQ)
	
