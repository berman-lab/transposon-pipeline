import glob, os
import pandas as pd
import numpy as np
import pysam
import itertools
import Shared


					
usage = '''CreateHitFile.py  
   -i  --in-dir     [str]   Input directory with .bam files to parse. Defaults to current directory if left unspecified.
   -o  --out-dir    [str]   Output directory to which the hit file will be writen. Defaults to current directory if left unspecified.
   -q  --min-mapq   [int]   Map Quality - hits to parse from the bam file (default is 20)
   -m  --merge-dist [int]   Hits to merge with at most x nt distance between two hits. Default is 2 
                                Example: Hits in positions 1 and 3  (3-1=2) will be united into a single hit
   -h  --help               Show this help message and exit 
'''


# TODO: move to config file.
ChrFile = Shared.get_dependency('albicans', 'reference genome', 'C_albicans_SC5314_version_A22-s07-m01-r08_chromosomes_HapA.fasta')
FeatureFName = Shared.get_dependency('albicans', 'reference genome', 'C_albicans_SC5314_version_A22-s07-m01-r08_chromosomal_feature.tab')

ChrFeatCols = ['FeatureName', 'GeneName','Aliases','FeatureType','Chromosome','StartCoord','StopCoord','Strand','PrimaryCGDID','SecondaryCGDID',\
        'Description','DateCreated','SeqCoordVerDate','Blank1','Blank2','GeneNameReserDate','ReservedIsstandardName','SC_ortholog']

def FindHitsPerSample(SamAlign, ChrFeatMap, Sep0N = 2,MapQ=10):
    """Goes through Sam file, checks for high confidence alignment, unites unique positions if they can be aligned with adjunct positions.

    Parameters
    ----------
        SamAlign    :   x
        ChrFeatMap  :   x
        SepON   :   integer
        MapQ    :   integer
            Setting for map quality. Default of 10 is equal to 1% chance of occurring in another position.

    Returns
    -------
        List of ??
            ChrHitList  :   x 
            TotalHits   :   x 
            TotalUniqueHits :   x 
            total_reads :   x 
    """
    Chromosomes = SamAlign.references
    hit_map = {}

    for line in SamAlign:
        assert line.tid == line.reference_id
        if line.mapq < MapQ:
            continue
        # Since start < end always, in alignments which are reversed (along the Crick strand) the start of the fragment is actually at the 'end' point. 
        if line.is_reverse:
            pos = line.reference_end
            source = 'C'
        else:
            # BAM files use 0-based indexing, and we work in 1-based indexing, so we have to add one.
            pos = line.reference_start + 1
            source = 'W'
        chrom = SamAlign.getrname(line.reference_id)
        if chrom not in hit_map:
            # We store the orientations of the fragments because we only sequence one end of the transposon,
            # so if there is a hit at a particular location, reported by both a Watson- and Crick-aligned fragment, that means that there were (at least) two insertion events, in which the transposon was "flipped".
            hit_map[chrom] = {'W': np.zeros(ChrLen[chrom]+1, dtype=np.int),
                              'C': np.zeros(ChrLen[chrom]+1, dtype=np.int)}
        hit_map[chrom][source][pos] += 1
    
    ChrHitList = {}
    TotalHits = 0
    for (Chr, source) in itertools.product(Chromosomes, ('W', 'C')):
        if Chr not in ChrFeatMap:
            print '{} not found in Chromosome feature file'.format(Chr)
            continue
            
        PosCount = hit_map[Chr][source]
        HitsPos = np.where(PosCount>0)[0] #returns only the positions in the chromosome which have at least one read
        TotalHits += len(HitsPos)
        #now we want to count how many hits are "unique" after merging close hits with at most N seperating non-reads between them.
        if Chr not in ChrHitList:
            ChrHitList[Chr] = []
        i=0
        
        while i < len(HitsPos): #checking all hits positions
            StartI=i #the index of the ith position with read in the chromosome
            #checking all following reads, if they are close to the current - unite them into one merged hit:
            while (i < len(HitsPos)-1) and (HitsPos[i+1] - HitsPos[i] <= Sep0N):
                i+=1
            ChrHitList[Chr].append((HitsPos[StartI], # position (on chromosome) of the first hit found
                                    HitsPos[i], # position of the last hit which was merged
                                    HitsPos[i] - HitsPos[StartI], # length of the merged section
                                    sum(PosCount[HitsPos[StartI]:HitsPos[i]+1]),  # number of reads found in that merged hit section
                                    Chr, # chromosome
                                    source)) # source orientation
            i+=1
    TotalUniqueHits = sum(len(chrom_hits) for chrom_hits in ChrHitList.values())
    total_reads = sum(hits.sum() for source in hit_map.values() for hits in source.values())
    return ChrHitList, TotalHits, TotalUniqueHits, total_reads
 
class NearestORF():
    """Class gets the position of a hit and finds its nearest ORF in the specified strand (W/C) up or down stream of that index

    Parameters
    ----------
        Ind : index of the ORF in the chromosome feature list
        Dist    :   distance between the ORF position and the hit position on the chromosome
        UpDown  : bool giving if up- or downstream ORF
    """

    def __init__(self,Ind, Dist, Strand, UpDown):
        self.Ind = Ind 
        self.Dist = Dist
        self.Strand = Strand
        self.UpDown = UpDown

    # Greater than function
    def __gt__(self, ORF2):
        return self.Dist > ORF2.Dist

    # Less than function
    def __lt__(self, ORF2):
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

def  FindNearestORFInStrand(StartI, StopI, ChrFeatMap, Strand):
    """Finds nearest open reading frame for a given position.

    Parameters
    ----------
        StartI  :   integer
            Start position on chromosome of the merged hit.
        StopI   :   integer
            Stop position on the chromosome of the merged hit.
        ChrFeatMap  :   vector
            Vector of chromsome length with feature index if exists, or 0 in any chromosome position
        Strand  :   boolean
            W (Watson) or C (Crick) strand

    Returns
    -------
        list of lists
            Closest ORF found on upstream (list in first element), or downstream (list in second element)
    """

    # ORF exists at hit position, first merged hit
    if not isinstance(ChrFeatMap[StartI], tuple):
        return NearestORF(ChrFeatMap[StartI],0, Strand,'Up'), NearestORF(ChrFeatMap[StartI],0, Strand,'Down')

    # ORF exists at hit position, last merged hit
    elif not isinstance(ChrFeatMap[StopI], tuple):
        return NearestORF(ChrFeatMap[StopI],0, Strand,'Up'), NearestORF(ChrFeatMap[StopI],0, Strand,'Down')

    else:
        # upstream is down the index
        upi=StartI-1 # upi will hold the chromosome position of the nearest ORF
        upF =-1      # upF will hold the feature index of the nearest ORF we found, if exists
        if not isinstance(ChrFeatMap[upi], tuple):
            upF = ChrFeatMap[upi]
        else:
            new_upi = ChrFeatMap[upi][0]
            if new_upi != -1:
                upi = new_upi
                upF = ChrFeatMap[upi]
            else:
                upi = 1

        # same procedure with downstream ORFS
        downi=StopI+1
        DownF = -1
        if not isinstance(ChrFeatMap[downi], tuple):
            DownF = ChrFeatMap[downi]
        else:
            new_downi = ChrFeatMap[downi][1]
            if new_downi != -1:
                downi = new_downi 
                DownF = ChrFeatMap[downi]
            else:
                downi = len(ChrFeatMap)

        # returns the closest ORF found up and down stream
        return NearestORF(upF,StartI - upi, Strand,'Up'), NearestORF(DownF, downi-StopI, Strand,'Down')

def ListHitProp(ChrHitList, FileName, ChrFeatC, ChrFeatW):
    """Takes mapped hits and for each hit checks nearest ORF and outputs data into hit file
    """
    Featuref = open(FileName,'w')
    Featuref.write('Chromosome\tSource\tUp feature type\tUp feature name\tUp gene name\tUp feature dist\tDown feature type\tDown feature name\tDown gene name\tDown feature dist\tIntergenicType\tHit position\tHit count\n')
    for Chr in ChrHitList.keys():
        for StartI,StopI,Len,Count,chr_,source in ChrHitList[Chr]:
            #find out if we are in ORF, other feature or intergenic
            CORFUp,CORFDown = FindNearestORFInStrand(StartI,StopI, ChrFeatC[Chr], 'C')
            WORFUp,WORFDown = FindNearestORFInStrand(StartI,StopI, ChrFeatW[Chr], 'W')
            UpORF = CORFUp if CORFUp < WORFUp else WORFUp
            DownORF = CORFDown if CORFDown< WORFDown else WORFDown
            #now returning a record of this type: 
            #   StartI, Count, Chr, UpORF, UpDist,DownORF,DownDist, Type (WW,WC,CW,CC)
            Featuref.write('\t'.join([Chr,source,UpORF.GetFeatureType(), UpORF.GetFeatureName(), str(UpORF.GetGeneName()), str(UpORF.Dist),\
                                      DownORF.GetFeatureType(), DownORF.GetFeatureName(), str(DownORF.GetGeneName()), str(DownORF.Dist),\
                                      (UpORF.GetStrand()+'-'+DownORF.GetStrand()),str(StartI), str(Count)])+'\n')
    Featuref.close()


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(usage=usage)
    
    parser.add_argument("-o", "--out-dir", default='.')
    parser.add_argument("-i", "--in-dir", default='.')
    parser.add_argument("-q", "--min-mapq", default=20, type=int)
    parser.add_argument("-m", "--merge-dist", default=2, type=int)
    
    args = parser.parse_args()
    
    SamFileDir = args.in_dir
    GeneListFileDir = args.out_dir
    MapQ = args.min_mapq
    MergeDist = args.merge_dist
    
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
            ChrFeatW[Chr][int(row[1].StartCoord): int(row[1].StopCoord)+1] = row[0]

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
            ChrFeatC[Chr][int(row[1].StopCoord): int(row[1].StartCoord)+1] = row[0]

    # Preprocess the ChrFeat maps to include (left, right) tuples of the nearest features for every index:
    for feat_map in (ChrFeatC, ChrFeatW):
        for chrom_name in feat_map:
            chrom = map(int, feat_map[chrom_name])
            prev_feat_ix = next_feat_ix = -1
            i = 1
            while i < len(chrom):
                if chrom[i] != -1:
                    prev_feat_ix = i
                    i += 1
                else:
                    # Find the next feature!
                    next_feat_ix = -1
                    j = i+1
                    while j < len(chrom):
                        if chrom[j] != -1:
                            next_feat_ix = j
                            break
                        j += 1
                    for k in range(i, j):
                        chrom[k] = (prev_feat_ix, next_feat_ix)
                    i = j
            feat_map[chrom_name] = chrom


    if len(GeneListFileDir) > 0 and not os.path.isdir(GeneListFileDir): 
        os.makedirs(GeneListFileDir)
    
    # SAM is 1-based and BAM is 0-based, so we work with BAM only out of convenience.
    fNames = glob.glob(os.path.join(SamFileDir, '*.bam'))
    #going over all sam files in a given directory and for each creating a hit file, by uniting close hits and check thier nearest ORFS
    for Name in fNames:
        BaseName = os.path.basename(Name)
        OutFileName = os.path.join(GeneListFileDir, BaseName[:-4] + '_Hits.txt')
        if os.path.isfile(OutFileName): 
        #check first if the file already exists as it is very time consuming. if we want to re create the hit file, just delete or rename it
            print 'ERROR: Hit file for %s already exists.' % (BaseName)
            continue;
        Sami = pysam.AlignmentFile(Name, "rb")
        print "Parsing file %s..." % (BaseName)
        
        MyList,TotalHits, TotalUniqueHits, TotalReads = FindHitsPerSample(Sami,ChrFeatW, MergeDist,MapQ) 
        ListHitProp(MyList, OutFileName, ChrFeatC, ChrFeatW) 
        UniqueHitPercent = round(float(TotalUniqueHits)/float(TotalHits)*100, 2)
        
        Log = '\r\n=== Finding hits ===\r\n%s reads found of map quality >= %s; in these:\r\n  %s hits were found; of these:\r\n    %s (%s%%) hit positions were found to be unique (Minimal distance = %s)\r\n' % (TotalReads, MapQ, TotalHits, TotalUniqueHits, UniqueHitPercent, MergeDist)
        print Log
        
        LogFile = open('log.txt', 'a')
        LogFile.write(Log)
        LogFile.close()