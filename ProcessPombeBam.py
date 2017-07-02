import pysam
import argparse
import os
import csv
import GenomicFeatures

# TODO: this eventually should be merged with MapFastq.

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bam", required=True)
    parser.add_argument("-q", type=int, default=20,
                        help="The minimum mapping quality to consider.")
    
    args = parser.parse_args()
    min_mapq = args.q
    bam = pysam.AlignmentFile(args.bam, "rb")
    
    pom_db = GenomicFeatures.default_pom_db()
    hit_map = {'I': {'W': {}, 'C': {}},
               'II': {'W': {}, 'C': {}},
               'III': {'W': {}, 'C': {}}}
    
    for line in bam:
        if line.mapq < min_mapq:
            continue
        
        chrom = bam.getrname(line.reference_id)
        if chrom not in 'III':
            continue
        
        # Since start < end always, in alignments which are reversed (along the
        # Crick strand) the start of the fragment is actually at the 'end' point. 
        if line.is_reverse:
            pos = line.reference_end
            strand = 'C'
        else:
            # BAM files use 0-based indexing, and we work in 1-based indexing,
            # so we have to add one.
            pos = line.reference_start + 1
            strand = 'W'
        
        hit_map[chrom][strand][pos] = hit_map[chrom][strand].get(pos, 0) + 1
    
    with open(os.path.splitext(args.bam)[0] + "_Hits.csv", "wb") as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["Chromosome", "Strand", "Position", "Reads", "Gene"])
        for chrom in sorted(hit_map.keys()):
            for strand in hit_map[chrom].keys():
                for pos in sorted(hit_map[chrom][strand].keys()):
#                     features = pom_db.get_features_at_location(chrom, pos)
                    features = [] # It appears finding the hit gene is not necessary at this point

                    if len(features) > 2:
                        print "More than 1 feature at position", chrom, pos
                    
                    writer.writerow([chrom, strand, pos,
                                     hit_map[chrom][strand][pos],
                                     "nan" if not features else features[0].standard_name])

if __name__ == "__main__":
    main()