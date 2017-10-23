from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import Shared
import os
import gzip
import pysam
from pprint import pprint
from RangeSet import RangeSet
import csv

def generate_reads(fasta_file, read_length):
    fasta_reader = SeqIO.parse(fasta_file, "fasta")
    
    letter_annotations = {"phred_quality": [30] * read_length}
    
    seq_id = 1
    for record in fasta_reader:
        for i in range(len(record) - read_length):
            yield SeqRecord(
                record.seq[i:i+read_length],
                description="",
                id=str(seq_id),
                letter_annotations=letter_annotations
            )
            seq_id += 1

def generate_fastq(fasta_file, read_length, output_file):    
    with gzip.open(output_file, "w") as out_file:
        SeqIO.write(generate_reads(fasta_file, read_length), out_file, "fastq-illumina")

def analyze_hom_regions(bam_file, fasta_file, out_file, threshold=50, feature_db=None):
    bam_reader = pysam.AlignmentFile(bam_file, "rb")
    chrom_names = []
    chrom_lens = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        chrom_names.append(record.id)
        chrom_lens[record.id] = RangeSet([(1, len(record))])
    
    seen = {chrom: [] for chrom in chrom_names}
    mapq_count = {i: 0 for i in range(20)}
    
    low_map_reads = 0
    for read in bam_reader.fetch():
        if read.mapq < 20:
            low_map_reads += 1
            mapq_count[read.mapq] += 1
        else:
            # We use reference_end because Crick reads can align beginning from the other side.
            seen[bam_reader.getrname(read.reference_id)].append((read.reference_start+1, read.reference_end))
            
    all_ranges = {chrom: chrom_lens[chrom] - RangeSet(seen[chrom]) for chrom in chrom_names}
    ranges = {chrom: [r for r in all_ranges[chrom] if r[1] - r[0] >= threshold] for chrom in chrom_names}
    write_ranges({chrom: RangeSet(ranges[chrom]) for chrom in chrom_names}, out_file)
    
    print low_map_reads
    pprint(mapq_count)
    pprint(ranges)
    
    for chrom in chrom_names:
        print chrom
        print "Ignored subranges:", len(all_ranges[chrom])
        print "Total length:", all_ranges[chrom].coverage
        print "Ignored long subranges:", len(ranges[chrom])
        print "Total length:", sum(r[1]-r[0]+1 for r in ranges[chrom])
        print "\n"
        
        if feature_db is not None:
            for r in ranges[chrom]:
                fs = feature_db.get_features_at_range(chrom, r)
                if fs:
                    print chrom, r, ", ".join(f.standard_name for f in fs)
                
            print "\n"
        
def analyze_deletions(bam_file, threshold=50):
    bam_reader = pysam.AlignmentFile(bam_file, "rb")
    
    fasta_file = Shared.get_dependency(os.path.join("albicans", "reference genome", "C_albicans_SC5314_version_A22-s07-m01-r08_chromosomes_HapA.fasta"))
    chrom_names = []
    chrom_lens = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        chrom_names.append(record.id)
        chrom_lens[record.id] = RangeSet([(1, len(record))])
    
    seen = {chrom: [] for chrom in chrom_names}
    
    for read in bam_reader.fetch():
        chrom_name = bam_reader.getrname(read.reference_id)
        if "chrM" in chrom_name:
            continue
        seen[chrom_name].append((read.reference_start+1, read.reference_end-1+1))
            
    unseen = {chrom: chrom_lens[chrom] - RangeSet(seen[chrom]) for chrom in chrom_names}
    write_ranges(unseen, "/Users/bermanlab/dev/transposon-pipeline/dependencies/albicans/deleted_regions.csv")
    ranges = {chrom: [r for r in unseen[chrom] if r[1] - r[0] >= threshold] for chrom in chrom_names}
    
    pprint(ranges)
    
    print "Total unseen:", sum(r.coverage for r in unseen.values())
    print "Total filtered unseen:", sum(sum(r[1]-r[0]+1 for r in rs) for rs in ranges.values())
    
    for chrom in chrom_names:
        print chrom
        print "Total subranges:", len(unseen[chrom])
        print "Total length:", unseen[chrom].coverage
        print "Ignored long subranges:", len(ranges[chrom])
        print "Total length:", sum(r[1]-r[0]+1 for r in ranges[chrom])
        print "\n"
        
        import GenomicFeatures
        alb_db = GenomicFeatures.default_alb_db()
        for r in unseen[chrom]:
            fs = alb_db.get_features_at_range(chrom, r)
            if fs:
                print chrom, r, ", ".join(f.standard_name for f in fs)
                
        print "\n"
        
def write_ranges(range_dict, out_filename):
    with open(out_filename, 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["Chromosome", "Start", "Stop"])
        for chrom in sorted(range_dict.keys()):
            for start, stop in range_dict[chrom]:
                writer.writerow([chrom, start, stop])

def get_read_lengths(bam_file):
    result = {}
    
    bam_reader = pysam.AlignmentFile(bam_file, "rb")
    for read in bam_reader.fetch():
        result[read.qlen] = result.get(read.qlen, 0) + 1
        
    pprint(result)

if __name__ == "__main__":
    import GenomicFeatures
#     generate_fastq(
#         Shared.get_dependency(os.path.join("albicans", "reference genome", "C_albicans_SC5314_version_A22-s07-m01-r08_chromosomes_HapA.fasta")),
#         108,
#         "/tmp/calb_108.fastq.gz"
#     )
#     generate_fastq(
#         Shared.get_dependency(os.path.join("pombe", "Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa")),
#         40,
#         "/tmp/spom_40.fastq.gz"
#     )
#     generate_fastq(
#         Shared.get_dependency(os.path.join("cerevisiae", "S288C_reference_sequence_R64-2-1_20150113.fsa")),
#         75,
#         "/tmp/scer_75.fastq.gz"
#     )
#     analyze_deletions("/Users/bermanlab/ngs-bench/Ella_haploids/BAMs/Galaxy249-[1081.bam_-_final].bam", 50)
#     analyze_hom_regions(
#         "/Users/bermanlab/ngs-bench/tn drawings/out_108.sorted.bam",
#         Shared.get_dependency(os.path.join("albicans", "reference genome", "C_albicans_SC5314_version_A22-s07-m01-r08_chromosomes_HapA.fasta")),
#         "/Users/bermanlab/dev/transposon-pipeline/dependencies/albicans/homologous_regions2.csv",
#         10,
#         GenomicFeatures.default_alb_db()
#     )
#     analyze_hom_regions(
#         "/Users/bermanlab/ngs-bench/tn drawings/Galaxy12-[cer_75].bam",
#         Shared.get_dependency(os.path.join("cerevisiae", "S288C_reference_sequence_R64-2-1_20150113.fsa")),
#         "/Users/bermanlab/dev/transposon-pipeline/dependencies/cerevisiae/homologous_regions.csv",
#         10,
#         GenomicFeatures.default_cer_db()
#     )
#     analyze_hom_regions(
#         "/Users/bermanlab/ngs-bench/Hermes/pom_40.bam",
#         Shared.get_dependency(os.path.join("pombe", "Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa")),
#         "/Users/bermanlab/dev/transposon-pipeline/dependencies/pombe/homologous_regions.csv",
#         10,
#         GenomicFeatures.default_pom_db()
#     )
#     get_read_lengths("/Users/bermanlab/ngs-bench/kornmann raw data/E-MTAB-4885.WT1.bam")
#     get_read_lengths("/Users/bermanlab/ngs-bench/kornmann raw data/E-MTAB-4885.WT2.bam")
#     get_read_lengths("/Volumes/Shared/Vladimir/Transposons/A22-s07-m01-r08/PostEvo/PostEvo03.bam")
#     get_read_lengths("/Volumes/Shared/Vladimir/Transposons/A22-s07-m01-r08/PostEvo/PostEvo07.bam")
#     get_read_lengths("/Volumes/Shared/Vladimir/Transposons/A22-s07-m01-r08/PostEvo/PostEvo11.bam")