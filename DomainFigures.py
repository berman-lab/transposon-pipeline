from math import floor, ceil
import cairocffi as cairo
import argparse
import SummaryTable
import glob
import os
import shlex
from SortedCollection import SortedCollection
from RangeSet import RangeSet
import Organisms
from Organisms import get_orths_by_name

def _get_organism(org_name):
    return {"Calb": Organisms.alb, "Scer": Organisms.cer, "Spom": Organisms.pom}[org_name]

GENES_ALL, GENES_HIGHLIGHTED, GENES_NONE = ("all", "highlighted", "none")

MAINUSEAGE = '''DomainFigures.py
    Draws genomic regions with hit locations from data file(s). Can also highlight essential domains and directions of genes.

    REQUIRES:   One positional argument (region, gene, or config_file) and any relevant subarguments. 
        Details via -h on desired argument. ('python DomainFigures.py [POSITIONAL] -h')
 
    OPTIONAL:   Non-positional arguments (placed before the positional argument)
        ('python DomainFigures.py --[NON-POSITIONAL] [POSITIONAL] --[SUBARGUMENT]')


    NON-POSITIONAL ARGUMENTS:
        --hits-dir          [str]   Input folder for hits file(s). Defaults to current directory.
        --output-dir        [str]   Output folder for image(s). Defaults to current directory.
        --domains           [str]   Which genes to draw ess domains on. 
                                    Choose between: all, highlighted, none (Default is 'highlighted') 
        --direction         [str]   Which genes to draw read direction for. 
                                    Choose between: all, highlighted, none (Default is 'highlighted') 
        --organism          [str]   Choose between: Calb, Scer, Spom (Default is 'Calb')
        --absolute-pixel-size [int] Draws figure length relative to length of region being drawn.
                                Use for multiple images with comparable sizes. (Default is off with all figures 250 px long)
'''

REGIONUSAGE = '''Domainfigures.py region
    Define drawn area via chromosomal coordinates.

    REQUIRES (all three):
        --chromosome    [str]   Which chromosome to draw. For Calb choose: number '1' through '7' or 'R'
                                (Uses SC5314 assembly 22, haplotype A)
        --start         [int]   bp position from which to start figure.
        --stop          [int]   bp position from which to stop figure.

    OPTIONAL:
        --genes     Choose gene(s) to highlight (for domains, direction arguments). Use standard name(s), 
                    or '*' for all genes. ('*' equivalent to 'all' for domains, direction arguments)
'''

GENEUSAGE = '''Domainfigures.py gene
    Define drawn area via gene name(s).

    REQUIRES:
        --genes     Choose gene(s) to draw and highlight. Use standard name(s), or '*' for all genes.

    PLUS ONE OF:    Defines flanking region to draw around drawn gene(s).
        --percent-of-length     [float]     Percent of each gene's length. Default is 0.2 (20 percent)
        --bps                   [int]       Basepairs before and after gene(s). Default is 20000.
'''

CONFIGUSAGE = '''Domainfigures.py config_file
    Talk to Yael and Vladimir.
'''

def main():
    ''' Sets arguments and subarguments for running the program, and reads in files for organism specified.
        If config_file being used, reads that in too.
    '''
    parser = argparse.ArgumentParser(usage=MAINUSEAGE)
    
    parser.add_argument("--hits-folder", default=".")
    parser.add_argument("--output-folder", default=".")
    parser.add_argument("--domains", default="highlighted", choices=[GENES_ALL, GENES_HIGHLIGHTED, GENES_NONE])
    parser.add_argument("--direction", default="highlighted", choices=[GENES_ALL, GENES_HIGHLIGHTED, GENES_NONE])
    parser.add_argument("--organism", default="Calb", choices=["Calb", "Scer", "Spom"])
    parser.add_argument("--absolute-pixel-size", type=int, default=0)
    
    gene_list_parser = lambda gs: [g for g in gs.split(',')]
    
    subparsers = parser.add_subparsers(dest="source_type")
    
    region_parser = subparsers.add_parser("region", usage=REGIONUSAGE)
    region_parser.add_argument("--chromosome", required=True)
    region_parser.add_argument("--start", type=int, required=True)
    region_parser.add_argument("--stop", type=int, required=True)
    region_parser.add_argument("--genes", type=gene_list_parser, default="")
    
    gene_parser = subparsers.add_parser("gene", usage=GENEUSAGE)
    gene_name = gene_parser.add_mutually_exclusive_group(required=True)
    # TODO: can we make --genes not be named, and instead come at the end of the parser?
    gene_name.add_argument("--genes", type=gene_list_parser)
    gene_region = gene_parser.add_mutually_exclusive_group()
    gene_region.add_argument("--percent-of-length", type=float, default=0.2)
    gene_region.add_argument("--bps", type=int)
    
    config_file_parser = subparsers.add_parser("config_file", usage=CONFIGUSAGE)
    config_file_parser.add_argument("config_file")
    
    args = parser.parse_args()
    
    if args.organism == "Calb":
        hits = SummaryTable.read_hit_files(glob.glob(os.path.join(args.hits_folder, "*_Hits.txt")))
    elif args.organism == "Scer":
        import cPickle
        all_track_files = glob.glob(os.path.join(args.hits_folder, "*.wig"))
            
        # We cache the hits because the hit reading process involves finding the
        # feature which got hit, for every hit, and that makes reading an O(n log n)
        # operation, which is a little slow. O(n) is better here.
        hit_cache = os.path.join(args.hits_folder, "cached_sc_track_hits.dat")
        if not os.path.exists(hit_cache):
            all_tracks = [SummaryTable.get_hits_from_wig(fname) for fname in all_track_files]
            with open(hit_cache, 'wb') as pickle_file:
                cPickle.dump(all_tracks, pickle_file)
        else:
            with open(hit_cache, 'rb') as pickle_file:
                all_tracks = cPickle.load(pickle_file)
                
        hits = all_tracks
    elif args.organism == "Spom":
        hits = [SummaryTable.read_pombe_hit_file("/Users/bermanlab/ngs-bench/Hermes/SRR327340.trimmed.trail_q_20.sorted_Hits.csv")]
    
    # Process hits for quicker access by chromosome name and position:
    new_hits = []
    db = _get_organism(args.organism).feature_db
    
    # TODO: we manipulate the chromosome names to reflect the standard names,
    # but this should be done in the hit-reading functions.
    chrom_names = db._chrom_names # TODO: don't do this, _chrom_names are protected!
    print chrom_names
    import sys; sys.exit()
    for hit_track in hits:
        new_hit_track = {chrom.name: SortedCollection(key=lambda h: h["hit_pos"]) for chrom in db}
        for hit in hit_track:
            chrom = chrom_names[hit["chrom"]]
            hit["chrom"] = chrom
            new_hit_track[chrom].insert(hit)
        new_hits.append(new_hit_track)
    hits = new_hits
    
    if not hits:
        raise Exception("No hit files were found in the hits folder: %s" % args.hits_folder)
    
    if args.source_type == "config_file":
        with open(args.config_file, "r") as in_file:
            # Read the file once, in case it changes in the middle of the run:
            for line in in_file.readlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                handle_args(parser.parse_args(shlex.split(line)), hits)
    else:
        handle_args(args, hits)
    
def handle_args(args, hits):
    '''Takes input arguments and data from hit file(s). Passes into correct drawing module.
    '''
    organism = _get_organism(args.organism) 
    db = organism.feature_db
    
    if args.source_type == "gene":
        if '*' in args.genes:
            gene_names = [f.standard_name for f in db.get_all_features() if f.is_orf]
        else:
            gene_names = args.genes
        
        features = map(db.get_feature_by_name, gene_names)
        
        missed_genes = [name for f, name in zip(features, gene_names) if f is None]
        if missed_genes:
            print "The following genes weren't found and will be skipped: %s" % ', '.join(missed_genes)
            features = filter(None, features)
            
        all_orthologs = map(get_orths_by_name, [f.standard_name for f in features])
        
        for feature, (alb_orth, cer_orth, pom_orth) in zip(features, all_orthologs):
            alb_name = None
            if alb_orth:
                alb_name_components = [alb_orth.standard_name]
                if alb_orth.common_name and alb_orth.common_name != alb_orth.standard_name:
                    alb_name_components.append(alb_orth.common_name)
                
                alb_name = "Calb-" + "_".join(alb_name_components)
            
            cer_name = None
            if cer_orth:
                cer_name_components = [cer_orth.feature_name]
                if cer_orth.common_name and cer_orth.common_name != cer_orth.feature_name:
                    cer_name_components.append(cer_orth.common_name)
                    
                cer_name = "Scer-" + "_".join(cer_name_components)
                    
            pom_name = None
            if pom_orth:
                pom_name_components = [pom_orth.standard_name]
                if pom_orth.common_name and pom_orth.common_name != pom_orth.standard_name:
                    pom_name_components.append(pom_orth.common_name)
                    
                pom_name = "Spom-" + "_".join(pom_name_components)
                    
            if isinstance(organism, Organisms.Calbicans):
                name = "-".join(filter(None, (alb_name, cer_name, pom_name))) 
            elif isinstance(organism, Organisms.Scerevisiae):
                name = "-".join(filter(None, (cer_name, alb_name, pom_name))) 
            elif isinstance(organism, Organisms.Spombe):
                name = "-".join(filter(None, (pom_name, alb_name, cer_name)))
            
            name += "-" + str(len(feature))
            
            label = "%s / %s" % (feature.name, len(feature))
            
            chromosome = feature.chromosome
            
            if args.bps is not None:
                gene_pad = args.bps
            else: # args.percent_of_length exists
                gene_pad = int(args.percent_of_length * len(feature))
            
            draw_gene(organism, feature, gene_pad, label, hits, args.domains, args.direction, name, args.output_folder, args.absolute_pixel_size)
    else: # "region"
        chromosome = args.chromosome
        start = args.start
        stop = args.stop
        name = "Chr%s_%d-%d.png" % (chromosome, start, stop)
        
        draw_genomic_region(
            organism,
            organism.feature_db.get_std_chrom_name(chromosome),
            start, stop, hits,
            args.domains,
            args.direction,
            os.path.join(args.output_folder, name),
            highlighted_genes=set(args.genes),
            absolute_pixel_size=args.absolute_pixel_size
        )
    
def draw_gene(organism, gene, gene_pad, label, hits, draw_domains, draw_directions, out_file_prefix, out_dir, absolute_pixel_size=None):  
    '''Takes drawing settings from handle_args module for gene source type, and draws figure.

    Parameters
    ----------
        organism    :       organism genome is of
        gene    :           feature being drawn
        gene_pad    :       amount of flanking region (defined by either percent of length, or bps)
        label   :           feature name and length (to label figure)
        hits    :           hit data for feature
        draw_domains    :   domain argument
        draw_directions :   direction argument
        out_file_prefix :   '-[feature length]'
        out_dir :           output directory
        absolute-pixel-size: absolute-pixel-size argument

    Writes
    ------
        Gene figure(s) to png image file(s) via draw_genomic_region
    '''
   
    region_start = max(floor(gene.start - gene_pad), 0)
    region_end = min(ceil(gene.stop + gene_pad), len(organism.feature_db[gene.chromosome]))

    # Output to PNG
    name = "%s.png" % out_file_prefix

    draw_genomic_region(
        organism,
        gene.chromosome, region_start, region_end,
        hits,
        draw_domains,
        draw_directions,
        os.path.join(out_dir, name),
        highlighted_genes=set([gene.standard_name]),
        label=label,
        absolute_pixel_size=absolute_pixel_size
    )

def draw_genomic_region(organism, chromosome, region_start, region_end, hits, draw_domains, draw_directions, out_file,
        highlighted_genes=frozenset(), label=None, absolute_pixel_size=0):
    '''Takes drawing settings from handle_args module for region or draw_gene, and draws figure.

    Parameters
    ----------
        organism    :       organism genome is of
        chromosome  :       chromosome for feature being drawn
        region_start:       bp position to start drawing figure on chromosome
        region_end  :       bp position to end drawing figure on chromosome
        hits    :           hit data for feature
        draw_domains    :   domain argument
        draw_directions :   direction argument
        out_file    :       output directory
        highlighted_genes:  set of genes to highlight in drawing
        label   :           label for figure
        absolute-pixel-size: absolute-pixel-size argument

    Writes
    ------
        Figure(s) to png image file(s)
    '''
    
    # TODO: there's an open issue of what to do when the label isn't provided.
    # Consider all usages and do it right.
    
    ignored_regions = organism.ignored_regions[chromosome]
    
    region_len = region_end - region_start
    
    track_height = 5
    feature_height = 20
#     label_height = 20 if label else 0
    label_height = 20
    
    width = 350 if absolute_pixel_size <= 0 else int(region_len / absolute_pixel_size)
    height = track_height * len(hits) + feature_height + label_height + track_height * 2
    
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context(surface)
    
    ctx.scale(width, height) # Normalizing the canvas
    
    # Draw white background
    ctx.rectangle(0, 0, 1, 1)
    ctx.set_source_rgb(1, 1, 1)
    ctx.fill()
    
    # Draw the tracks
    track_height_scaled = float(track_height) / height
    for track_ix, track_hits in enumerate(hits):
        track_y = track_ix * track_height / float(height)
        chrom_track = track_hits[chromosome]
        try:
            left_hit_ix = chrom_track.index(chrom_track.find_ge(region_start))
        except ValueError:
            left_hit_ix = 0
        try:
            right_hit_ix = chrom_track.index(chrom_track.find_gt(region_end))
        except ValueError:
            right_hit_ix = len(chrom_track)
        relevant_hits = chrom_track[left_hit_ix:right_hit_ix]
        
        for pixel in range(width):
            left_bp = region_start + int(ceil(pixel * (region_len / float(width))))
            right_bp = region_start + int(floor((pixel+1) * (region_len / float(width))))
            reads_in_pixel = sum(h["hit_count"] for h in relevant_hits if left_bp <= h["hit_pos"] <= right_bp)
            if reads_in_pixel:
                ctx.rectangle((pixel / float(width)),
                              track_y,
                              1.0 / width,
                              track_height_scaled)
                color = max(0, 0.9 - reads_in_pixel / 100.0)
                ctx.set_source_rgb(color, color, color)
                ctx.fill()
    
    # Draw the feature:
    features_in_range = organism.feature_db.get_features_at_range(chromosome, (region_start, region_end))
    
    feature_track_y = float(height - feature_height - label_height - track_height*2) / height
    feature_track_height_scaled = float(feature_height) / height
    scale_bp = lambda bp: min(1.0, max(0.0, float(bp - region_start)) / region_len) 
    for feature in features_in_range:
        # Draw the feature as a blue rectangle:
        feature_start = scale_bp(feature.start)
        feature_end = scale_bp(feature.stop)
        
        ctx.rectangle(feature_start, feature_track_y, feature_end - feature_start, feature_track_height_scaled)
        ctx.set_source_rgb(0.0, 100/255.0, 180/255.0)
        ctx.fill()
        
        # Draw introns:
        intron_track_y = feature_track_y + feature_track_height_scaled
        for intron_start, intron_end in feature.exons.complement(feature.start, feature.stop):
            intron_start = scale_bp(intron_start)
            intron_end = scale_bp(intron_end)
            ctx.rectangle(intron_start, intron_track_y, intron_end - intron_start, track_height_scaled)
            ctx.set_source_rgb(1, 1, 0)
            ctx.fill()
        
        # Draw its domains:
        if (feature.standard_name in highlighted_genes and draw_domains == GENES_HIGHLIGHTED) or \
            draw_domains == GENES_ALL:
            domain_pad = 1.0/6 * feature_height
            domain_y = feature_track_y + domain_pad / height
            domain_height = 2.0/3 * feature_track_height_scaled
            for domain in feature.domains:
                domain_start = scale_bp(domain[0])
                domain_end = scale_bp(domain[1])
                
                ctx.rectangle(domain_start, domain_y,
                              domain_end - domain_start, domain_height)
                ctx.set_source_rgb(0, 0, 0)
                ctx.fill()
                
        # Chop off a rectangle to indicate directionality:
        if (feature.standard_name in highlighted_genes and draw_directions == GENES_HIGHLIGHTED) or \
            draw_directions == GENES_ALL:
            arrow_width = min(10.0 / width, feature_end - feature_start) # 10 pixels or the length of the gene
            if feature.strand == 'W':
                end = scale_bp(feature.stop)
                ctx.move_to(end - arrow_width, feature_track_y)
                ctx.line_to(end, feature_track_y)
                ctx.line_to(end, feature_track_y + 0.5 * feature_track_height_scaled)
                ctx.close_path()
                ctx.set_source_rgb(1, 1, 1)
                ctx.fill()
                
                ctx.move_to(end - arrow_width, feature_track_y + feature_track_height_scaled)
                ctx.line_to(end, feature_track_y + feature_track_height_scaled)
                ctx.line_to(end, feature_track_y + 0.5 * feature_track_height_scaled)
                ctx.close_path()
                ctx.set_source_rgb(1, 1, 1)
                ctx.fill()
            else:
                start = scale_bp(feature.start)
                ctx.move_to(start, feature_track_y)
                ctx.line_to(start + arrow_width, feature_track_y)
                ctx.line_to(start, feature_track_y + 0.5 * feature_track_height_scaled)
                ctx.close_path()
                ctx.set_source_rgb(1, 1, 1)
                ctx.fill()
                
                ctx.move_to(start, feature_track_y + feature_track_height_scaled)
                ctx.line_to(start + arrow_width, feature_track_y + feature_track_height_scaled)
                ctx.line_to(start, feature_track_y + 0.5 * feature_track_height_scaled)
                ctx.close_path()
                ctx.set_source_rgb(1, 1, 1)
                ctx.fill()
                
        if label is None:# and feature.standard_name in highlighted_genes:
            feature_name = feature.name
            ctx.set_font_matrix(cairo.Matrix(xx=15/float(width), yy=15/float(height)))
            (_x, _y, label_text_width, label_text_height, _dx, _dy) = ctx.text_extents(feature_name)
            label_track_y = (height - label_height) / float(height)
            ctx.move_to(feature_start + (feature_end - feature_start - label_text_width) / 2, label_track_y + label_text_height * 1.25)
            ctx.set_source_rgb(0, 0, 0)
            ctx.show_text(feature_name)
    
    # Draw ignored regions:
    ignored_track_y = feature_track_y + feature_track_height_scaled + track_height_scaled
    for ignored_start, ignored_stop in ignored_regions & RangeSet([(region_start, region_end)]):
        ignored_start = scale_bp(ignored_start)
        ignored_stop = scale_bp(ignored_stop)
        ctx.rectangle(ignored_start, ignored_track_y, ignored_stop - ignored_start, track_height_scaled)
        ctx.set_source_rgb(1, 0, 0)
        ctx.fill()
    
    # Draw the text:
    if label:
        ctx.set_font_matrix(cairo.Matrix(xx=15/float(width), yy=15/float(height)))
        (_x, _y, label_text_width, label_text_height, _dx, _dy) = ctx.text_extents(label)
        label_track_y = (height - label_height) / float(height)
        ctx.move_to((1.0 - label_text_width) / 2.0, label_track_y + label_text_height * 1.25)
        ctx.set_source_rgb(0, 0, 0)
        ctx.show_text(label)
    
    # Output to PNG
    surface.write_to_png(out_file)


if __name__ == "__main__":
        main()