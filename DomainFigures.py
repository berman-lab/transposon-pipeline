from math import floor, ceil
import cairocffi as cairo
import GenomicFeatures
import argparse
import SummaryTable
import glob
import os
from argparse import ArgumentParser
import shlex

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--hits-folder", default=".")
    parser.add_argument("--output-folder", default=".")

    subparsers = parser.add_subparsers(dest="source_type")
    
    region_parser = subparsers.add_parser("region")
    region_parser.add_argument("chromosome")
    region_parser.add_argument("start", type=int)
    region_parser.add_argument("stop", type=int)
    
    gene_list_parser = lambda gs: [g.upper() for g in gs.split(',')]
    
    gene_parser = subparsers.add_parser("gene")
    gene_name = gene_parser.add_mutually_exclusive_group(required=True)
    gene_name.add_argument("--alb", type=gene_list_parser)
    gene_name.add_argument("--cer", type=gene_list_parser)
    gene_region = gene_parser.add_mutually_exclusive_group()
    gene_region.add_argument("--percent-of-length", type=float, default=0.15)
    gene_region.add_argument("--bps", type=int)
    
    config_file_parser = subparsers.add_parser("config_file")
    config_file_parser.add_argument("config_file")
    
    args = parser.parse_args()
    
    hits = SummaryTable.read_hit_files(glob.glob(os.path.join(args.hits_folder, "*_Hits.txt")))
    
    if not hits:
        raise Exception("No hit files were found in the hits folder: %s" % args.hits_folder)
    
    if args.source_type == "config_file":
        with open(args.config_file, "r") as in_file:
            for line in in_file:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                handle_args(parser.parse_args(shlex.split(line)), hits)
    else:
        handle_args(args, hits)
    
def handle_args(args, hits):
    alb_db = GenomicFeatures.default_alb_db()
    
    if args.source_type == "gene":
        if args.alb is not None:
            features = map(alb_db.get_feature_by_name, args.alb)
            gene_names = args.alb
        else: # args.cer exists
            features = map(alb_db.get_feature_by_cerevisiae_ortholog, args.cer)
            gene_names = args.cer
            
        missed_genes = [name for f, name in zip(features, gene_names) if f is None]
        if missed_genes:
            print "The following genes weren't found and will be skipped: %s" % ', '.join(missed_genes)
            features = filter(None, features)
        
        for feature in features:
            if feature.cerevisiae_orthologs:
                cer_ortholog = list(feature.cerevisiae_orthologs)[0]
            else:
                cer_ortholog = ""
            
            # TODO: ommit Scer if no ortholog exists?
            if feature.name == feature.standard_name:
                name = "Calb-%s-Scer-%s" % (feature.name, cer_ortholog)
            else:
                name = "Calb-%s_%s-Scer-%s" % (feature.name, feature.standard_name, cer_ortholog)
                
            chromosome = feature.chromosome
            start = feature.start
            stop = feature.stop
            
            if args.bps is not None:
                start -= args.bps
                stop += args.bps
            else: # args.percent_of_length exists
                start -= args.percent_of_length * len(feature)
                stop += args.percent_of_length * len(feature)
                
            start = int(max(1, start))
            stop = int(min(stop, len(alb_db[chromosome])))
            
            draw_gene({"cerevisiae_name": cer_ortholog, "albicans_name": feature.name}, hits, args.output_folder)
#             print name, chromosome, start, stop
    else: # "region"
        chromosome = args.chromosome
        start = args.start
        stop = args.stop
        name = "Chr%s_%d-%d" % (chromosome, start, stop)
        
        draw_genomic_region(chromosome, start, stop, hits, args.output_folder)
#         print name, chromosome, start, stop
    
def draw_gene(config, hits, out_dir):
    alb_db = GenomicFeatures.default_alb_db()
    
    target_gene = alb_db.get_feature_by_cerevisiae_ortholog(config["cerevisiae_name"]) or \
        alb_db.get_feature_by_name(config["albicans_name"])
    gene_pad = 0.2 * len(target_gene)
    region_start, region_end = floor(target_gene.start - gene_pad), ceil(target_gene.stop + gene_pad)
    region_len = region_end - region_start
    
    track_height = 5
    feature_height = 20
    
    width = 350
    height = track_height * len(hits) + feature_height
    
    surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context (surface)
    
    ctx.scale (width, height) # Normalizing the canvas
    
    # Draw white background
    ctx.rectangle(0, 0, 1, 1)
    ctx.set_source_rgb(1, 1, 1)
    ctx.fill()
    
    # Draw the tracks
    pixel_width = 1.0 / width
    track_height_scaled = float(track_height) / height
    for track_ix, track_hits in enumerate(hits):
        track_y = track_ix * track_height / float(height)
        relevant_hits = [h for h in track_hits if h["chrom"] == target_gene.chromosome
                         and region_start <= h["hit_pos"] <= region_end]
        
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
    features_in_range = alb_db.get_features_at_range(target_gene.chromosome, (region_start, region_end))
    
    feature_track_y = float(height - feature_height) / height
    feature_track_height_scaled = float(feature_height) / height
    scale_bp = lambda bp: min(1.0, max(0.0, float(bp - region_start)) / region_len) 
    for feature in features_in_range:
        start = scale_bp(feature.start)
        end = scale_bp(feature.stop)
        
        ctx.rectangle(start, feature_track_y, end - start, feature_track_height_scaled)
        ctx.set_source_rgb(0.0, 100/255.0, 180/255.0)
        ctx.fill()
        
        # Draw its domains:
        if feature == target_gene:
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
            arrow_width = 0.05
            if feature.strand == 'W':
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
    
    # Output to PNG
    name = "sc-%s_ca-%s.png" % (config["cerevisiae_name"], target_gene.name)
    surface.write_to_png(os.path.join(out_dir, name))

def draw_genomic_region(chromosome, region_start, region_end, hits, out_dir):
    alb_db = GenomicFeatures.default_alb_db()
    
    region_len = region_end - region_start
    
    track_height = 5
    feature_height = 20
    
    width = 500
    height = track_height * len(hits) + feature_height
    
    surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context (surface)
    
    ctx.scale (width, height) # Normalizing the canvas
    
    # Draw white background
    ctx.rectangle(0, 0, 1, 1)
    ctx.set_source_rgb(1, 1, 1)
    ctx.fill()
    
    # Draw the tracks
    pixel_width = 1.0 / width
    track_height_scaled = float(track_height) / height
    for track_ix, track_hits in enumerate(hits):
        track_y = track_ix * track_height / float(height)
        relevant_hits = [h for h in track_hits if h["chrom"] == chromosome
                         and region_start <= h["hit_pos"] <= region_end]
        
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
    features_in_range = alb_db.get_features_at_range(chromosome, (region_start, region_end))
    
    feature_track_y = float(height - feature_height) / height
    feature_track_height_scaled = float(feature_height) / height
    scale_bp = lambda bp: min(1.0, max(0.0, float(bp - region_start)) / region_len) 
    for feature in features_in_range:
        start = scale_bp(feature.start)
        end = scale_bp(feature.stop)
        
        ctx.rectangle(start, feature_track_y, end - start, feature_track_height_scaled)
        ctx.set_source_rgb(0.0, 100/255.0, 180/255.0)
        ctx.fill()
    
    # Output to PNG
    name = "%s_%d-%d.png" % (chromosome, region_start, region_end)
    surface.write_to_png(os.path.join(out_dir, name))

if __name__ == "__main__":
    # The raw_input() at the end is useful when this is run within a batch file
    # in Windows, so that we can see the output before it disappears.
    try:
        main()
    except Exception:
        import traceback
        print traceback.format_exc()
        raw_input("An error occurred! Please tell Vlad.")
    else:
        raw_input("Finished! Press Enter to close the window.")
