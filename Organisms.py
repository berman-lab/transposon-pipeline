import GenomicFeatures
from RangeSet import RangeSet
import csv
import Shared
import os
import pandas as pd

# Feature DB
# Ignored regions
# Deleted regions
# Homologous regions
# Literature essentials/non-essentials
# Ignored features at coverage 
class Organism(object):
    _KEYS = ("feature_db", "homologous_regions", "deleted_regions",
             "ignored_regions", "ignored_features", "literature_essentials",
             "literature_non_essentials", "genes_with_paralogs")
    
    def __init__(self):
        self._attr_cache = {}
        
        self._ignore_region_threshold = 50
        self._min_feature_coverage = 0.95
        
        # These are set here for PyDev's auto-completion:
        # TODO: if we do it like this, __getattr__ won't get called.
        # Switch to properties?
        self.feature_db = None
        self.homologous_regions = None
        self.deleted_regions = None
        self.ignored_regions = None
        self.ignored_features = None
        self.literature_essentials = None
        self.literature_non_essentials = None
        self.genes_with_paralogs = None
    
    @Shared.memoized
    def _read_range_data(self, file_name):
        result = {}
        
        with open(file_name, 'r') as in_file:
            reader = csv.reader(in_file)
            next(reader)
            for chrom, start, stop in reader:
                chrom = self.feature_db.get_std_chrom_name(chrom)
                if chrom not in result:
                    result[chrom] = []
                result[chrom].append((int(start), int(stop)))
            
        return {chrom: RangeSet(ranges) for chrom, ranges in result.items()}
    
    def __getattribute__(self, attr):
        if attr in Organism._KEYS:
            if attr not in self._attr_cache:
                self._attr_cache[attr] = getattr(self, "_get_%s" % attr)()
            return self._attr_cache[attr]
        else:
            return super(Organism, self).__getattribute__(attr)
        
    def _get_ignored_regions(self):
        deleted_regions = self.deleted_regions
        homologous_regions = self.homologous_regions
        result = {}
        for chrom in set(deleted_regions.keys() + homologous_regions.keys()):
            result[chrom] = deleted_regions.get(chrom, RangeSet()) + homologous_regions.get(chrom, RangeSet())
        
        return result
    
    def _get_ignored_features(self):
        result = set()
    
        for f in self.feature_db.get_all_features():
            chrom_ignored = self.ignored_regions.get(f.chromosome, RangeSet())
            if (chrom_ignored & RangeSet([(f.start, f.stop)])).coverage / float(len(f)) > (1.0 - self._min_feature_coverage) or \
                not f.is_orf:# or "dubious" in f.type.lower():
                result.add(f.standard_name)
        
        return result
    
    def _get_genes_with_paralogs(self, paralog_filename):
        with open(paralog_filename, 'r') as in_file:
            features = (self.feature_db.get_feature_by_name(f.strip()) for f in in_file.readlines())
            return set(f.standard_name for f in features if f is not None) 

class Calbicans(Organism):
    def _get_feature_db(self):
        return GenomicFeatures.default_alb_db()
    
    def _get_homologous_regions(self):
        ranges = self._read_range_data(Shared.get_dependency(os.path.join("albicans", "homologous_regions.csv")))
    
        return {chrom: RangeSet(r for r in range_set if r[1] - r[0] + 1 >= self._ignore_region_threshold) for chrom, range_set in ranges.iteritems()}
    
    def _get_deleted_regions(self):
        ranges = self._read_range_data(Shared.get_dependency(os.path.join("albicans", "deleted_regions.csv")))
    
        return {chrom: RangeSet(r for r in range_set if r[1] - r[0] + 1 >= self._ignore_region_threshold) for chrom, range_set in ranges.iteritems()}
    
    def _get_literature_essentials(self):
        return None
    
    def _get_literature_non_essentials(self):
        return None
    
    def _get_genes_with_paralogs(self):
        return Organism._get_genes_with_paralogs(self, Shared.get_dependency(os.path.join("albicans", "hasParalogs_ca.txt")))

class Scerevisiae(Organism):
    def _get_feature_db(self):
        return GenomicFeatures.default_cer_db()
    
    def _get_homologous_regions(self):
        ranges = self._read_range_data(Shared.get_dependency(os.path.join("cerevisiae", "homologous_regions.csv")))
    
        return {chrom: RangeSet(r for r in range_set if r[1] - r[0] + 1 >= self._ignore_region_threshold) for chrom, range_set in ranges.iteritems()}
    
    def _get_deleted_regions(self):
        # TODO: data unavailable, but there's a reason to suspect there are
        # deleted regions - during analysis, some features have no hits in
        # their neighborhoods.
        return {}
    
    def _get_literature_essentials(self):
        annotated_as_viable, annotated_as_inviable = self.get_nominal_annotations()
        return annotated_as_viable - annotated_as_inviable
    
    def _get_literature_non_essentials(self):
        annotated_as_viable, annotated_as_inviable = self.get_nominal_annotations()
        return annotated_as_inviable - annotated_as_viable
    
    @Shared.memoized
    def get_nominal_annotations(self):
        """Get nominal essentials and non-essentails in cerevisiae.
        
        Returns
        -------
        (set, set)
            A pair sets, denoting the essential and non-essential genes, using
            their standard names.
        """
        
        viable_filepath = Shared.get_dependency("cerevisiae", "cerevisiae_viable_annotations.txt")
        inviable_filepath = Shared.get_dependency("cerevisiae", "cerevisiae_inviable_annotations.txt")
        
        viable_table = pd.read_csv(viable_filepath, skiprows=8, delimiter="\t")
        inviable_table = pd.read_csv(inviable_filepath, skiprows=8, delimiter="\t")
        
        annotated_as_viable = set(self.feature_db.get_feature_by_name(f) for f in set(viable_table[viable_table["Mutant Information"] == "null"]["Gene"])) - set([None])
        annotated_as_inviable = set(self.feature_db.get_feature_by_name(f) for f in set(inviable_table[inviable_table["Mutant Information"] == "null"]["Gene"])) - set([None])
        
        # TODO: the dubious genes shouldn't be filtered here.
        consensus_viable_orfs = [f for f in annotated_as_viable if f.is_orf and f.feature_qualifier != "Dubious"]
        consensus_inviable_orfs = [f for f in annotated_as_inviable if f.is_orf and f.feature_qualifier != "Dubious"]
        
        return (set(f.standard_name for f in consensus_inviable_orfs),
                set(f.standard_name for f in consensus_viable_orfs))
        
    @property
    @Shared.memoized
    def conflicting_essentials(self):
        return self.literature_essentials & self.literature_non_essentials
    
    def _get_genes_with_paralogs(self):
        return Organism._get_genes_with_paralogs(self, Shared.get_dependency(os.path.join("cerevisiae", "hasParalogs_sc.txt")))
    
class Spombe(Organism):
    def _get_feature_db(self):
        return GenomicFeatures.default_pom_db()
    
    def _get_homologous_regions(self):
        ranges = self._read_range_data(Shared.get_dependency(os.path.join("pombe", "homologous_regions.csv")))
    
        return {chrom: RangeSet(r for r in range_set if r[1] - r[0] + 1 >= self._ignore_region_threshold) for chrom, range_set in ranges.iteritems()}
    
    def _get_deleted_regions(self):
        return {}
    
    def _get_literature_essentials(self):
        return self._get_spom_essentials()[0]
    
    def _get_literature_non_essentials(self):
        return self._get_spom_essentials()[1]

    @Shared.memoized
    def _get_spom_essentials(self):
        viability_table = pd.read_csv(Shared.get_dependency("pombe/FYPOviability.tsv"),
                                      header=None,
                                      delimiter='\t',
                                      names=["pombe standard name", "essentiality"])
        
        return set(r[0] for _ix, r in viability_table.iterrows() if r[1] == "inviable"), \
            set(r[0] for _ix, r in viability_table.iterrows() if r[1] == "viable")
            
    def _get_genes_with_paralogs(self):
        return Organism._get_genes_with_paralogs(self, Shared.get_dependency(os.path.join("pombe", "hasParalogs_sp.txt")))
        
# The lazy singletons:    
alb = Calbicans()
cer = Scerevisiae()
pom = Spombe()

@Shared.memoized
def get_genes_by_name(name):
    return (alb.feature_db.get_feature_by_name(name), cer.feature_db.get_feature_by_name(name), pom.feature_db.get_feature_by_name(name))

@Shared.memoized
def get_orths_by_name(name):
    alb_feature, cer_feature, pom_feature = get_genes_by_name(name)
    
    # If the name was not a Calb name, attempt to recover it:
    if alb_feature is None:
        alb_name = None
        if cer_feature:
            alb_name = cer_to_alb_map().get(cer_feature.standard_name)
        if not alb_name and pom_feature:
            alb_name = pom_to_alb_map().get(pom_feature.standard_name)
            
        if not alb_name:
            return (alb_feature, cer_feature, pom_feature)
        else:
            alb_feature = alb.feature_db.get_feature_by_name(alb_name)
    
    cer_ortholog = None
    if alb_feature.cerevisiae_orthologs:
        cer_ortholog = cer.feature_db.get_feature_by_name(list(alb_feature.cerevisiae_orthologs)[0])
    
    pom_ortholog = pom.feature_db.get_feature_by_name(get_calb_orths_in_sp().get(alb_feature.standard_name, ""))
    
    return (alb_feature, cer_ortholog, pom_ortholog)

@Shared.memoized
def cer_to_alb_map():
    result = {}
    
    for f in alb.feature_db.get_all_features():
        if not f.cerevisiae_orthologs:
            continue
        
        cer_orth = list(f.cerevisiae_orthologs)[0]
        cer_feature = cer.feature_db.get_feature_by_name(cer_orth)
        if cer_feature is None:
            print "WARNING: albicans ortholog %s doesn't exist in cerevisiae database!" % f.standard_name
            continue
        
        result[cer_feature.standard_name] = f.standard_name
        
    return result

@Shared.memoized
def pom_to_alb_map():
    return {pom.feature_db.get_feature_by_name(pom_name).standard_name: alb_name for alb_name, pom_name in get_calb_orths_in_sp().iteritems()}

@Shared.memoized
def get_calb_orths_in_sp():
    pom_db = GenomicFeatures.default_pom_db()
    
    ortholog_table = pd.read_csv(Shared.get_dependency("albicans", "C_albicans_SC5314_S_pombe_orthologs.txt"),
                                 skiprows=8,
                                 delimiter='\t',
                                 header=None,
                                 usecols=['albicans standard name', 'pombe standard name'],
                                 names=['albicans standard name', 'albicans common name', 'albicans alb_db id',
                                        'pombe standard name', 'pombe common name', 'pombe alb_db id'])
    
    # TODO: we probably don't want to use the hit table, though the InParanoid
    # table is very stringent.
    best_hit_table = pd.read_csv(Shared.get_dependency("albicans", "C_albicans_SC5314_S_pombe_best_hits.txt"),
                                 skiprows=8,
                                 delimiter='\t',
                                 header=None,
                                 usecols=['albicans standard name', 'pombe standard name'],
                                 names=['albicans standard name', 'albicans common name', 'albicans alb_db id',
                                        'pombe standard name', 'pombe common name', 'pombe alb_db id'])
    
    joined_table = pd.concat([ortholog_table, best_hit_table])
    
    result = {}
    for alb_feature in GenomicFeatures.default_alb_db().get_all_features():
        ortholog_row = joined_table[joined_table["albicans standard name"] == alb_feature.standard_name]
        if ortholog_row.empty:
            continue
        
        pom_feature = pom_db.get_feature_by_name(ortholog_row["pombe standard name"].iloc[0])
        if pom_feature:
            result[alb_feature.standard_name] = pom_feature.name
    
    return result

@Shared.memoized
def get_all_orthologs():
    result = []
    
    for f in alb.feature_db.get_all_features():
        orths = get_orths_by_name(f.standard_name)
        if None not in orths:
            result.append(orths)
        
    return result

def write_ignored_genes_table(organism, out_file):
    columns = ["Standard name", "Common name", "Type", "Deleted fraction", "Duplicated fraction", "Reason for exclusion"]
    ignored_features = [organism.feature_db.get_feature_by_name(n) for n in sorted(organism.ignored_features)]
    data = [
        [f.standard_name,
         f.common_name,
         f.type,
         (f.exons & organism.deleted_regions[f.chromosome]).coverage / float(len(f)),
         (f.exons & organism.homologous_regions[f.chromosome]).coverage / float(len(f))]
        for f in ignored_features
    ]
    
    for row in data:
        ignore_reasons = []
        
        if "ORF" not in row[3-1].upper():
            ignore_reasons.append("Not ORF")
#         if "dubious" in f.type.lower():
#             ignore_reasons.append("Dubious ORF")
        if row[-2] > 0.05:
            ignore_reasons.append("More than 5% deleted")
        if row[-1] > 0.05:
            ignore_reasons.append("More than 5% duplicated")
        
        row.append("; ".join(ignore_reasons))
    
    with pd.ExcelWriter(out_file) as excel_writer:
        result = pd.DataFrame(columns=columns, data=data)
        result.to_excel(excel_writer, sheet_name="Ignored genes", index=False)

if __name__ == "__main__":
    write_ignored_genes_table(alb, os.path.join(Shared.get_script_dir(), "output", "ignored_alb_genes.xlsx"))