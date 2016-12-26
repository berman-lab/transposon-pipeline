"""A module for wrapping genome feature tables.

Can work with multiple gene aliases, find features by location (and range),
compute the inter-feature region at a given range, and generally allow easier
access to the features.

Attributes
----------
cer_db : CerevisiaeFeatureDB
    The default cerevisiae feature database.

alb_db : AlbicansFeatureDB
    The default albicans feature database.

Notes
-----

All indices are 1-based.
"""

import Shared

class _FeatureDB(object):
    """The base class for genome feature databases.
    
    Wraps the tab-delimited tables, such as those used by both CGD and SGD.
    
    Attributes
    ----------
    _features : sequence of _Feature objects
        A sequence of all of the genome features.
    _name_map : dict of str to _Feature
        A dictionary mapping all known feature names (including aliases) to the
        feature instance itself. Names that map to more than one feature are
        ignored.
    _chrom_features : dict of str to list of _Feature objects
        A map of chromosome names to the features that reside on them. The
        features are sorted by their .start attribute.
    _max_feature_lengths : dict of str to int
        A map of the longest features in each chromosome. This gives a maximum
        search window when searching for features in `get_features_at_range`.
    _chrom_names : dict of obj to str
        A map of chromosome name aliases to their "canonical" names. Useful
        when the canonical names are very long, and we'd like to use common
        shorthands. Note that aliases can really be anything, for example
        integers as well as strings.
    """
    
    def __init__(self, features):
        self._features = tuple(features) # Read-only copy
        self._name_map = {}
        self._chrom_features = {}
        
        duplicate_names = set()
        for feature in features:
            for name in feature.all_names:
                # We we will only keep aliases that are unique to the feature:
                if name not in duplicate_names:
                    if name in self._name_map:
                        duplicate_names.add(name)
                        del self._name_map[name] 
                    else:
                        self._name_map[name] = feature
            
            chrom_name = feature.chromosome
            if chrom_name not in self._chrom_features:
                self._chrom_features[chrom_name] = []
            self._chrom_features[chrom_name].append(feature)
            
        for features in self._chrom_features.values():
            features.sort(key=lambda f: f.start)
            
        # Longest feature lengths are required for getting features in ranges:
        self._max_feature_lengths = {}
        for chrom, features in self._chrom_features.items():
            self._max_feature_lengths[chrom] = max(len(f) for f in features)
            
        # Create basis for chromosome name maps (subclasses can add their own
        # chromosome aliases):
        self._chrom_names = {}
        for chrom in self._chrom_features:
            self._chrom_names[chrom] = chrom
            
        self._create_chrom_cache()
    
    # Returns the features sorted by their .start attribute.
    def get_features_at_location(self, chrom, location):
        """
        Get all the features at a particular location in the genome. Can be
        more than one, as features can overlap.
        
        Parameters
        ----------
        chrom : str or int
            The target chromosome.
        location : int
            The target location in the chromosome. 1-based.
            
        Returns
        -------
        list of _Feature objects
            A list of features at a given location.
        """
        
        return self.get_features_at_range(chrom, (location, location))
    
    # Returns the features sorted by their .start attribute.
    def get_features_at_range(self, chrom, (start, stop)):
        """
        Get features that intersect a given range.
        
        Parameters
        ----------
        chrom : str or int
            The target chromosome.
        start : int
            The left border of the range. 1-based.
        stop : int
            The right border of the range. Inclusive. 1-based.
                
        Returns
        -------
        list of _Feature objects
            A list of the features (as _Feature instances), sorted by their
            .start attributes. 
        """
        chrom = self._chrom_names[chrom]
        features = self._chrom_features[chrom]
        # The left border is the furthest point at which we can find features
        # overlapping with the required range, which can deduce as we know the
        # what the length of the longest feature is.
        left_border = start - self._max_feature_lengths[chrom]
        
        last_feature_ix = _binary_search(features, stop, lambda f: f.start)
        first_feature_ix = last_feature_ix - 1
        
        result = []
        
        # Walk left:
        while first_feature_ix >= 0:
            f = features[first_feature_ix]
            if f.start < left_border:
                break
            if f.start <= stop and f.stop >= start:
                result.append(f)
            first_feature_ix -= 1
            
        # Walk right:
        while last_feature_ix < len(features):
            f = features[last_feature_ix]
            if f.start > stop:
                # Because the features are sorted by their .start attribute,
                # once we find a feature that starts after the requested range
                # stops, we can be sure that no succeeding feature will overlap
                # the range.
                break
            if f.start <= stop and f.stop >= start:
                result.append(f)
            last_feature_ix += 1
            
        result.sort(key=lambda f: f.start)
        
        return result

    def get_last_effective_chrom_index(self, chrom):
        """Return the last chromosome position that is accessible."""
        
        # TODO: we need real chromosome lengths.
        chrom = self._chrom_names[chrom]
        return max(f.stop for f in self._chrom_features[chrom])

    def get_feature_by_name(self, feature_name):
        return self._name_map.get(feature_name)

    def get_all_features(self):
        return self._features

    def get_interfeature_range(self, chrom_name, (start, stop)):
        """Return a collection of ranges within a given region that contain no
        features.
        
        The returned collection is maximal.
        
        Returns
        -------
        DisjointRange
            A special object that simplifies working with a union of disjoint intervals.
        """
        
        chrom = self[chrom_name]
        start = max(start, 1)
        stop = min(len(chrom), stop)
        features = chrom[start:stop]
        if len(features) == 0:
            return DisjointRange([(start, stop)])
        
        # First, compute the union of intervals of all features within the
        # given region. This is necessary because features can overlap. 
        joined_range = [[features[0].start, features[0].stop]]
        for f in features[1:]:
            _last_start, last_stop = last_range = joined_range[-1]
            # We know that f.start >= last_start.
            if f.start > last_stop + 1:
                # We need to add a new range
                joined_range.append([f.start, f.stop])
            elif f.stop > last_stop:
                # We only extend the previous one
                last_range[1] = f.stop

        # Get the inverse of the features:
        inverse_range = []
        if start < joined_range[0][0]:
            inverse_range.append((start, joined_range[0][0]-1))
        for i in range(len(joined_range)-1):
            current_stop, next_start = joined_range[i][1], joined_range[i+1][0]
            assert next_start - current_stop > 1
            inverse_range.append((current_stop+1, next_start-1))
        if stop > joined_range[-1][1]:
            inverse_range.append((joined_range[-1][1]+1, stop))

        return DisjointRange(inverse_range)

    def _create_chrom_cache(self):
        self._cached_chroms = {}
        for name in self._chrom_features:
            self._cached_chroms[name] = Chromosome(name, self)

    def __getitem__(self, key):
        # Use caching because creating the chromosome may trigger expensive
        # computation, e.g. compute the length.
        return self._cached_chroms.get(self._chrom_names.get(key))

    def __iter__(self):
        for chrom in self._cached_chroms.values():
            yield chrom

class AlbicansFeatureDB(_FeatureDB):
    A21, A22 = (21, 22)
    
    def __init__(self, feature_filename):
        features = []
        with open(feature_filename, "r") as in_file:
            for line in in_file:
                if not line:
                    continue
                if line.startswith('!'):
                    if "Genome version" in line:
                        if "A21" in line:
                            self.assembly = AlbicansFeatureDB.A21
                        elif "A22" in line:
                            self.assembly = AlbicansFeatureDB.A22
                    continue

                feature = AlbicansFeature(line)
                if ("chrM" in feature.chromosome or  # Ignore the mitochondria
                    feature.primary_name.endswith("_B") or  # If A22, only keep the A haplotype
                    feature.start == -1): # Unmapped features aren't useful
                    continue
                
                features.append(feature)
                
        super(AlbicansFeatureDB, self).__init__(features)
        
        # Enrich the chromosome name maps:
        if self.assembly == AlbicansFeatureDB.A21:
            chrom_template = "Ca21chr%s_C_albicans_SC5314"
        else:
            chrom_template = "Ca22chr%sA_C_albicans_SC5314"
        for chrom_shorthand in list(range(1, 8)) + ['R']:
            chrom_name = chrom_template % chrom_shorthand
            # We allow numeric and string identifiers for indexed chromosomes:
            self._chrom_names[chrom_shorthand] = chrom_name
            self._chrom_names[str(chrom_shorthand)] = chrom_name

class CerevisiaeFeatureDB(_FeatureDB):
    def __init__(self, feature_filename):
        features = []
        with open(feature_filename, 'r') as in_file:
            for line in in_file:
                feature = CerevisiaeFeature(line)
                
                if (feature.type in ("telomere", "ARS") or # Ignore unwanted feature types
                    feature.chromosome == '17' or # Ignore mitochondria
                    feature.start == -1 or # Unmapped features begone!
                    "chromosome" not in feature.parent_feature): # We don't care about sub-annotations
                    continue
                
                features.append(feature)
                
        super(CerevisiaeFeatureDB, self).__init__(features)
        
        self._chrom_names.update({"chrI" : "1", "chrII" : "2", "chrIII" : "3", "chrIV" : "4",
                                  "chrV" : "5", "chrVI" : "6", "chrVII" : "7", "chrVIII" : "8",
                                  "chrIX" : "9", "chrX" : "10", "chrXI" : "11", "chrXII" : "12",
                                  "chrXIII" : "13", "chrXIV" : "14", "chrXV" : "15", "chrXVI" : "16",})

class Chromosome(object):
    """A chromosome in the genome.
    
    Mostly for convenience, as it wraps _FeatureDB methods."""
    
    def __init__(self, name, alb_db):
        self._name = self.name = name
        self._db = alb_db
        self._len = alb_db.get_last_effective_chrom_index(name)

    def __len__(self):
        return self._len

    def __getitem__(self, key):
        if type(key) == slice:
            return self._db.get_features_at_range(self._name,
                                                  (key.start, key.stop))
        else:
            return self._db.get_features_at_location(self._name, key)
        
    def __iter__(self):
        return iter(self.get_features())
        
    def get_features(self):
        # TODO: optimize to use the DB's _chrom_features?
        return sorted((f for f in self._db.get_all_features() if f.chromosome == self._name),
                      key=lambda f: f.start)


class _Feature(object):
    """A base class for a genomic feature.
    
    Attributes
    ----------
    strand : {'W', 'C'}
        'W' stands for the Watson/forward/+ strand, 'C' stands for the
        Crick/reverse/- strand.
    start : int
        The 1-based leftmost position of the feature.
    stop : int
        The 1-based rightmost position of the feature.
    """
    
    def __init__(self, components, primary_name_col, common_name_col,
                 type_col, chrom_col, start_col, stop_col,
                 strand_col, description_col):
        self.primary_name = self.standard_name = components[primary_name_col]
        self.common_name = components[common_name_col]
        self.name = self.common_name or self.standard_name
        self.all_names = set((self.name, self.primary_name))
        
        self.type = components[type_col]
        self.chromosome = components[chrom_col]
        self.strand = components[strand_col]
        self.description = components[description_col]
        
        if components[start_col] == '':
            self.start = self.stop = -1
        else:
            start = int(components[start_col])
            stop = int(components[stop_col])
            self.start, self.stop = min(start, stop), max(start, stop)
            
        self._len = self.stop - self.start + 1
        
        self.is_orf = "ORF" in self.type

    def __len__(self):
        return self._len

class AlbicansFeature(_Feature):
    def __init__(self, feature_line):
        # TODO: the file should be read as a CSV file, so we shouldn't
        # have to split it ourselves.
        components = feature_line.split("\t")
        
        super(AlbicansFeature, self).__init__(components, 0, 1, 3, 4, 5, 6, 7, 10)
        
        self.primary_cgdid = components[9-1] # Example: CAL0000184345. Used in GO annotations.
        self.all_names |= set(filter(None, components[2].split('|')) + [self.primary_cgdid])
        self.cerevisiae_orthologs = \
            frozenset(c.strip() for c in components[17].split('|') if c.strip())

class CerevisiaeFeature(_Feature):
    def __init__(self, feature_line):
        components = feature_line.split("\t")
        
        super(CerevisiaeFeature, self).__init__(components, 1-1, 5-1, 2-1, 9-1, 10-1, 11-1, 12-1, 16-1)
        
        self.parent_feature = components[7-1] # This is for filtering, as we don't care about sub-annotations.
        
        # The "feature name" seems to be a legacy SGD identifier of the from YDRXXX[W/C] (for ORFs).
        self.feature_name = components[4-1]
        self.all_names.add(self.feature_name)
        self.name = self.common_name or self.feature_name or self.standard_name
        self.feature_qualifier = components[3-1]

# 
# 
def _binary_search(seq, query, key_func=lambda x: x):
    """Ye olde binary search.
    
    Returns an index mid at which seq[mid] == query (not necessarily the first
    one). If such an index doesn't exist, returns the first index at which
    seq[mid] > query. If such an index doesn't exist, returns the last index of
    seq.
    
    Parameters
    ----------
    seq : indexable sequence
        The sequence to search.
    query : object
        The value to look for.
    key_func : function
        A key funciton, similar to the `key` argument in `sorted`, that will be
        applied to each item in `seq` for comparison with `query`.
    """
    
    # TODO: consider replacing with bisect. We can preserve the key_func
    # functionality by having a class that wraps the sequence and calls
    # the key_func in __getitem__.
    
    # Search boundaries:
    low = 0
    high = len(seq) - 1
    while low < high:
        mid = (low + high) / 2
        cmp_result = cmp(key_func(seq[mid]), query)
        if cmp_result == 0:
            return mid
        elif cmp_result < 0:
            # seq[mid] < query
            if low != mid:
                low = mid
            else:
                low = mid + 1
        else:
            # query < seq[mid]
            if high != mid:
                high = mid
            else:
                high = mid - 1
    return low


class DisjointRange(object):
    """A disjoint union of closed intervals.
    
    Used to represent inter-feature regions.
    
    Attributes
    ----------
    coverage : int
        The sum of lengths of the intervals composing this object.
    """
    
    def __init__(self, ranges):
        """
        Parameters
        ----------
        ranges : sequence of (int, int)
            A sequence of closed intervals. Assumes they are already disjoint
            and sorted.
        """
        
        self._ranges = ranges
        self.start = ranges[0][0]
        self.stop = ranges[-1][1]
        self.coverage = sum(r[1] - r[0] + 1 for r in ranges)

    def __contains__(self, pos):
        ix = _binary_search(self._ranges, pos, lambda r: r[0])
        start, stop = self._ranges[ix]
        if pos >= start and pos <= stop:
            # Either pos is exactly equal to the start of an interval, or we
            # overshot the range, but are still contained within the last
            # interval.
            return True
        elif ix > 0 and pos <= self._ranges[ix-1][1]:
            # We know that pos < self._ranges[ix][0], so the pos could be
            # contained in the previous interval.
            return True
        
        return False

# A cache for the pombe genes, so we don't have to parse the file each time:
_POMBE_GENES = None
def get_pombe_genes():
    """Get a list of pombe genes.
    
    S. pombe doesn't have a feature table similar to albicans or cerevisiae,
    but it suffices to use a simple standard name -> common name map instead.
    
    Returns
    -------
    dict of str to str
        A mapping between the standard name and the common name.
    """
    
    global _POMBE_GENES
    if _POMBE_GENES is not None:
        return _POMBE_GENES
    
    import re
    
    pombe_gene_extractor = "(?:I|II|III)\tPomBase\tgene\t.*?ID=gene:(\S+?);.*Name=(\S+?);"
    with open(Shared.get_dependency("pombe/schizosaccharomyces_pombe.chr.gff3"), 'r') as pombe_gff_file:
        pombe_genes = re.findall(pombe_gene_extractor, pombe_gff_file.read())
        
    _POMBE_GENES = {standard_name: common_name for (standard_name, common_name) in pombe_genes}
    return _POMBE_GENES


# The default feature DBs, for convenience.
cer_db = CerevisiaeFeatureDB(Shared.get_dependency("cerevisiae/SGD_features.tab"))
alb_db = AlbicansFeatureDB(Shared.get_dependency("albicans/reference genome/C_albicans_SC5314_version_A22-s07-m01-r08_chromosomal_feature.tab"))