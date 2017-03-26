from SortedCollection import SortedCollection
from operator import itemgetter

class RangeSet(object):
    # TODO: currently doesn't handle the null range set very well. Should
    # introduce a NULL static singelton somehow. 
    
    def __init__(self, ranges=tuple()):
        # Sort by the start of every range:
        self._ranges = SortedCollection(ranges, itemgetter(0))
        
        if ranges:
            self._consolidate()
        
            self.begin = self.start = self._ranges[0][0]
            self.end = self.stop = self._ranges[-1][1]
            
            self.span = self.end - self.begin + 1
            self.coverage = sum(end - begin + 1 for (begin, end) in self._ranges)
        else:
            self.begin = self.start = self.end = self.stop = None
            self.span = self.coverage = 0
    
    def __len__(self):
        return len(self._ranges)
    
    def __iter__(self):
        return iter(self._ranges)
    
    def __getitem__(self, key):
        return self._ranges[key]
    
    def __contains__(self, pos):
        try:
            begin, end = self._ranges.find_le(pos)
            return pos >= begin and pos <= end
        except ValueError:
            return False
        
    def __add__(self, other):
        return RangeSet(list(self) + list(other))
    
    def __or__(self, other):
        return self + other
    
    def __and__(self, other):
        leftmost = min(self.start, other.start)
        rightmost = max(self.stop, other.stop)
        
        return (self.complement(leftmost, rightmost) | other.complement(leftmost, rightmost)).complement(leftmost, rightmost)
    
    def __sub__(self, other):
        return self & other.complement(self.start, self.stop)
    
    def __eq__(self, other):
        return len(self) == len(other) and all(b1 == b2 and e1 == e2 for (b1, e1), (b2, e2) in zip(self, other))
    
    def __str__(self):
        return str(list(self._ranges))
    
    def _consolidate(self):
        new_ranges = SortedCollection(key=itemgetter(0))
        prev_begin, prev_end = self._ranges[0]
        for begin, end in self._ranges[1:]:
            if prev_end >= begin - 1:
                # Consolidate the previous and current ranges:
                prev_end = max(prev_end, end)
            else:
                # Add the previous range, and continue with the current range
                # as the seed for the next iteration:
                new_ranges.insert((prev_begin, prev_end))
                prev_begin = begin
                prev_end = end
                
        new_ranges.insert((prev_begin, prev_end))
        
        self._ranges = new_ranges
    
    def complement(self, begin, end):
        if not self:
            return RangeSet([(begin, end)])
        
        inter_complement = [(e1+1, b2-1) for (b1, e1), (b2, e2) in zip(self._ranges, self._ranges[1:])]
        
        if begin < self.start:
            inter_complement.append((begin, self.start - 1))
        if end > self.end:
            inter_complement.append((self.end + 1, end))
            
        return RangeSet(inter_complement)
    
    def intersects(self, begin, end):
        # TODO: can be optimized
        return bool(self & RangeSet([(begin, end)]))
    
    def cut_to(self, begin, end):
        # TODO: can be optimized
        return self & RangeSet([(begin, end)])
    
if __name__ == "__main__":
    # Ersatz testing:
    r1 = RangeSet([(10, 20)])
    r2 = RangeSet([(30, 40)])
    r3 = RangeSet([(70, 80), (50, 60)])
    
    r = r1 + r2 + r3
    rc = r.complement(1, 100)
    rc2 = rc - RangeSet([(1, 50)])
    
    assert r3.begin == r3.start == 50
    assert r3.end == r3.stop == 80
    assert r3.span == 31
    assert r3.coverage == 22
    assert list(r) == [(10, 20), (30, 40), (50, 60), (70, 80)]
    assert list(rc) == [(1, 9), (21, 29), (41, 49), (61, 69), (81, 100)]
    assert list(rc2) == [(61, 69), (81, 100)]
    assert rc2.begin == rc2.start == 61
    assert rc2.end == rc2.stop == 100
    assert rc2.span == 40
    assert rc2.coverage == 29
    assert 65 in rc2
    assert rc2.intersects(75, 81)
    assert 105 not in rc2
    assert not rc2.intersects(50, 60)
    
    r4 = RangeSet([(61, 61), (61, 69), (81, 90), (85, 95), (90, 100)])
    assert rc2 == r4
    
    null = RangeSet()
    assert list(null) == []
    assert list(null.complement(1, 10)) == [(1, 10)]