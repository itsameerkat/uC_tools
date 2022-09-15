

class PairsSnipper:
    
    def __init__(self, filepath, chromsizes, strands=None):
        self.filepath = filepath
#         self.columns = ['chrom1', 'chrom2', 'pos1', 'pos2', 'strand1', 'strand2']
        self.strands = strands
        self.chromsizes = chromsizes
    
    def select(self, region1, region2):
        print(region1, region2)
        chrom1, start1, end1 = region1
        chrom2, start2, end2 = region2
        
        columns = ['chrom1', 'chrom2', 'pos1', 'pos2', 'strand1', 'strand2', 'pair_type']
        pairq = SaneParquet(self.filepath, columns)
        df = pairq[chrom1, chrom2]
        if self.strands is not None:
            idx_s = (df.strand1 == self.strands[0]) & (df.strand2 == self.strands[1])
            df = df[idx_s]
            
        m, n = self.chromsizes[chrom1], self.chromsizes[chrom2]
        matrix = sp.csr_matrix((
            np.ones(len(df), dtype=int), 
            (df['pos1'].values, df['pos2'].values)
        ), (m, n))
        del df
        return matrix
    

    def snip(self, data, region1, region2, tup):

        s1, e1, s2, e2 = tup
        if s1 > s2:
            s1, s2 = s2, s1
            e1, e2 = e2, e1

        matrix = data
        m, n = matrix.shape
        dm, dn = e1 - s1, e2 - s2
        assert e1 >= 0
        assert e2 >= 0

        out_of_bounds = False
        pad_left = pad_right = pad_bottom = pad_top = None
        if s1 < 0:
            pad_bottom = -s1
            out_of_bounds = True
        if s2 < 0:
            pad_left = -s2
            out_of_bounds = True
        if e1 > m:
            pad_top = dm - (e1 - m)
            out_of_bounds = True
        if e2 > n:
            pad_right = dn - (e2 - n)
            out_of_bounds = True

        if out_of_bounds:        
            i0 = max(s1, 0)
            i1 = min(e1, m)
            j0 = max(s2, 0)
            j1 = min(e2, n)
            snippet = sp.csr_matrix((dm, dn))
        else:
            snippet = matrix[s1:e1, s2:e2]
        return snippet