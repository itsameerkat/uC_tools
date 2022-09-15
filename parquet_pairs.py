import pyarrow as pa
import pyarrow.parquet



class SaneParquet:
    def __init__(self, filepath, columns=None, **kwargs):
        self.f = pa.parquet.ParquetFile(filepath, **kwargs)
        self.columns = columns
        self.row_groups = [self.f.metadata.row_group(i) 
                               for i in range(self.f.metadata.num_row_groups)]
        self.index = pd.DataFrame(([
                rg.column(1).statistics.min.strip(), 
                rg.column(3).statistics.min.strip(), 
                rg.column(2).statistics.min, 
                rg.column(2).statistics.max
            ] for rg in self.row_groups),
            columns=['chrom1', 'chrom2', 'pos1.min', 'pos1.max'])
        
        self.index['chrom1'] =self.index['chrom1'].str.decode('utf-8')
        self.index['chrom2'] =self.index['chrom2'].str.decode('utf-8')
        
        # row groups happen to be chromosome blocks in this case, so we can do this
        self.blocks = defaultdict(list)
        for ix, row in self.index.iterrows():
            self.blocks[row.chrom1, row.chrom2].append(ix)
            if row.chrom1 != row.chrom2:
                self.blocks[row.chrom2, row.chrom1].append(ix)
    
    def __getitem__(self, key):
        return pd.concat(
            (self.f.read_row_group(rg, columns=self.columns).to_pandas() 
                for rg in self.blocks[key]),
            axis=0)