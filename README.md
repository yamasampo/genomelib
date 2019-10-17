# genomelib

Retrieve data/information of gene, CDS and/or intron from an annotation table and genomic sequences. Allele frequency spectrum can also be stored so that you can filter genes based on evolutionary history (e.g. retrieve only gene which shows GC preference).

## Data Structure (`EvoGenDatabase`)

`EvoGenDatabase` class inherit basic functions for filtering or searching from `Database` class. In addition to table searching, `EvoGenDatabase` can store DNA sequences and add basic statistics such as GC content in table. In the future, SFS will also be processable in the class.

### Attributes

- df: `pd.DataFrame`
- description: `str`

