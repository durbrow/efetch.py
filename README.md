# efetch.py
A Python module for retrieving sequence data as FASTA from NCBI using eutils

Requires the ability to communicate with eutils.ncbi.nlm.nih.gov

The following example will fetch the sequence of chromosome 8 of GRCh37:
```bash
python efetch.py CM000670.1
```

This can also be imported as a module, like so:
```python
import efetch

defline, lines = efetch.FASTA('CM000670.1')
assert(defline == '>gi|224384761|gb|CM000670.1| Homo sapiens chromosome 8, GRCh37 primary reference assembly')
seqlen = 0
for i in lines:
    seqlen += len(i)
assert(seqlen == 146364022)
```
