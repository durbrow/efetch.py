# efetch.py
A Python module for retrieving sequence data as FASTA from NCBI using eutils

Requires the ability to communicate with eutils.ncbi.nlm.nih.gov

The following example will fetch the sequence of chromosome 8 of GRCh37:
```bash
python efetch.py CM000670.1
```
