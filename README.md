# pyfastg-oop
A minimal Python library for parsing SPAdes FASTG files

Object oriented  programming version of pyfastg (https://github.com/fedarko/pyfastg)

A demo and the validation is included in the pyfast-oop_demo.ipynb jupyter notebook

### Testing:
I tested using SPAdes assembly graphs from RNAseq data from m. musculus and c. elegans<br>
Unfortunately, due to the file size limitation, I was only able to upload a compressed<br>
version of the c. elegans assembly graph.

### Extending pyfastg-oop
The object oriented approach of pyfastg-oop makes it easy to extend.
I already added an additional reader for graphs from BCALM, allowing to convert .fa-files
from BCALM to .fastg-files
