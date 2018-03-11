# circM

## Introduction

The circM.py script for merging circRNA prediction output. 



## Usage

**circM.py**
 
``` python
Usage:
  circM.py [ARGUMENTS] > [OUTPUT]
  
Options:
  -f <files>        Input files
  -a <string>       Algorithm 
                    [acsf, circexplorer, circexplorer2, 
                    ciri, ciri2, circrna_finder, dcc 
                    find_circ, knife, mapsplice, uroborus] 
  -cf <int>         Expression threshold for each sample, back-splice spanning reads (default=2)
  -cs <int>         Merged expression threshold, backsplice-spanning reads (default=2)
```

## Example

To merge output from find_circ:

```bash
python circM.py -f *.circ_candidates.bed -a find_circ
```




## Citation

**Hansen TB. Improved circRNA Identification by Combining Prediction Algorithms. Front. Cell Dev Biol. 2018**


## License

Copyright (C) 2018 ncRNALab.  
