# Vicinity Analyzer
created 2023 by gaenssle
written in Python 3.8

## Process 
- It accepts three types of inputs
  * A KEGG Orthology (KO) ID (e.g. K22276)
  * A (list of) KEGG gene ID(s) (e.g. cak:Caul_3276), list with ',' and no spaces
  * A file (table) with gene IDs (.txt or .csv), with the column header ID'
- If the input is a KO ID, all associated gene IDs are downloaded first from KEGG
- Determine the gene label increments (1,5 or 10) for each corresponding KEGG genome
- Download all neighbouring genes within the given range (default= +/-5)
- Count the occurence of the provided input targets
  * KO ID
  * Pfam domain
  * Keyword in assigned name
  * May be prodived as file with target-type pairs (type=[KO-ID, Domain, Name])
- Export accumulated neighbours and occurence count

***

## Data
- Downloaded sequence data:
  * Assigned description
  * Organism and taxonomy
  * Sequence
  * Domain architecture
  * Available IDs from other databases   
- Count each target:
 * Per entry (occurences for each gene ID)
 * Per range (occurences at each range position)
- Save data in the following files:
  * Gene IDs (provided or downloaded via entered KO-ID)
  * Gene details for each neighbor (organism, architecture, sequence, etc)
  * Count files (entry and range)

***

## Dependencies

- The program used python 3.8 and the following modules:
  * pandas
  * argparse
  * Bio (Bio.KEGG)
  * ssl, urllib.request
  * os, re, io

***

## How to use

```
Main.py [-h] [-ask]
        [-ti TARGETID] [-td TARGETDOMAIN] [-tn TARGETNAME]
        [-tf TARGETFILE] [-a ACTION] [-r RANGE]
        [-n NAME] [-f FOLDER]
        [-cs CLUSTERSIZE] [-ft FILETYPE] [-sep SEPARATOR]
        input

VICINITY ANALYZER This program downloads neighboring genes from KEGG genomes
via KEGG.API The input is either a KO-ID or a list of gene IDs (in file or as
list) All gene IDs within the given range of the provided ID(s) are obtained
from KEGG It then downloads relevant details for each gene ID, e.g. organism
and domain architecture The occurrences of each provided target is counted per
entry and per position

positional arguments:
  input                 KO ID, KEGG gene ID(s) or file containing KEGG gene
                        IDs (e.g. blb:BBMN68_1454,blf:BLIF_1909)

optional arguments:
  -h, --help            show this help message and exit
  -ask, --askoverwrite  ask before overwriting files
  -ti TARGETID, --targetID TARGETID
                        target KO ID(s) used for filtering (e.g. K21572)
  -td TARGETDOMAIN, --targetDomain TARGETDOMAIN
                        target domain name(s) (Pfam) used for filtering (e.g.
                        RagB,SusD-like)
  -tn TARGETNAME, --targetName TARGETNAME
                        target keword(s) in name/annotation used for filtering
                        (e.g. RagB,SusD)
  -tf TARGETFILE, --targetFile TARGETFILE
                        File containing target;type pairs used for filtering
                        (types=[KO-ID, Name, Domain], sep=-sep)
  -a ACTION, --action ACTION
                        add actions to be conducted: a=all, i=retrieve IDs,
                        g=get neighbors, f=filter with target(default: a)
  -r RANGE, --range RANGE
                        +/- range in which genes will be searched (default: 5)
  -n NAME, --name NAME  name of files (default: same as 'input')
  -f FOLDER, --folder FOLDER
                        name of the parent folder (default: same as 'input')
  -cs CLUSTERSIZE, --clustersize CLUSTERSIZE
                        entries/frament files (default: 25)
  -ft FILETYPE, --filetype FILETYPE
                        type of the generated files (default: .csv)
  -sep SEPARATOR, --separator SEPARATOR
                        separator between columns in the output files
                        (default: ;)
```
