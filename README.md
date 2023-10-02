# Vicinity Analyzer
created 2023 by gaenssle
written in Python 3.8

## Process 
- It accepts three types of inputs
  * A KEGG Orthology (KO) ID (e.g. K22276)
  * A KEGG gene ID (e.g. cak:Caul_3276)
  * A file with gene IDs (.txt or .csv)
- If the input is a KO ID, all associated gene IDs are downloaded first from KEGG
- Determine if gene IDs of each included KEGG genome are in 1, 5, or 10 increments
- Download all neighbouring genes within the given range (default= +/-5)
- Count the occurence of an input target
  * KO ID
  * Pfam domain
  * Keyword in assigned name
- Export accumulated neighbours and occurence count

***

## Data
- Downloaded sequence data:
  * Assigned description
  * Organism and taxonomy
  * Sequence
  * Domain architecture
  * Available IDs from other databases   
- Count:
 * Taxonomy (Phylum and organism)
 * Domain architecture
- Save data in the following files:
  * Gene IDs
  * Gene details (organism, architecture, sequence, etc)
  * Summary domain architecture
  * Summary of taxonomic distribution
  * Fasta files (containing entire sequence)
  * Fasta files (only containing the target domain sequence)  

***

## Dependencies

- The program used python 3.8 and the following modules:
  * pandas
  * argparse
  * multiprocessing (optional for Linux)

***

## How to use

```
python3 Main.py name [-h] [-m] [-ask] 
               [-db DBLIST] [-a ACTION] [-st SEARCHTYPE]
               [-c CUTOFF] [-sam SAMPLESIZE] [-f FOLDER] [-cs CLUSTERSIZE]
               [-ft FILETYPE] [-sep SEPARATOR]


positional arguments:

  name                  name of the domain

optional arguments:
  -h, --help            show this help message and exit
  -m, --multiprocess    turn on mutltiprocessing (only for Linux)
  -ask, --askoverwrite  ask before overwriting files
  -db DBLIST, --dblist DBLIST
                        list databases to be searched, separated by ','
                        (default: UniProt;KEGG;PDB;swissprot)
  -a ACTION, --action ACTION
                        add actions to be conducted: a=all, i=entry IDs,
                        d=protein data, m=KEGG motif, e=extract (default: a)
  -st SEARCHTYPE, --searchtype SEARCHTYPE
                        type of the searched id (default: pf)
  -c CUTOFF, --cutoff CUTOFF
                        min E-Value of Pfam domains (default: 0.0001)
  -sam SAMPLESIZE, --samplesize SAMPLESIZE
                        max number of downloaded entries (default: 0)
  -f FOLDER, --folder FOLDER
                        name of the parent folder (default: same as 'name')
  -cs CLUSTERSIZE, --clustersize CLUSTERSIZE
                        entries/frament files (default: 100)
  -ft FILETYPE, --filetype FILETYPE
                        type of the generated files (default: .csv)
  -sep SEPARATOR, --separator SEPARATOR
                        separator between columns in the output files
                        (default: ;)
```
