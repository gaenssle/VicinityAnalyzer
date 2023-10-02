#!/usr/bin/python
# Written in Python 3.8 in 2023 by A.L.O. Gaenssle

import os
import re
from multiprocessing import Pool
import pandas as pd
import argparse

# Own modules
import Import_Export as IE
import Download_KEGG as KEGG


## ------------------------------------------------------------------------------------------------
## INPUT ARGUMENTS --------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="VICINITY ANALYZER"
	"\nThis program downloads neighboring genes from KEGG genomes via KEGG.API"
    "\nThe input is an id (single or list)"
    "\nAll gene IDs within the given range of the provided ID(s) are obtained from KEGG"
    "\nIt then downloads relevant details for each gene ID, e.g. organism and domain architecture"
    "\nThe data can be counted regarding e.g. taxonomy, domain architecture and sequence length")
parser.add_argument("input", 
	help="KO ID, KEGG gene ID or file containing KEGG gene IDs")
parser.add_argument("-ask", "--askoverwrite", 
	help="ask before overwriting files",
	action="store_true")
parser.add_argument("-ti", "--targetID", 
	help="target KO ID used for filtering (e.g. K21572)")
parser.add_argument("-td", "--targetDomain", 
	help="target domain(s) (Pfam) used for filtering (e.g. RagB,SusD-like)")
parser.add_argument("-tk", "--targetKeyword", 
	help="target keword(s) in annotation used for filtering (e.g. RagB,SusD)")
parser.add_argument("-a", "--action", 
	help="add actions to be conducted: "
	"a=all, i=retrieve IDs, g=get neighbors, f=filter with target(default: %(default)s)", 
	default="a")
parser.add_argument("-r", "--range", 
	help="+/- range in which genes will be searched (default: %(default)s)", 
	default=5, 
	type=int)
parser.add_argument("-f", "--folder", 
	help="name of the parent folder (default: same as 'name')")
parser.add_argument("-ft", "--filetype", 
	help="type of the generated files (default: %(default)s)", 
	default=".csv")
parser.add_argument("-sep", "--separator", 
	help="separator between columns in the output files (default: %(default)s)", 
	default=";")

# Set folder name to searched ID if not set
args = parser.parse_args()
if args.folder == None:
    args.folder = args.input

# Check if the given action is valid and replace with the list if == 'a'
while all(ch in "aigf" for ch in args.action) == False:
    args.action = input("\nWhich action do you want to conduct?"
        "\n- a\tconduct all actions\n- i\tdownload sequence IDs"
        "\n- g\tget neighbors\n- f\tfilter with target\n"
        "\nPlease enter any or multiple of letters (e.g 'a' or 'igf' [without ''])\n")
if "a" in args.action:
    args.action = "igf"


## ------------------------------------------------------------------------------------------------
## MAIN FUNCTIONS ---------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------
## ================================================================================================
## Get index list of neighbors and retrieves protein data
def GetNeighbors(GeneList, FilePath, FileType, Sep, Ask):
	ManualList = []
	Complete = []
	Incomplete = []
	for GeneID in GeneList:
		ProteinList = []
		if GeneID.split(":",1)[1].isdigit() == False:
			try:
				List, ManualList = KEGG.GetNeighborIndices(GeneID, ManualList)
				Data = KEGG.DownloadProteinEntries(List, GeneID)
				if len(Data) < 6:
					List, ManualList = KEGG.GetNeighborIndices(GeneID, ManualList, Multi=5)
					Data = KEGG.DownloadProteinEntries(List, GeneID)
			except:
				try:
					List, ManualList = KEGG.GetNeighborIndices(GeneID, ManualList, Multi=5)
					Data = KEGG.DownloadProteinEntries(List, GeneID)
					if len(Data) < 6:
						List, ManualList = KEGG.GetNeighborIndices(GeneID, ManualList, Multi=10)
						Data = KEGG.DownloadProteinEntries(List, GeneID)
				except:
					ManualList.append(GeneID)
			for Entry in Data:
				ProteinList.append(KEGG.GetDetailedData(Entry, GeneID, GeneID.split(":",1)[0]))
			if len(Data) == 10:
				Complete.extend(ProteinList)
			else:
				Incomplete.extend(ProteinList)
		else:
			ManualList.append(GeneID)
	print("Neighbors found:", len(Complete), "complete,", len(Incomplete), "incomplete; of", len(GeneList))
	print("Manual search required for:", len(ManualList))
	if Complete != []:
		DataFrame = pd.DataFrame(Complete)
		IE.ExportDataFrame(DataFrame, FilePath + "_Neighbors", 
			FileType=FileType, Sep=Sep, Ask=Ask)
	if Incomplete != []:
		DataFrame = pd.DataFrame(Incomplete)
		IE.ExportDataFrame(DataFrame, FilePath + "_Neighbors_incomplete", 
			FileType=FileType, Sep=Sep, Ask=Ask)
	if ManualList != []:
		DataFrame = pd.DataFrame(ManualList, columns=["ID"])
		IE.ExportDataFrame(DataFrame, FilePath + "_ManualDownload", 
			FileType=FileType, Sep=Sep, Ask=Ask)



## ------------------------------------------------------------------------------------------------
## SCRIPT -----------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------

IE.CreateFolder(os.path.join(args.folder, "Input"))
IE.CreateFolder(os.path.join(args.folder, "Output"))


# Retrieve Sequence IDs from input, file or KEGG
if any(s in ["i", "g"] for s in args.action):
	OutputPath = os.path.join(args.folder, "Input", args.folder)
	if re.search(r"^K\d+$", args.input):
		print("KO ID provided as input")
		IDList, DataFrame = KEGG.DownloadOrthology(args.input)
	elif os.path.exists(args.input):
		print("File provided as input")
		DataFrame = pd.read_csv(args.input, sep=args.separator)
		IDList = DataFrame[DataFrame.columns[0]].to_list()
	elif ":" in args.input:
		print("Gene ID provided as input")
		if "," in args.input:
			IDList = args.input.split(",")
		else:
			IDList = [args.input]
		DataFrame = pd.DataFrame(IDList, columns=["ID"])
	else:
		print("Please provide a valid input"
			"\n- a KEGG Orthology (KO) ID (e.g. K22276)"
			"\n- KEGG gene ID(s) (e.g. cak:Caul_3276,stax:MC45_14985)"
			"\n- an existing file with gene IDs (.txt or .csv)")
		quit()
	if DataFrame.empty == False:
		IE.ExportDataFrame(DataFrame, OutputPath, 
			FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)

	# print(IDList)


# Find neighboring genes on KEGG
	if "g" in args.action:
		OutputPath = os.path.join(args.folder, "Output", args.folder + "_NeigbourIDs")
		FragmentFolder = IE.CreateFolder(OutputPath + "Fragments")
		FragmentFile = os.path.join(FragmentFolder, args.folder + "_NeigbourIDs")
		GetNeighbors(IDList, OutputPath, args.filetype, args.separator, args.askoverwrite)

