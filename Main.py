#!/usr/bin/python
# Written in Python 3.8 in 2023 by A.L.O. Gaenssle

import os
import re
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
    "\nThe input is either a KO-ID or a list of gene IDs (in file or as list)"
    "\nAll gene IDs within the given range of the provided ID(s) are obtained from KEGG"
    "\nIt then downloads relevant details for each gene ID, e.g. organism and domain architecture"
    "\nThe occurrences of each provided target is counted per entry and per position")
parser.add_argument("input", 
	help="KO ID, KEGG gene ID(s) or file containing KEGG gene IDs (e.g. blb:BBMN68_1454,blf:BLIF_1909)")
parser.add_argument("-ask", "--askoverwrite", 
	help="ask before overwriting files",
	action="store_true")
parser.add_argument("-ti", "--targetID", 
	help="target KO ID(s) used for filtering (e.g. K21572)")
parser.add_argument("-td", "--targetDomain", 
	help="target domain name(s) (Pfam) used for filtering (e.g. RagB,SusD-like)")
parser.add_argument("-tn", "--targetName", 
	help="target keword(s) in name/annotation used for filtering (e.g. RagB,SusD)")
parser.add_argument("-tf", "--targetFile", 
	help="File containing target;type pairs used for filtering (types=[KO-ID, Name, Domain], sep=-sep)")
parser.add_argument("-a", "--action", 
	help="add actions to be conducted: "
	"a=all, i=retrieve IDs, g=get neighbors, f=filter with target(default: %(default)s)", 
	default="a")
parser.add_argument("-r", "--range", 
	help="+/- range in which genes will be searched (default: %(default)s)", 
	default=5, 
	type=int)
parser.add_argument("-n", "--name",
	help="name of files (default: same as 'input')")
parser.add_argument("-f", "--folder", 
	help="name of the parent folder (default: same as 'input')")
parser.add_argument("-cs", "--clustersize", 
	help="entries/frament files (default: %(default)s)", 
	default=25, 
	type=int)
parser.add_argument("-ft", "--filetype", 
	help="type of the generated files (default: %(default)s)", 
	default=".csv")
parser.add_argument("-sep", "--separator", 
	help="separator between columns in the output files (default: %(default)s)", 
	default=";")



args = parser.parse_args()

# Check if the given action is valid and replace with the list if == 'a'
while all(ch in "aigf" for ch in args.action) == False:
    args.action = input("\nWhich action do you want to conduct?"
        "\n- a\tconduct all actions\n- i\tdownload sequence IDs"
        "\n- g\tget neighbors\n- f\tfilter with target\n"
        "\nPlease enter any or multiple of letters (e.g 'a' or 'igf' [without ''])\n")
if "a" in args.action:
    args.action = "igf"

# Check if targets were given
if all(target == None for target in [args.targetID, args.targetDomain,
	args.targetName, args.targetFile]):
	while True:
		if "f" not in args.action:
			break
		Continue = input("\nNo target for filtering were given\t->Do you want to continue without?"
			"\n(y=yes, n=no)\n")
		if Continue == "y":
			break
		elif Continue == "n":
			print("Add the argmuments -ti (--targetID), -td (--targetDomain), -tn (--targetName) or -tf (--targetFile)"
				"\nEnter 'Main.py --help' for more info")
			quit()
		else:
			print("Please enter 'y' or 'no'!")

# Set file and folder name
if args.name == None:
	if os.path.isfile(args.input):
		args.folder, args.name = os.path.split(args.input)
		args.name = args.name.rsplit(".",1)[0]
	else:
		args.name = input("\nPlease enter a name for the created files (e.g. Test)\n")
		if args.folder == None:
			args.folder = args.name



# Check for tab separator
if args.separator in ["\\t", "tab", "'\\t'", "{tab}"]:
	args.separator = "\t"

## ------------------------------------------------------------------------------------------------
## MAIN FUNCTIONS ---------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------
## ================================================================================================
## Get index list of neighbors and retrieves protein data
def GetNeighbors(IDList, FilePath, Range, FileType, Sep, Ask, ClusterSize):
	Organisms = None
	print("Download protein data for", len(IDList), "IDs . . .")

	# Create clusters of sequences to generate smaller files (in case the download crashes)
	ClusteredList = [IDList[x:x+ClusterSize] for x in range(0, len(IDList), ClusterSize)]
	for ClusterID in range(len(ClusteredList)):
		print("Download cluster", ClusterID+1, "of", len(ClusteredList))
		FragmentFile = FilePath + "_" + str(ClusterID+1)
		print(FragmentFile)

		# Ignore all files that have already been downloaded
		if os.path.exists(FragmentFile + FileType):
			print("File already exists, skip to next cluster\n")

		# Download all files that have not yet been saved
		else:
			Neighbors = []
			Count = 0
			for GeneID in ClusteredList[ClusterID]:
				try:
					# Cycle through step size until the correct one is found
					ProteinSet = KEGG.DownloadNeighbors(GeneID, Range, Step=1)
					if len(ProteinSet) < Range + 1:
						ProteinSet =  KEGG.DownloadNeighbors(GeneID, Range, Step=5)
					if len(ProteinSet) < Range + 1:
						ProteinSet =  KEGG.DownloadNeighbors(GeneID, Range, Step=10)

					# Check if all entries for the range have been found 
					if len(ProteinSet) == Range*2:
						Status = "Complete"
						Count += 1
					else:
						Status = "Incomplete"
				except:
					Status = "Error"
					ProteinSet = [{"ID":GeneID}]
				for Protein in ProteinSet:
					Protein["Status"] = Status
					Neighbors.append(Protein)

			# Only download the list of organisms on KEGG if needed and add to dataframe
			if Organisms is None:
				Organisms = KEGG.DownloadOrganismsTemp()
			ProteinTable = pd.DataFrame(Neighbors)
			ProteinTable = pd.merge(ProteinTable, Organisms, on=["orgID"],  how="left")
			IE.ExportDataFrame(ProteinTable, FragmentFile, FileType=FileType, Sep=Sep, Ask=Ask)
			print("Done!\n->Neighbors found:", Count, "of", len(IDList), "complete")

	# After all entries have been downloaded, combine all fragments into one dataframe
	DataFrame = IE.CombineFiles(os.path.split(FragmentFile)[0], Sep, FileType)
	DataFrame["Length"] = DataFrame["Length"].fillna(0).astype(int)
	return(DataFrame)


def GetTargets(targetID, targetDomain, targetName, targetFile, Sep):
	TargetDict = {}
	TargetDict = {"KO-ID": targetID, "Name": targetName, "Domain": targetDomain}
	for Target in TargetDict:
		if TargetDict[Target]:
			TargetDict[Target] = TargetDict[Target].split(",")
		else:
			TargetDict[Target] = []
	if targetFile != None:
		with open(targetFile) as File:
			for Line in File:
				(Target, TargetType) = Line.strip().split(Sep)
				try:
					TargetDict[TargetType].append(Target)
				except KeyError:
					if TargetType != "Type":
						print(f"\nDetected {TargetType} not in type list, will be ignored")
	print(f"\nThe input targets are:\n{TargetDict}\n")
	return(TargetDict)


## ------------------------------------------------------------------------------------------------
## SCRIPT -----------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------

# Print program header
print('\n{:=<70}'.format(''))
print('{:=^70}'.format('  VICINITY ANALYZER  '))
print('{:=<70}'.format(''))
print('{: ^70}\n\n'.format('2024, by A.L.O. Gaenssle'))

IE.CreateFolder(os.path.join(args.folder, "VicinityAnalysis"))
OutputName =  os.path.join(args.folder, "VicinityAnalysis", args.name)


# Retrieve Sequence IDs from input, file or KEGG
if any(s in ["i", "g"] for s in args.action):
	OutputPath = os.path.join(args.folder, args.name)
	if re.search(r"^K\d+$", args.input):
		print("The input is a KO ID")
		IDList, DataFrame = KEGG.DownloadOrthology(args.input)
	elif os.path.exists(args.input):
		print("The input is a file")
		DataFrame = pd.read_csv(args.input, sep=args.separator)
		IDList = DataFrame[DataFrame.columns[0]].to_list()
	elif ":" in args.input:
		print("The input is a (list of) gene ID(s)")
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


# Get data from neighboring genes on KEGG
if "g" in args.action:
	OutputPath = OutputName + "_Neighbors"
	FragmentFolder = IE.CreateFolder(OutputPath + "_Fragments")
	FragmentFile = os.path.join(FragmentFolder, args.name + "_Neighbors")
	Detailed = GetNeighbors(IDList, FragmentFile, args.range, 
		args.filetype, args.separator, args.askoverwrite, args.clustersize)
	IE.ExportDataFrame(Detailed, OutputPath, 
		FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)

# Get 
if "f" in args.action:
	OutputPath = OutputName + "_Neighbors"
	ProteinData = pd.read_csv(OutputPath + args.filetype, sep=args.separator)

	# Set up dictionary of targets with Input:Type
	TargetDict = GetTargets(args.targetID, args.targetDomain, args.targetName, 
		args.targetFile, args.separator)

	# Cycle through all target and add boolean column
	TargetColumns = []
	for TargetType in TargetDict:
		for Target in TargetDict[TargetType]:
			NewColumn = TargetType[:2] + "-" + Target
			SearchColumn = TargetType
			ProteinData[NewColumn] = ProteinData[SearchColumn].str.contains(Target)
			TargetColumns.append(NewColumn)

	# Count occcurences of each target at each range position
	RangeCount = ProteinData.groupby('Pos')[TargetColumns] \
		.apply(sum).reset_index()
	print(f"\nFOUND TARGETS\n{RangeCount}\n\n")
	IE.ExportDataFrame(RangeCount, OutputPath + "_RangeCount", 
		FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)

	# Count occcurences of each target for each entry
	EntryCount = ProteinData.copy().groupby('Ref')[TargetColumns] \
		.apply(sum).reset_index()
	IE.ExportDataFrame(EntryCount, OutputPath + "_EntryCount", 
		FileType=args.filetype, Sep=args.separator, Ask=args.askoverwrite)

print('{:=^70}'.format('  End of program  '))