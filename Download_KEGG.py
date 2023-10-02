#!/usr/bin/python
# Written in Python 3.8 in 2023 by A.L.O. Gaenssle

# MODULE: DOWNLOAD PROTEIN DATA from KEGG
# -> downloads information of protein entries by ID in chunks
# -> downloads all organism ids available on KEGG and their taxonomic classification
# -> downloads all motifs (domain architecture) of each given protein ID

import pandas as pd
from io import StringIO
import re
import Bio
from Bio.KEGG import REST
import urllib.request
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

##-------------------------------------------------------------------------------------------------
## DOWNLOAD FUNCTIONS -----------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------
## ================================================================================================
## Download all gene IDs associated with the supplied KEGG Orthology (KO)
def DownloadOrthology(Input):
	Download = REST.kegg_find("genes",Input).read()
	List = Download.strip().split("\n")
	ListOfList = [i.split("\t") for i in List]
	DataFrame = pd.DataFrame(ListOfList, columns=["ID", "Description"])
	return(DataFrame["ID"].to_list(), DataFrame)


## ================================================================================================
## Get list of indexes +/- range of the reference gene ID for KEGG
def GetNeighborIndices(Gene, ManualList, Multi=1, Range=5, Size=4):
	try:
		List = []
		if "_" in Gene and Gene.rsplit("_",1)[1].isdigit():
			Label = Gene.rsplit("_",1)[0] + "_"
			Index = int(Gene.rsplit("_",1)[1])
			Fill = len(Gene.rsplit("_",1)[1])
		else:
			Label = Gene[:-Size]
			Index = int(Gene[-Size:])
			Fill = Size
		for i in range(Index-Range*Multi,Index+Range*Multi+1, Multi):
			if i != Index:
				NewID = Label + str(i).zfill(Fill)
				List.append(NewID)
	except:
		ManualList.append(Gene)
	return(List, ManualList)


## ================================================================================================
## Download protein entries from KEGG -> in chunks of 10 gene IDs --> KEGG-get
def DownloadProteinEntries(List, GeneID):
	Data = []
	Entry = []
	print("Download protein info for", GeneID, ". . .")
	Download = REST.kegg_get(List).read()
	Download = Download.split("\n")
	for Line in Download:
		if Line.startswith("///"):
			Data.append(Entry)
			Entry = []
		else:
			Entry.append(Line)
	return(Data)

## ================================================================================================
## Download Info for each protein from KEGG
def GetDetailedData(Entry, GeneID, orgID):
	Dict = {"Ref": GeneID,"ID": orgID, "orgID": orgID,"Sequence": ""}
	inAASeq = False
	for Line in Entry:
		Line = re.sub("\s\s+" , " ", Line)
		Line = Line.strip()
		if inAASeq == True:
			if Line.startswith("NTSEQ"):
				break
			else:
				Dict["Sequence"] += Line.strip()
		if inAASeq == False:
			if Line.startswith("ENTRY"):
				Dict["ID"] += ":" + Line.split(" ",2)[1]
			elif Line.startswith("NAME"):
				Dict["Description"] = Line.split(" ",1)[1].replace("(GenBank)", "").strip()
			elif Line.startswith("ORTHOLOGY"):
				Dict["KO-ID"] = Line.split(" ",2)[1].strip()
			elif Line.startswith("ORGANISM") or Line.startswith("VIRUS"):
				Line = Line.split(" ",1)[1].strip()
				Dict["Organism"] = Line.split(" ",1)[1]
			elif Line.startswith("MOTIF"):
				Dict["Domains"] = Line.split(" ",1)[1].replace("Pfam:", "").strip()
			elif "UniProt" in Line:
				Dict["UniProt"] = Line.split(" ",1)[1]
			elif Line.startswith("AASEQ"):
				Dict["Length"] = Line.split(" ",1)[1]
				inAASeq = True
	return(Dict)