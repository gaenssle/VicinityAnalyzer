#!/usr/bin/python
# Written in Python 3.8 in 2023 by A.L.O. Gaenssle

# MODULE: DOWNLOAD PROTEIN DATA from KEGG
# -> downloads information of protein entries by ID in chunks
# -> downloads all organism ids available on KEGG and their taxonomic classification
# -> downloads all neighbors within the given range of each given protein ID

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
	GeneList = Download.strip().split("\n")
	ListOfList = [i.split("\t") for i in GeneList]
	DataFrame = pd.DataFrame(ListOfList, columns=["ID", "Description"])
	return(DataFrame["ID"].to_list(), DataFrame)


## ================================================================================================
## Download all genome taxonomy from KEGG --> KEGG-list
def DownloadOrganismsTemp(Name="organism"):
	print("Download organism taxonomy. . .")
	Entry = REST.kegg_list(Name).read()
	Entry = Entry.replace(";" , "\t")
	ColList = ["ID long", "orgID", "Organism", "Kingdom", "Phylum", "Class", "Order"]
	DataFrame = pd.read_csv(StringIO(Entry), sep="\t", names=ColList)
	DataFrame["Taxonomy"] = DataFrame["Kingdom"] + "-" + DataFrame["Phylum"]
	DataFrame = DataFrame[["orgID","Taxonomy"]]
	return(DataFrame)


##-------------------------------------------------------------------------------------------------
## SUB-FUNCTIONS OF DownloadNeighbors--------------------------------------------------------------
##-------------------------------------------------------------------------------------------------
## ================================================================================================
## Get list of indexes +/- range of the reference gene ID for KEGG
def GetNeighborIndices(Gene, Range, Step, Size=4):
	IndexList = []
	IndexDict = {}
	if "_" in Gene and Gene.rsplit("_",1)[1].isdigit():
		Label = Gene.rsplit("_",1)[0] + "_"
		Index = int(Gene.rsplit("_",1)[1])
		Fill = len(Gene.rsplit("_",1)[1])
	else:
		Label = Gene[:-Size]
		Index = int(Gene[-Size:])
		Fill = Size
	for i in range(Index-Range*Step,Index+Range*Step+1, Step):
		if i != Index:
			NewID = Label + str(i).zfill(Fill)
			List.append(NewID)
			Dict[NewID] = int((i - Index)/Step)
	return(List, Dict)


## ================================================================================================
## Download protein entries from KEGG -> in chunks of 10 gene IDs --> KEGG-get
def DownloadProteinEntries(List, GeneID):
	Data = []
	Entry = []
	print("Download neighbors of", GeneID, ". . .")
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
				Dict["Name"] = Line.split(" ",1)[1].replace("(GenBank)", "").strip()
			elif Line.startswith("ORTHOLOGY"):
				Dict["KO-ID"] = Line.split(" ",2)[1].strip()
			elif Line.startswith("ORGANISM") or Line.startswith("VIRUS"):
				Line = Line.split(" ",1)[1].strip()
				Dict["Organism"] = Line.split(" ",1)[1]
			elif Line.startswith("MOTIF"):
				Dict["Domain"] = Line.split(" ",1)[1].replace("Pfam:", "").strip()
			elif "UniProt" in Line:
				Dict["UniProt"] = Line.split(" ",1)[1]
			elif Line.startswith("AASEQ"):
				Dict["Length"] = Line.split(" ",1)[1]
				inAASeq = True
	return(Dict)


##-------------------------------------------------------------------------------------------------
## MAIN FUNCTION ----------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------
## ================================================================================================
## Main function to download neighbors
def DownloadNeighbors(GeneID, Range, Step=1):
	print(f"Download neighbors of {GeneID} (Increment={Step}) . . .")
	Data = []
	ProteinSet = []
	IDList, RangeDict = GetNeighborIndices(GeneID, Range, Step)
	if Range > 5:
		# Create chunks of clusters since data of 10 proteins can be downloaded from KEGG at once
		ClusteredList = [IDList[x:x+10] for x in range(0, len(IDList), 10)]
	else:
		ClusteredList = [IDList]
	for Cluster in ClusteredList:
		Data.extend(DownloadProteinEntries(Cluster, GeneID))
	for Entry in Data:
		ProteinSet.append(GetDetailedData(Entry, GeneID, GeneID.split(":",1)[0]))
	for Entry in ProteinSet:
		Entry["Pos"] = RangeDict[Entry["ID"]]
	return(ProteinSet)