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
## SUBFUNCTION: of DownloadProteinEntries()
## Download Info for each protein from KEGG
def GetDetailedData(Entry, ID):
	Dict = {"ID": ID, "orgID": ID.split(":",1)[0],"Sequence": ""}
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
			if Line.startswith("ORGANISM") or Line.startswith("VIRUS"):
				Line = Line.split(" ",1)[1].strip()
				Dict["Organism"] = Line.split(" ",1)[1]
			elif "UniProt" in Line:
				Dict["UniProt"] = Line.split(" ",1)[1]
			elif Line.startswith("AASEQ"):
				Dict["Length"] = Line.split(" ",1)[1]
				inAASeq = True
	return(Dict)

## ================================================================================================
## Download protein entries from KEGG -> in chunks of 10 gene IDs --> KEGG-get
def DownloadProteinEntries(Chunked_List):
	Data = []
	Entry = []
	print("Download protein info, set of", Chunked_List[0], ". . .")
	Download = REST.kegg_get(Chunked_List).read()
	Download = Download.split("\n")
	Count = 0
	for Line in Download:
		if Line.startswith("///"):
			ProteinData = GetDetailedData(Entry, Chunked_List[Count])
			Data.append(ProteinData)
			Entry = []
			Count += 1
		else:
			Entry.append(Line)
	return(Data)

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
	print(DataFrame.head())
	return(DataFrame)

## ================================================================================================
## Download domain motifs (architecture) for each gene ID from KEGG
def DownloadMotif(ID):
	Data = []
	Domains = []
	Row = ""
	inDomains = False
	Index = 1
	url = "https://www.kegg.jp/ssdb-bin/ssdb_motif?kid=" + ID
	with urllib.request.urlopen(url) as File:
		for Line in File:
			Line = Line.decode("utf-8").strip().replace(" : ", "")
			Line = re.sub('<[^>]*>', '|', Line)
			Line = Line.replace("&nbsp;", "-")
			Line = Line.split("|")
			Line = [i for i in Line if i]
			if Line:
				Data.append(Line)
		for Line in Data:
			if inDomains:
				if Line[0].startswith("pf:"):
					if Row:
						Index += 1
						Domains.append(Row)
					Row = [ID, Index, Line[0].split(":",1)[1]] + Line[1:]
				else:
					if Line[0].startswith("["):
						Domains.append(Row)
						break
					Row += Line
			elif Line[0] == "Motif id":
				inDomains = True
	print(ID, "downloaded")
	return(Domains)
