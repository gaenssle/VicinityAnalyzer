#!/usr/bin/python
# Written in Python 3.7 in 2022 by A.L.O. Gaenssle


import re
import Bio
from Bio.KEGG import REST
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

##------------------------------------------------------

import Import_Export as IE
import Download_KEGG as KEGG

##------------------------------------------------------
## FUNCTIONS
##------------------------------------------------------

# Get list of indexes +/- range of the reference gene ID for KEGG
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


# Download protein entries from KEGG -> in chunks of 10 gene IDs --> KEGG-get
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


# Download Info for each protein from KEGG
def GetDetailedData(Entry):
	NCBI = ""
	UniProt = ""
	Pfam = ""
	KEGGID = ""
	ProteinName = ""
	for Line in Entry:
		Line = re.sub("\s\s+" , " ", Line)
		Line = Line.strip()
		if Line.startswith("AASEQ"):
			break
		elif Line.startswith("ENTRY"):
			Entry = Line.split(" ",2)[1]
		elif Line.startswith("NAME"):
			Line = Line.split(" ",1)[1]
			Name = Line.replace("(GenBank) ", "")
		elif Line.startswith("ORTHOLOGY"):
			Line = Line.split(" ",1)[1].strip()
			KEGGID, ProteinName = Line.split(" ",1)
		elif Line.startswith("ORGANISM"):
			Entry = Line.split(" ",2)[1] + ":" + Entry
		elif Line.startswith("MOTIF"):
			Line = Line.split(" ",1)[1]
			Pfam = Line.replace("Pfam: ", "")
		elif "NCBI-ProteinID:" in Line:
			Line = Line.split(" ",1)[1]
			NCBI = Line.replace("NCBI-ProteinID: ", "")
		elif "UniProt" in Line:
			UniProt = Line.split(" ",1)[1]
	ProteinData = [Entry, Name, KEGGID, ProteinName, NCBI, UniProt, Pfam]
	return(ProteinData)

# Get the list of indices for the needed columns
def GetColumnList(Header, TagList, Offset=0):
	IndexList = []
	HeaderList = Header.strip().split("\t")
	for Tag in TagList:
	    try:
	        IndexList.append(HeaderList.index(Tag)+Offset)
	    except:
	        print("Error", Tag)
	return(IndexList)

def ConvertToOffset(Index, Range=5):
	Offset = Index - Range + 1
	if Index < Range:
		Offset -= 1
	return(Offset)

##------------------------------------------------------
## MAIN FUNCTION
##------------------------------------------------------

# Main function: Gets index list of neighbors and retrieves protein data
def GetNeighbors(GeneList, InputFile, Ask=False):
	ManualList = []
	Complete = {}
	Incomplete = {}
	for GeneID in GeneList:
		ProteinList = []
		if GeneID.split(":",1)[1].isdigit() == False:
			try:
				List, ManualList = GetNeighborIndices(GeneID, ManualList)
				Data = DownloadProteinEntries(List, GeneID)
				if len(Data) < 6:
					List, ManualList = GetNeighborIndices(GeneID, ManualList, Multi=5)
					Data = DownloadProteinEntries(List, GeneID)
			except:
				try:
					List, ManualList = GetNeighborIndices(GeneID, ManualList, Multi=5)
					Data = DownloadProteinEntries(List, GeneID)
					if len(Data) < 6:
						List, ManualList = GetNeighborIndices(GeneID, ManualList, Multi=10)
						Data = DownloadProteinEntries(List, GeneID)
				except:
					ManualList.append(GeneID)
			if len(Data) == 10:
				for Entry in Data:
					ProteinList.append(GetDetailedData(Entry))
				Complete[GeneID] = ProteinList
			else:
				for Entry in Data:
					ProteinList.append(GetDetailedData(Entry))
				Incomplete[GeneID] = ProteinList
		else:
			ManualList.append(GeneID)
	print("Neighbors found:", len(Complete), "complete,", len(Incomplete), "incomplete; of", len(GeneList))
	print("Manual search required for:", len(ManualList))
	Header = "Domain\tEntry\tName\tKEGGID\tProteinName\tNCBI\tUniProt\tPfam\n"
	IE.ExportDoubleNestedDictionary(Complete, InputFile, Header, Add="_Neighbors", Ask=Ask)
	if Incomplete != {}:
		IE.ExportDoubleNestedDictionary(Incomplete, InputFile, Header, Add="_Neighbors_incomplete", Ask=Ask)
	if ManualList != []:
		IE.ExportList(ManualList, InputFile, Add="_ManualDownload", Ask=Ask)


# Count neighbors based on KEGG-ID, Name and Domain
def FindNeighbors(Data, ColumnList, IDList, NameList, DomainList, InputFile, Range=5):
	VicinityAll = {}
	for GeneID in Data:
		VicinityList = []
		for Neighbor in Data[GeneID]:
			Label = ""
			try:
				if Neighbor[ColumnList[0]] in IDList:
					Label = IDList[Neighbor[ColumnList[0]]] + "-ID"
					# print(GeneID, ConvertToOffset(Data[GeneID].index(Neighbor)), Label)
				else:
					for Name in NameList:
						if Name in Neighbor[ColumnList[1]]:
							Label = NameList[Name] + "-N"
							break
					if Label == "":
						for Domain in DomainList:
							if Domain in Neighbor[ColumnList[2]]:
								Label = DomainList[Domain] + "-D"
								break
			except IndexError:
				# print(Neighbor)
				pass
			VicinityList.append(Label)
		if len(VicinityList) != 10:
			print(GeneID, len(VicinityList), VicinityList)
		VicinityAll[GeneID] = VicinityList
		# print(VicinityList)
	RangeLow = "\t".join(str(x) for x in list(range(-Range,0)))
	RangeHigh = "\t".join(str(x) for x in list(range(1,Range+1)))
	HeaderOutput = "GeneID\t" + RangeLow + "\t" + RangeHigh + "\n"
	IE.ExportNestedDictionary(VicinityAll, InputFile, HeaderOutput, Add="_HitList", Ask=False)
	return(VicinityAll)

def CountNeighbors(Vicinity):
	HitList = {}
	HitCount = 0
	for GeneID in Vicinity:
		Hit = False
		for Offset in Vicinity[GeneID]:
			if Offset != "":
				Hit = True
		if Hit:
			HitCount += 1
	print(HitCount, "of", len(Vicinity), "found")


##------------------------------------------------------
## SCRIPT
##------------------------------------------------------

Folder = "Vicinity/"
Domain = "BACON_2"
# InputFile = Folder + Domain + "_KEGG_Cutoff-0.0001_List.txt"
# InputFile = Folder + Domain + "_KEGG_List_Test2.txt"

# GeneList = IE.ImportList(InputFile)
# GetNeighbors(GeneList, InputFile)

InputFile = Folder + Domain + "_KEGG_Cutoff-0.0001_List_Neighbors.txt"
# InputFile = Folder + Domain + "_KEGG_List_Test2_Neighbors.txt"

TagList = ["KEGGID", "Name", "Pfam"]
ReferenceFile = Folder + "Transporter_Names"
IDList = IE.ImportDictionary(ReferenceFile + "_IDs.txt")
NameList = IE.ImportDictionary(ReferenceFile + ".txt")
DomainList = IE.ImportDictionary(ReferenceFile + "_Domains.txt")

Data, Header = IE.ImportDoubleNestedDictionary(InputFile, getHeader=True)
ColumnList = GetColumnList(Header, TagList, Offset=-1)
# print(ColumnList)
VicinityAll = FindNeighbors(Data, ColumnList, IDList, NameList, DomainList, InputFile)
CountNeighbors(VicinityAll)
