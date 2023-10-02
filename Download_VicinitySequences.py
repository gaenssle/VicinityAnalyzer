#!/usr/bin/python
# Written in Python 3.7 in 2023 by A.L.O. Gaenssle

import Import_Export as IE	# Own module



##------------------------------------------------------
## FUNCTIONS
##------------------------------------------------------

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


def SplitFile(InputFile, TagList):
	Export = {}
	Data, Header = IE.ImportNestedList(InputFile, getHeader=True)
	IndexList = GetColumnList(Header, TagList)
	for Gene in Data:
		try:
			Type = Gene[IndexList[1]]
			if Type in Export:
				Export[Type].append(Gene[IndexList[0]])
			else:
				Export[Type] = [Gene[IndexList[0]]]
		except IndexError:
			pass
	for Type in Export:
		IE.ExportList(Export[Type], InputFile, Add="_" + Type, Ask=False)


def GetNeighbors(InputFile, Identifier):
	Data, Header = IE.ImportDoubleNestedDictionary(InputFile, getHeader=True)
	ColumnList = GetColumnList(Header, TagList, Offset=-1)
	Extract = {}
	Extract_List = []
	Extract_Failed = []
	Count = 0
	for Gene in Data:
		for Neighbor in Data[Gene][2:8]:
			for Tag in range(len(ColumnList)):
				try:
					if Identifier[Tag] in Neighbor[ColumnList[Tag]]:

						List = [ConvertToOffset(Data[Gene].index(Neighbor)), TagList[Tag], Neighbor[0], Neighbor[ColumnList[Tag]]]
						if Gene in Extract:
							Extract[Gene].append(List)
						else:
							Extract[Gene] = [List]
						Extract_List.append(Neighbor[0])
						break
				except IndexError:
					pass
		if Gene not in Extract:
			Extract[Gene] = [[""] * len(List)]
			Count += 1
			Extract_Failed.append(Gene)
		# else:
		# 	if len(Extract[Gene]) > 1:
		# 		print(Extract[Gene])
	# print(Extract)
	print(Count, "of", len(Data), "not found")
	OutputHeader = "GeneID\tPos\tHit\tNeighbor\tHit_Text\n"
	IE.ExportDoubleNestedDictionary(Extract, InputFile, OutputHeader, Add="_Extract", Ask=False)
	IE.ExportList(Extract_List, InputFile, Add="_Extract_IDs", Ask=False)
	IE.ExportList(Extract_Failed, InputFile, Add="_Extract_Failed", Ask=False)




##------------------------------------------------------
## MAIN SCRIPT
##------------------------------------------------------

# TagList = ["KEGGID", "Name", "Pfam"]
# Identifier = ["K21572", "SusD", "SusD"] # ID, Name, Domain

Folder = "Vicinity/"
Domain = "DUF1735"


# InputFile = Folder + Domain + "_KEGG_Cutoff-0.0001_List_Neighbors.txt"
# InputFile = Folder + Domain + "_KEGG_List_Test2_Neighbors.txt"
# GetNeighbors(InputFile, Identifier)

TagList = ["Neighbor", "Hit"]
InputFile = Folder + Domain + "_KEGG_Cutoff-0.0001_List_Neighbors_Extract.txt"
SplitFile(InputFile, TagList)