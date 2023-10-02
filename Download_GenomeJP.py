#!/usr/bin/python
# Written in Python 3.8 in 2023 by A.L.O. Gaenssle

# MODULE: DOWNLOAD GENE/PROTEIN DATA from Genome.jp
# -> downloads all IDs associated with the input (domain) id
# -> exports the otained data into a table
# -> download protein data from UniProt (e.g. organims, domains and sequence)

import pandas as pd
import re
import math
import urllib.request
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

##-------------------------------------------------------------------------------------------------
## DOWNLOAD FUNCTIONS -----------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------
## ================================================================================================
## Download gene list from Genome.jp
def DownloadGeneList(url, getAmount=False):
	Amount = 0
	with urllib.request.urlopen(url) as File:
		Add = False
		List = []
		for Line in File:
			Line = Line.decode("utf-8").strip()
			Line = re.sub('<[^>]*>', '', Line)
			if Line:
				if Add == False:
					if Line.startswith("Hits:"):
						Amount = Line.split(" ",2)[1]
						Amount = math.ceil(int(Amount)/1000)

					# Start of list of entries
					elif Line.startswith("----------"):
						Add = True
				else:
					# End of list of entries
					if Line.startswith("DBGET integrated"):
						break
					else:
						List.append(Line.strip())
	print(url, "downloaded")
	print(len(List), "gene IDs added")
	if getAmount:
		return(List, Amount) 	# Required for the first call to get the number of pages
	else:
		return(List)	# For all subsequent pages

## ================================================================================================
## Download UniProt protein entry
def DownloadEntryUniProt(ID):
	url = "https://www.genome.jp/entry/up:" + str(ID)
	with urllib.request.urlopen(url) as Input:
		Dict = {"ID": ID, "Sequence": ""}
		TaxString = ""
		InSequence = False
		InDomain = False
		DomainID = 0
		for Line in Input:
			Line = Line.decode("utf-8").strip()
			Line = re.sub('<[^>]*>', '', Line)
			if Line:
				try:
					if Line.startswith("ID"):
						Dict["Length"] = Line.rsplit(";",1)[1].rsplit(" ",1)[0].strip()
					elif Line.startswith("OS"):
						Dict["Organism"] = Line.split(" ",1)[1].replace(".", "").strip()
					elif Line.startswith("OC"):
						TaxString += Line.split(" ",1)[1].replace(".", "").strip()
					elif Line.startswith("DR"):
						if "KEGG" in Line:
							Dict["KEGG"] = Line.split(";")[1].strip()

					# Each domain is stored in 3 columns: Name-D[n], Start-D[n], End-D[n]
					elif Line.startswith("FT"):
						if "DOMAIN" in Line:
							InDomain = True
							DomainID += 1
							DID = "-" + "D" + str(DomainID)
							String = Line.rsplit(" ",1)[1]
							Dict["Start" + DID], Dict["End" + DID] = String.split("..",1)
						elif InDomain and "/note=" in Line:
							Dict["Name" + DID] = Line.split("\"")[1]
							InDomain = False

					# The amino acid sequence is the last information on the page
					elif InSequence:
						# Separator between data sets of different proteins
						if Line.startswith("//"):
							break
						else:
							Dict["Sequence"] += Line.replace(" ", "")
					elif Line.startswith("SQ"):
						InSequence = True
				except IndexError:
					pass

	# Extract and store the taxonomic classification and Kingdom-Phylum
	try:
		Kingdom, Phylum = TaxString.split("; ",3)[:2]
		Dict["Taxonomy"] = Kingdom + "-" + Phylum
	except ValueError:
		Dict["Taxonomy"] = "-".join(TaxString.split("; ",))
	print("ID", ID, "downloaded")
	return(Dict)


##-------------------------------------------------------------------------------------------------
## CLEANUP DATA FUNCTIONS -------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------
## ================================================================================================
## Convert downloaded KEGG gene text to table
def CleanKEGG(Data):
	ListOfDicts = []
	for Line in Data:
		Dict = {}
		try:
			Dict["ID"], String = Line.strip().split(" ",1)
			if "no KO assigned" in String:
				Dict["Description"] = String.split("|")[1].split(") ")[1]
			else:
				if "|" in String:
					String = String.split("|")[0]
				if "[EC" in String:
					String, Dict["#EC"] = String.rsplit(" [EC:")
					Dict["#EC"] = Dict["#EC"].replace("]", "").strip()
				Dict["KO ID"], Dict["Description"] = String.strip().split(" ", 1)
		except ValueError:
			pass
		ListOfDicts.append(Dict)
		DataFrame = pd.DataFrame(ListOfDicts)
		DataFrame = DataFrame[["ID", "KO ID", "#EC", "Description"]]
	return(DataFrame)

## ================================================================================================
## Convert downloaded UniProt/SWISS-Prot gene text to table
def CleanUniProt(Data):
	ListOfDicts = []
	for Line in Data:
		Dict = {}
		try:
			Dict["ID"], String = Line.split(" ",1)
			String = String.strip().split("Full=",1)[1]
			Dict["Description"], IDString = String.split(" {",1)
			IDList = IDString.split(",",1)[0].split("}",1)[0].split("|")
			for ID in IDList:
				if len(ID) > 10:
					Type, ID = ID.strip().split(":",1)
					Dict[Type] = ID
			ListOfDicts.append(Dict)
		except ValueError:
			pass
	DataFrame = pd.DataFrame(ListOfDicts)
	return(DataFrame)

## ================================================================================================
## Convert downloaded linked gene text to table (e.g. PDB)
def CleanPDB(Data):
	ListOfDicts = []
	for Line in Data:
		Dict = {}
		try:
			Dict["ID"], Dict["Description"] = Line.strip().split(" ",1)
		except ValueError:
			pass
		ListOfDicts.append(Dict)
		DataFrame = pd.DataFrame(ListOfDicts)
	return(DataFrame)
