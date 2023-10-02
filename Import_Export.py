#!/usr/bin/python
# Written in Python 3.8 in 2023 by A.L.O. Gaenssle

# MODULE: EXPORT DATA AND CHECK FILES
# -> create new folder
# -> check if files exists and give an option to rename the file
# -> Combine files (sets of entries) into a large dataframe
# -> export pandas dataframe to defined filetype with set separator (Main.py)

import pandas as pd
import os
import re

##-------------------------------------------------------------------------------------------------
## HELPER FUNCTIONS -------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------
## ================================================================================================
## Create new Folder
def CreateFolder(NewPath):
	if not os.path.exists(NewPath):
		os.makedirs(NewPath)
		print("Created folder:", NewPath)
	else:
		print("Files will be added to:", NewPath)
	return(NewPath)

## ================================================================================================
## Check if the file already exists and if it should be replaced
def CheckFileExists(FileName, Ask):
	if Ask:
		Replace = "n"
	else:
		Replace = "y"
	while Replace == "n":
		if not os.path.exists(FileName):
			break
		else:
			Replace = input("\nFile " + FileName + " already exits -> should it be replaced?"
				"\n(y=yes, n=no)\n")
			while Replace not in ("y", "n"):
				Replace = input("\nPlease enter 'y' or 'n'!\n")
		if Replace == "n":
			FileName = input("\nEnter a new filename\n")
	return(FileName)

## ================================================================================================
## Copy all data from Fragment files (250-500 genes/file) into one large file
def CombineFiles(Folder, Sep, FileType):
	FileListAll = os.listdir(Folder)
	FileList = list(filter(lambda File: File.endswith(FileType), FileListAll))
	DataList = []
	for File in FileList:
		FilePath = os.path.join(Folder, File)
		DataFrame = pd.read_csv(FilePath, sep=Sep)
		DataList.append(DataFrame)
	DataFrame = pd.concat(DataList, axis=0, ignore_index=True)
	DataFrame = DataFrame.sort_values(DataFrame.columns[0])
	return(DataFrame)


##-------------------------------------------------------------------------------------------------
## EXPORT FILE FUNCTION ---------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------
## ================================================================================================
## Export pandas dataframe
def ExportDataFrame(DataFrame, FileName, Add="", Columns="", 
	FileType=".csv", Sep=";", Ask=True, Header=True):
	FileName = FileName + Add + FileType
	FileName = CheckFileExists(FileName, Ask)
	if Columns == "":
		Columns = list(DataFrame)
	DataFrame.to_csv(FileName, sep=Sep, columns = Columns, index=False, header=Header)
	print("File saved as:", FileName, "\n")


## ================================================================================================
## Export Fasta of either the full sequence or only the domain
def CreateFasta(df, FilePath, only=False):
    df_fasta = df.copy()
    df_fasta = df_fasta.reset_index()   
    if only:
        FilePath = FilePath + "_only"
        df_fasta['ID'] = df_fasta['ID'] + '_' + df_fasta['Domain']
    FilePath = FilePath  + ".fasta"

    df_fasta["fasta"] = df_fasta.agg(lambda x: 
    	f">{x['ID']} [{x['Organism']}] {x['Taxonomy']}\n{x['Sequence']}\n", axis=1)
    with open(FilePath, "w") as f_out:
        f_out.write("\n".join(df_fasta["fasta"]))