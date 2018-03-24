#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict
import string
import shlex
import time
import subprocess
import glob
import shutil

#function that detect the python version
def python_version():
	if(sys.version_info >= (3,0,0)):
		return(3)
	else:
		return(2)

#parse gene location from a file formatted: GeneName \t ContigName \t GeneStart \t GeneEnd
#use boolGeneBased = False when mapping when contigs are the reference
#checked
def getReferenceDict(geneLocationFileName, boolGeneBased=True):

	dictReference2geneLocation = defaultdict(list)
	geneLocationFile = open(geneLocationFileName, "r")

	for strLine in geneLocationFile:
		strLine = strLine.rstrip('\n')
		arrLine = strLine.split('\t')

		strGeneName = arrLine[0]
		strContigName = arrLine[1]
		strGeneStart = int(arrLine[2])
		strGeneEnd = int(arrLine[3])

		geneInfoTuple = (strGeneName, strGeneStart, strGeneEnd)

		if (boolGeneBased):
			dictReference2geneLocation[strGeneName].append(geneInfoTuple)
		else:
			dictReference2geneLocation[strContigName].append(geneInfoTuple)
	return(dictReference2geneLocation)


#parses AS flag, score=0 if there is no flag
#checked
#AS:i:64
#Could add another function in case AS flag isnt found and let this return logical if AS found or not
def parseSAMAlignmentScore(strSAMinfo):

	listSamInfo = strSAMinfo.split("\t")
	alignmentScore = 0
	alignmentScoreAS = -1
	for addOn in listSamInfo:
		if (addOn.startswith("AS:i:")):
			alignmentScoreString = addOn.split(":")[2]
			alignmentScoreAS = int(alignmentScoreString)
	if (alignmentScoreAS != -1):
		alignmentScore = alignmentScoreAS
	else:
		alignmentScore = 0

	return(alignmentScore)


#calculate if and by how much the alignment overlaps with the gene/exon
#the overlap count here does not account for gaps within the alignment
#checked
def checkAlignGeneOverlap(alignStart, alignEnd, geneStart, geneEnd):
	overlap = False
	overlapBaseCount = 0
	if (alignStart >= geneStart and alignStart <= geneEnd and alignEnd >= geneStart and alignEnd <= geneEnd):
		overlap = True
		overlapBaseCount = alignEnd - alignStart + 1
	elif (alignEnd >= geneStart and alignEnd <= geneEnd and alignStart < geneStart):
		overlap = True
		overlapBaseCount = alignEnd - geneStart + 1
	elif(alignStart >= geneStart and alignStart <= geneEnd and alignEnd >= geneEnd):
		overlap = True
		overlapBaseCount = geneEnd - alignStart + 1
	elif(alignStart < geneStart and alignEnd > geneEnd):
		overlap = True
		overlapBaseCount = geneEnd - geneStart+ 1

	return(overlap, overlapBaseCount)



##
#adds fields "overlapCount" and "overlap" to dictSAMline
#checked

#parseCigar
#checkAlignGeneOverlap
def calculateOverlap(dictSAMline, dictReference2geneLocation):

	(listCigarTypes, listNumbers, seqLength, alignLength, alignCigarMismatches, alignLengthOnRef) = parseCigar(dictSAMline["cigar"])
	alignStart = int(dictSAMline["alignstart"])
	alignEnd = alignStart + int(alignLengthOnRef) - 1

	dictSAMline["overlapCount"] = []
	dictSAMline["overlap"] = False

	listRefGenes = dictReference2geneLocation[dictSAMline["refname"]]
	for geneTuple in listRefGenes:
		(geneName, geneStart, geneEnd) = geneTuple

		overlap, overlapCount = checkAlignGeneOverlap(alignStart, alignEnd, geneStart, geneEnd)

		if (overlap and overlapCount > 0):
			dictSAMline["overlap"] = True
			dictSAMline["overlapCount"].append((geneName,overlapCount))

		elif(overlap and overlapCount == 0):
			sys.stderr.write("Overlap of 0?\n")
			dictSAMline["overlap"] = False
		else:
			dictSAMline["overlap"] = False



def filterInsert_getHighestScoringHit(listInsertSAMdicts, dictRef2alignmentScore):
	#find best alignmentScore
	listInsertSAMdicts_besthit = []
	bestAlignmentScore = 0
	if (len(listInsertSAMdicts) >= 1):
		for dictSAMline in listInsertSAMdicts:
			currRefname = dictSAMline["refname"]
			alignmentScore = dictRef2alignmentScore[currRefname]
			if ( alignmentScore > bestAlignmentScore ):
				listInsertSAMdicts_besthit = []
			if ( alignmentScore >= bestAlignmentScore ):
				bestAlignmentScore = alignmentScore
				listInsertSAMdicts_besthit.append(dictSAMline)

	return(listInsertSAMdicts_besthit)


def filterInsert_getPreferedAlignments(listInsertSAMdicts, setOtherReadsRefs):
	#prefered alignments are those that for which both pairs map on the same reference, only genes with overlap to at least one gene are used
	#setOtherReadsRefs: ref for the other direction of the insert (forward/reverse)
	normalSAMdicts = []
	preferedSAMdicts = []
	preferedRefs = set()

	for dictSAMline in listInsertSAMdicts:
		if (dictSAMline["overlap"]):
			currRefname = dictSAMline["refname"]
			if (currRefname in setOtherReadsRefs):
				preferedSAMdicts.append(dictSAMline)
				preferedRefs.add(currRefname)
			else:
				normalSAMdicts.append(dictSAMline)

	return(preferedSAMdicts, preferedRefs, normalSAMdicts)


#parse one line of SAM alignment to a dictionary
#checked
#469616.genomes.fasta:1:1:1:91#0/1       16      metaMG0013354.COG0541   98      255     75M     *       0       0       TTTTTAGAATGGAAATTTTCCACCCTTTCCACCTTTAAATCCACCAAATGATGGTAATGAAGGTATTTTTCCACC     hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh        AS:i:1000       NM:i:0
def parseSamLine(strSamline):

	try:

		dictSAMline = defaultdict(str)

		strSamline = strSamline.strip("\n")
		listSamline = strSamline.split("\t")

		strReadName = listSamline[0]
		strSAMFlag = listSamline[1]
		strRefName = listSamline[2]

		strAlignStart = listSamline[3]
		strMapQ = listSamline[4]
		strCIGAR = listSamline[5]
		strSAMinfo = '\t'.join(listSamline[11:])

		dictSAMline["readname"] = strReadName
		dictSAMline["samflag"] = strSAMFlag
		dictSAMline["refname"] = strRefName
		dictSAMline["alignstart"] = strAlignStart
		dictSAMline["mapq"] = strMapQ
		dictSAMline["cigar"] = strCIGAR
		dictSAMline["saminfo"] = strSAMinfo
		dictSAMline["samline"] = strSamline
		dictSAMline["alignmentScore"] = 0

	except:
		dictSAMline = None

	return(dictSAMline)


#parse SAM flag to a dictionary
#checked
def parseSAMflag(strSAMflag):

	dictSAMflag = defaultdict(bool)
	intSAMflag = int(strSAMflag)

	flagCodes = [(0x1, "paired"), (0x2, "properPaired"),(0x4, "unmapped"),(0x8, "mateUnmapped"),(0x10, "reverseStrand"),(0x20, "mateReverseStrand"),(0x40, "firstRead"),(0x80, "secondRead"),(0x100, "secondaryAlignment"),(0x200, "qualFail"),(0x400, "pcrDuplicate")]

	for intCode, name in flagCodes:
		if intSAMflag & intCode:
			dictSAMflag[name] = True
		else:
			dictSAMflag[name] = False

	return(dictSAMflag)



#this function filters a list of SAM dicts representing one insert for pseudo proper paired ends (both reads mapping to same ref)
#assumes that one reference represents one gene, as multiple alignments to the same reference are being kept (as unique mappers)
#this function does the heavy lifting. I.e. most complicated function in this script, maybe it can be simplyfied
#Change for Insert splitting
#checked, but still some open issues
def filterInsert_uniq(listInsertSAMdicts):

	insertName = ""
	listInsertSAMdicts_firstRead = []
	listInsertSAMdicts_secondRead = []
	listInsertInfos = []

	listInsertSAMdicts_filtered_paired = []

	mapped = False
	mappedUnique = False
	setFirstReadsRefs = set()
	setSecondReadsRefs = set()

	dictFirstReadsRef2alignmentScore = defaultdict(int)
	dictSecondReadsRef2alignmentScore = defaultdict(int)
	dictRef2alignmentScore = defaultdict(int)

	counter = 0
	#this loop divides all alignments of an insert into forward and reverse reads
	#if neither the firstRead nor the secondRead is set, the read is assigned to the forward reads
	#forward read is the default so both inserts or single ended reads will be stored as such
	unsignedReadFound = False
	signedReadFound = False
	for dictSAMline in listInsertSAMdicts:
		dictSAMflag = dictSAMline["dictSAMflag"]
		currInsert = dictSAMline["insertName"]
		alignmentScore = dictSAMline["alignmentScore"]

		if (insertName == ""):
			insertName = currInsert
		if (insertName != currInsert):
			sys.stderr.write("Warning: filterInsert received alignments from different inserts. Please check your code. Violating inserts: " + insertName + " " + currInsert+"\n")

		if (dictSAMflag["firstRead"]):
			listInsertSAMdicts_firstRead.append(dictSAMline)
			signedReadFound = True
		elif(dictSAMflag["secondRead"]):
			listInsertSAMdicts_secondRead.append(dictSAMline)
			signedReadFound = True
		else:
			listInsertSAMdicts_firstRead.append(dictSAMline)
			unsignedReadFound = True

	if (signedReadFound and unsignedReadFound):
		sys.stderr.write("Warning: Found signed (forward and reverse) and unsigned (not forward and reverse) reads for following insert: " + insertName)


	#these loops get the best score for each reference
	for dictSAMline in listInsertSAMdicts_firstRead:
		alignmentScore = dictSAMline["alignmentScore"]
		currRefname = dictSAMline["refname"]
		setFirstReadsRefs.add(currRefname)
		dictFirstReadsRef2alignmentScore[currRefname] = max(alignmentScore, dictFirstReadsRef2alignmentScore[currRefname])

	for dictSAMline in listInsertSAMdicts_secondRead:
		alignmentScore = dictSAMline["alignmentScore"]
		currRefname = dictSAMline["refname"]
		setSecondReadsRefs.add(currRefname)
		dictSecondReadsRef2alignmentScore[currRefname] = max(alignmentScore, dictSecondReadsRef2alignmentScore[currRefname])


	#calc paired alignmentScores
	allRefs = list(setFirstReadsRefs.union(setSecondReadsRefs))
	for ref in allRefs:
		dictRef2alignmentScore[ref] = dictFirstReadsRef2alignmentScore[ref] + dictSecondReadsRef2alignmentScore[ref]

	#find alignments with highest pairedScore
	listInsertSAMdicts_besthit_firstRead = filterInsert_getHighestScoringHit(listInsertSAMdicts_firstRead, dictRef2alignmentScore)
	listInsertSAMdicts_besthit_secondRead = filterInsert_getHighestScoringHit(listInsertSAMdicts_secondRead, dictRef2alignmentScore)

	#prefered alignments are those that for which both pairs map on the same reference, only genes with overlap to at least one gene are used
	(preferedSAMdicts_firstRead, preferedRefs_firstRead, normalSAMdicts_firstRead) = filterInsert_getPreferedAlignments(listInsertSAMdicts_besthit_firstRead, setSecondReadsRefs)
	(preferedSAMdicts_secondRead, preferedRefs_secondRead, normalSAMdicts_secondRead) = filterInsert_getPreferedAlignments(listInsertSAMdicts_besthit_secondRead, setFirstReadsRefs)


	#we only care about uniq
	listInsertSAMdicts_uniq = []
	commonPreferedRefs = list(preferedRefs_firstRead.intersection(preferedRefs_secondRead))
	if (len(commonPreferedRefs) == 1):
		uniqRef = commonPreferedRefs[0]
		for dictSAMline in preferedSAMdicts_firstRead:
			currRefname = dictSAMline["refname"]
			if (currRefname == uniqRef):
				listInsertSAMdicts_uniq.append(dictSAMline)
		for dictSAMline in preferedSAMdicts_secondRead:
			currRefname = dictSAMline["refname"]
			if (currRefname == uniqRef):
				listInsertSAMdicts_uniq.append(dictSAMline)

	elif(len(commonPreferedRefs) == 0 and (len(preferedSAMdicts_firstRead) >= 1 or (len(preferedSAMdicts_secondRead) >= 1))):
		if (len(list(preferedRefs_firstRead))== 1):
			listInsertSAMdicts_uniq += preferedSAMdicts_firstRead
		if (len(list(preferedRefs_secondRead))== 1):
			listInsertSAMdicts_uniq += preferedSAMdicts_secondRead

	elif(len(commonPreferedRefs) == 0 and (len(normalSAMdicts_firstRead) >= 1 or (len(normalSAMdicts_secondRead) >= 1))):
		if (len(list(normalSAMdicts_firstRead))== 1):
			listInsertSAMdicts_uniq  += normalSAMdicts_firstRead
		if (len(list(normalSAMdicts_secondRead))== 1):
			listInsertSAMdicts_uniq  += normalSAMdicts_secondRead


	return(listInsertSAMdicts_uniq)



#function to parse cigar strings for a number of values
#checked
def parseCigar(cigarString):
	dictCigarCounts = defaultdict(int)

	listTypes = [i for i in cigarString if not str(i).isdigit()]

	#inTable = ["M","I","D","S","H","P","N","X","="]
	#outTable = ["<>","<>","<>","<>","<>","<>","<>","<>","<>"]
	inTable = "MIDSHPNX="
	outTable = "_"*len(inTable)
	if python_version() == 2:
		translationTable = string.maketrans(inTable, outTable)
		cigarString = cigarString.encode('utf-8')
	else:
		translationTable = str.maketrans(inTable, outTable)


	temp = cigarString.translate(translationTable)
	listNumbers = temp.split("_")

	for position, type in enumerate(listTypes):
		count = listNumbers[position]
		dictCigarCounts[type] += int(count)

	seqLength = 0.0
	alignLength = 0.0
	alignCigarMismatches = 0.0
	alignLengthOnRef = 0.0
	for type in dictCigarCounts:
		if type in "MISH=X":
			seqLength += dictCigarCounts[type]
		if type in "MI=X":
			alignLength += dictCigarCounts[type]
		if type in "IDX":
			alignCigarMismatches += dictCigarCounts[type]
		if type in "MD=X":
			alignLengthOnRef += dictCigarCounts[type]

	return(listTypes, listNumbers, seqLength, alignLength, alignCigarMismatches, alignLengthOnRef)


#please note that samLines is a handle/iteratable
#function that runs the bwa output parsing and calls all the filtering steps
#checked
def parseBWA_SAMoutput(samLines, geneLocationFileName):
	boolGeneBased = True
	dictReference2geneLocation = getReferenceDict(geneLocationFileName, boolGeneBased)

	listInsertSAMdicts = []
	intMappedInserts = 0
	filtered_sam_lines = list()

	prevInsert = ''
	for strSamline in samLines:
		# we dont have lines that starts with @
		#if (strSamline.startswith('@')):
		#	filtered_sam_lines.append(strSamline)
		#	continue

		strSamline = strSamline.strip()

		dictSAMline = parseSamLine(strSamline)
		dictSAMflag = parseSAMflag(dictSAMline["samflag"])
		alignmentScore = parseSAMAlignmentScore(dictSAMline["saminfo"])

		#find insert name
		currInsert = dictSAMline["readname"]
		if (currInsert.endswith('/1') and not dictSAMflag["secondRead"]):
			currInsert = currInsert.rstrip('/1')
			dictSAMflag["firstRead"] = True
			dictSAMflag["secondRead"] = False
		elif (currInsert.endswith('/2') and not dictSAMflag["firstRead"]):
			currInsert = currInsert.rstrip('/2')
			dictSAMflag["secondRead"] = True
			dictSAMflag["firstRead"] = False

		dictSAMline["insertName"] = currInsert
		dictSAMline["dictSAMflag"] = dictSAMflag
		dictSAMline["alignmentScore"] = alignmentScore

		#process the previous Insert
		if (currInsert != prevInsert and prevInsert != ''):
			listInsertSAMdicts_uniq = filterInsert_uniq(listInsertSAMdicts)

			for curr_dictSamline in listInsertSAMdicts_uniq:
				strSamline = curr_dictSamline["samline"]
				filtered_sam_lines.append((strSamline + "\n"))

			listInsertSAMdicts = []

		if (not dictSAMflag["unmapped"]):
			calculateOverlap(dictSAMline, dictReference2geneLocation)
			listInsertSAMdicts.append(dictSAMline)

			prevInsert = currInsert

	listInsertSAMdicts_uniq = filterInsert_uniq(listInsertSAMdicts)
	for curr_dictSamline in listInsertSAMdicts_uniq:
		strSamline = curr_dictSamline["samline"]
		filtered_sam_lines.append((strSamline + "\n"))

	return filtered_sam_lines
