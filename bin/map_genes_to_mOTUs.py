#!/usr/bin/env python

from __future__ import division
import os
import sys
import argparse
from collections import defaultdict
import string
import shlex
import time
import subprocess
import tempfile
import shutil

#function that detect the python version
def python_version():
    if(sys.version_info >= (3,0,0)):
        return(3)
    else:
        return(2)


# ------------------------------------------------------------------------------
# function to check if a specific tool exists
def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

#function to print a python dictionary to a tab delimited file
#checked
def printDictToFile(dictData, header, outfileName, return_dictionary):
    if return_dictionary:
        res_dict = dict()
        for key in dictData:
            value = dictData[key]
            res_dict[key] = value
        return res_dict

    if outfileName!="":
        #outfile = open(outfileName,"w")
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(outfile.name, 0o644)
        if (header):
            outfile.write(header + "\n")

        for key in dictData:
            value = dictData[key]
            outfile.write(str(key) + "\t" + str(value) + "\n")

        # flush,sync and close
        try:
            outfile.flush()
            os.fsync(outfile.fileno())
            outfile.close()
        except:
            sys.stderr.write("[E::main] Error: failed to save the profile\n")
            sys.exit(1)

        try:
            #os.rename(outfile.name,outfileName) # atomic operation
            shutil.move(outfile.name,outfileName) #It is not atomic if the files are on different filsystems.
        except:
            sys.stderr.write("[E::calc_motu] Error: failed to save the intermediate mgc table\n")
            sys.stderr.write("[E::calc_motu] you can find the file here:\n"+outfile.name+"\n")
            sys.exit(1)

    else: # print to stdout
        if (header):
            print(header)

        for key in dictData:
            value = dictData[key]
            print(str(key) + "\t" + str(value))


#parse gene location from a file formatted: GeneName \t ContigName \t GeneStart \t GeneEnd
#checked
def getReferenceDict(geneLocationFileName):

    try:
        dictReference2geneLocation = defaultdict(list)
        geneLocationFile = open(geneLocationFileName, "r")

        for strLine in geneLocationFile:
            strLine = strLine.rstrip('\n')
            arrLine = strLine.split('\t')

            strGeneName = arrLine[0]
            strContigName = arrLine[1]  # not used
            strGeneStart = int(arrLine[2])
            strGeneEnd = int(arrLine[3])

            geneInfoTuple = (strGeneName, strGeneStart, strGeneEnd)

            dictReference2geneLocation[strGeneName].append(geneInfoTuple)
    except:
        sys.stderr.write("[E::calc_mgc] Error loading file: "+geneLocationFileName+"\n[E::calc_mgc] Try to download again the motus profiler\n\n")
        sys.exit(1)

    return(dictReference2geneLocation)


#parse mOTU or any other category for the genes in the reference. Counts will be for these categories
#format: GeneName \t category/mOTU
#checked
def getGene2mOTUdict(gene2mOTUfileName):
    dictGene2mOTUs = defaultdict(str)

    try:
        gene2mOTUfile = open(gene2mOTUfileName, "r")
        for strLine in gene2mOTUfile:
            strLine = strLine.rstrip('\n')
            arrLine = strLine.split('\t')

            strGeneName = arrLine[0]
            str_mOTU = arrLine[3]

            dictGene2mOTUs[strGeneName] = str_mOTU
    except:
        sys.stderr.write("[E::calc_mgc] Error loading file: "+gene2mOTUfileName+"\n[E::calc_mgc] Try to download again the motus profiler\n\n")
        sys.exit(1)
    return(dictGene2mOTUs)


#generic function to parse 2 column file with identifier \t value (value is int and parsed as such)
#format: GeneName \t category/mOTU
#checked
def parse2columnFile_int(infileName):
    dictIn = defaultdict(int)

    try:
        infile = open(infileName, "r")

        for strLine in infile:
            strLine = strLine.rstrip('\n')
            arrLine = strLine.split('\t')


            strColumn_1_entry = arrLine[0]
            strColumn_2_entry = int(arrLine[1])

            dictIn[strColumn_1_entry] = strColumn_2_entry
    except:
        sys.stderr.write("[E::calc_mgc] Error loading file: "+infileName+"\n[E::calc_mgc] Try to download again the motus profiler\n\n")
        sys.exit(1)


    return(dictIn)


#function to read in read in SAM or BAM formatted mapping file with initial filtering by msamtools
#checked
def readSAMfile(strSAMfilename, msamPercID, msamminLength, msam_script, msamOverlap):

    try:
        from subprocess import DEVNULL
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')

    #we can understand if it is SAM or BAM reading the first line:
    # in the new version of samtools it is detected automatically if we have a bam or sam file
    isBam = False
    try:
        location = open(strSAMfilename,'r')
    except:
        sys.stderr.write("[E::calc_mgc] Error: failed to open "+strSAMfilename+"\n")
        sys.exit(1)
    try:
        first_line = location.readline()
    except UnicodeDecodeError:
        # is bam
        isBam = True
    location.close()
    # here we keep it for compatibility with old versions of samtools
    samtoolsCMD = "samtools view -Sh -F 4 " + strSAMfilename
    if (isBam):
        samtoolsCMD = "samtools view -h -F 4 " + strSAMfilename

    #sys.stderr.write(samtoolsCMD+"\n")
    samtools_popenCMD = shlex.split(samtoolsCMD)

    if is_tool("samtools"):
        # execute samtools view 2 times, first to check errors and second to save output. Cannot find better solution
        samtools_cmd_t = subprocess.Popen(samtools_popenCMD,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout_s,stderr_s = samtools_cmd_t.communicate()
        stderr_samtool = list()
        for k in stderr_s:
            stderr_samtool.append(k)
        if len(stderr_samtool) != 0:
            sys.stderr.write("[E::calc_mgc] Error from samtools view for file "+strSAMfilename+":\n")
            subprocess.call(samtools_popenCMD,stdout=DEVNULL)
            sys.exit(1)
        samtools_cmd = subprocess.Popen(samtools_popenCMD,stdout=subprocess.PIPE,)
    else:
        sys.stderr.write("[E::calc_mgc] Error: samtools is not in the path. Cannot load the files\n")
        sys.exit(1)


    msamtoolsCMD = "python "+msam_script+" "+str(msamPercID/100)+" "+str(msamminLength)+" "+str(msamOverlap) + " None"

    #sys.stderr.write(msamtoolsCMD+"\n")
    msamtools_popenCMD = shlex.split(msamtoolsCMD)
    msamtools_cmd = subprocess.Popen(msamtools_popenCMD,stdin=samtools_cmd.stdout,stdout=subprocess.PIPE,)

    output = msamtools_cmd.stdout
    return(output)

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


#function to parse mismatches from a sam misc info field
#checked
def getMismatchesSamInfo(strSamInfo):

    listSamInfo = strSamInfo.split("\t")
    mismatches = -1
    mismatchesNM = -1
    mismatchesXM = -1
    for addOn in listSamInfo:
        if (addOn.startswith("NM:i:")):
            misMatchString = addOn.split(":")[2]
            mismatchesNM = int(misMatchString)
        if (addOn.startswith("XM:i:")):
            misMatchString = addOn.split(":")[2]
            mismatchesXM = int(misMatchString)

    if (mismatchesXM != -1):
        mismatches = mismatchesXM
    elif (mismatchesNM != -1):
        mismatches = mismatchesNM
    else:
        mismatches = -1

    return(mismatches)


# this function checks if the ends of an alignment are clipped and if this clipping is within the gene region
#checked
def checkClippedEnds(listCigarTypes):

    unclippedStart = True
    unclippedEnd = True

    if (listCigarTypes[0] == "S" or listCigarTypes[0] == "H"):
        unclippedStart = False
#        if (alignStart < geneStart):
#            unclippedStart = True

    if (listCigarTypes[-1] == "S" or listCigarTypes[-1] == "H"):
        unclippedEnd = False
#        if (alignEnd > geneEnd):
#            unclippedEnd = True

    endsOK = unclippedStart and unclippedEnd
    return(endsOK)


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


#parse one line of SAM alignment to a dictionary
#checked
#469616.genomes.fasta:1:1:1:91#0/1       16      metaMG0013354.COG0541   98      255     75M     *       0       0       TTTTTAGAATGGAAATTTTCCACCCTTTCCACCTTTAAATCCACCAAATGATGGTAATGAAGGTATTTTTCCACC     hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh        AS:i:1000       NM:i:0
def parseSamLine(strSamline):

    #print(strSamline)

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


#filter alignment (single sam lines(in dict format) according to perc id and alignment length)
#accepts 2 different min lengths: 1 for clipped alignments (local alignment) and one for fully aligned reads
#checked
def filterAlignment(dictSAMline, minClippedAlignLength, minAlignLength, minPercID):

    keepAlignment = False
    alignPercID = 0.0

    (listCigarTypes, listNumbers, seqLength, alignLength, alignCigarMismatches, alignLengthOnRef) = parseCigar(dictSAMline["cigar"])
    alignStart = int(dictSAMline["alignstart"])
    alignEnd = alignStart + int(alignLengthOnRef) - 1  # not used

    intMissmatches = getMismatchesSamInfo(dictSAMline["saminfo"])
    if (intMissmatches == -1):
        intMissmatches = alignCigarMismatches
        #noXM_noNM = True

    endsOK = checkClippedEnds(listCigarTypes)

    if (not endsOK):
        minAlignLenRead = minClippedAlignLength
    else:
        minAlignLenRead = minAlignLength

    if (alignLength >= minAlignLenRead):
        alignPercID = 100.0*((alignLength - intMissmatches)/float(alignLength))
        if (alignPercID >= minPercID):
            keepAlignment = True

    return(keepAlignment)


#adds fields "overlapCount" and "overlap" to dictSAMline
#checked
def calculateOverlap(dictSAMline, dictReference2geneLocation):

    (listCigarTypes, listNumbers, seqLength, alignLength, alignCigarMismatches, alignLengthOnRef) = parseCigar(dictSAMline["cigar"])
    alignStart = int(dictSAMline["alignstart"])
    alignEnd = alignStart + int(alignLengthOnRef) - 1

    dictSAMline["overlapCount"] = []

    listRefGenes = dictReference2geneLocation[dictSAMline["refname"]]
    for geneTuple in listRefGenes:
        (geneName, geneStart, geneEnd) = geneTuple

        overlap, overlapCount = checkAlignGeneOverlap(alignStart, alignEnd, geneStart, geneEnd)

        if (overlap and overlapCount > 0):
            dictSAMline["overlap"] = True
            dictSAMline["overlapCount"].append((geneName,overlapCount))

        elif(overlap and overlapCount == 0):
            sys.stderr.write("[W::calc_mgc] Warning: Overlap of 0?\n")
            dictSAMline["overlap"] = False
        else:
            dictSAMline["overlap"] = False


# start of filterInsert helper functions
def filterInsert_getBestAlignmentPer_mOTU(listInsertSAMdicts, dictGene2mOTUs, verbose):
    #first read
    dictBestClusterHit_score = defaultdict(int)
    dictBestClusterHit_overlap = defaultdict(int)
    dictBestClusterHit_hits = defaultdict(list)

    for dictSAMline in listInsertSAMdicts:
        currRefname = dictSAMline["refname"]

        mOTUname = "no_mOTU"
        if (currRefname in dictGene2mOTUs):
            mOTUname = dictGene2mOTUs[currRefname]
        else:
            if verbose>=2: sys.stderr.write("[W::calc_mgc] Warning: "+str(currRefname) + " has no mOTU information\n")

        if (mOTUname == "no_mOTU"):
            mOTUname = currRefname

        alignmentScore = dictSAMline["alignmentScore"]
        dictSAMline["refcluster"] = mOTUname
        dictSAMline["refname_gene"] = currRefname
        dictSAMline["refname"] = mOTUname

        currOverlapCount = 0
        overlapCounts = dictSAMline["overlapCount"]
        for (ref, overlapCount_ref) in overlapCounts:
            if (ref == currRefname and currOverlapCount < overlapCount_ref):
                currOverlapCount = overlapCount_ref

        #this looks for best alignment to this specific mOTU
        if (alignmentScore > dictBestClusterHit_score[mOTUname]):
            dictBestClusterHit_score[mOTUname] = alignmentScore
            dictBestClusterHit_hits[mOTUname] = []
            dictBestClusterHit_hits[mOTUname].append(dictSAMline)
            dictBestClusterHit_overlap[mOTUname] = currOverlapCount

        elif (alignmentScore == dictBestClusterHit_score[mOTUname] and currOverlapCount > dictBestClusterHit_overlap[mOTUname]):
            dictBestClusterHit_hits[mOTUname] = []
            dictBestClusterHit_hits[mOTUname].append(dictSAMline)
            dictBestClusterHit_overlap[mOTUname] = currOverlapCount

    #this could be streamlined
    listInsertSAMdicts = []
    for mOTUname in dictBestClusterHit_hits:
        currentListDictSAMlines = dictBestClusterHit_hits[mOTUname]
        for dictSAMline in currentListDictSAMlines:
            listInsertSAMdicts.append(dictSAMline)

    return(listInsertSAMdicts)


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



#this function filters a list of SAM dicts representing one insert for pseudo proper paired ends (both reads mapping to same ref)
#assumes that one reference represents one gene, as multiple alignments to the same reference are being kept (as unique mappers)
#this function does the heavy lifting. I.e. most complicated function in this script, maybe it can be simplyfied
#Change for Insert splitting
#checked, but still some open issues
def filterInsert(listInsertSAMdicts, dictGene2mOTUs,verbose):

    insertName = ""
    listInsertSAMdicts_firstRead = []
    listInsertSAMdicts_secondRead = []
    listInsertInfos = []

    dictGene2basecount_insert_paired = defaultdict(int)
    listInsertSAMdicts_filtered_paired = []

    mapped = False
    mappedUnique = False
    setFirstReadsRefs = set()
    setSecondReadsRefs = set()

    dictFirstReadsRef2alignmentScore = defaultdict(int)
    dictSecondReadsRef2alignmentScore = defaultdict(int)
    dictRef2alignmentScore = defaultdict(int)

    #this loop divides all alignments of an insert into forward and reverse reads
    #if neither the firstRead nor the secondRead is set, the read is assigned to the forward reads
    #forward read is the default so inserts or single ended reads will be stored as such
    unsignedReadFound = False
    signedReadFound = False
    for dictSAMline in listInsertSAMdicts:
        dictSAMflag = dictSAMline["dictSAMflag"]
        currInsert = dictSAMline["insertName"]
        alignmentScore = dictSAMline["alignmentScore"]

        if (insertName == ""):
            insertName = currInsert
        if (insertName != currInsert):
            if verbose>=2: sys.stderr.write("[W::calc_mgc] Warning: filterInsert received alignments from different inserts. Please check your code. Violating inserts: " + insertName + " " + currInsert + "\n")

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
        if verbose>=2: sys.stderr.write("[W::calc_mgc] Warning: Found signed (forward and reverse) and unsigned (not forward and reverse) reads for following insert: " + insertName +"\n")

    #execute this loop to establish preferential treatment of inserts for which both reads map to the same mOTU/cluster
    #this keeps only the strongest hit per mOTU and replaces the refname with the mOTU name

    #first read
    listInsertSAMdicts_firstRead = filterInsert_getBestAlignmentPer_mOTU(listInsertSAMdicts_firstRead, dictGene2mOTUs,verbose)
    #second read
    listInsertSAMdicts_secondRead = filterInsert_getBestAlignmentPer_mOTU(listInsertSAMdicts_secondRead, dictGene2mOTUs,verbose)


    #these loops get the best score for each reference
    for dictSAMline in listInsertSAMdicts_firstRead:
        alignmentScore = dictSAMline["alignmentScore"]
        currRefname = dictSAMline["refname"]
        #currRefCluster = dictSAMline["refcluster"]
        hasOverlap = False
        for (geneName, overlapCount) in dictSAMline["overlapCount"]:
            if (overlapCount > 0):
                hasOverlap = True
        if (hasOverlap):
            setFirstReadsRefs.add(currRefname)
            dictFirstReadsRef2alignmentScore[currRefname] = max(alignmentScore, dictFirstReadsRef2alignmentScore[currRefname])

    for dictSAMline in listInsertSAMdicts_secondRead:
        alignmentScore = dictSAMline["alignmentScore"]
        currRefname = dictSAMline["refname"]
        hasOverlap = False
        for (geneName, overlapCount) in dictSAMline["overlapCount"]:
            if (overlapCount > 0):
                hasOverlap = True
        if (hasOverlap):
            setSecondReadsRefs.add(currRefname)  ## -> only if we have overlap with the gene (not counting the padded)
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

    #either we have paired mapping or not
    type = ''
    if ((len(preferedSAMdicts_firstRead) >= 1) or (len(preferedSAMdicts_secondRead) >= 1)):
        ## for now not used
        if (len(preferedSAMdicts_firstRead) >= 1):
            listInsertSAMdicts_filtered_paired += preferedSAMdicts_firstRead
        if (len(preferedSAMdicts_secondRead) >= 1):
            listInsertSAMdicts_filtered_paired += preferedSAMdicts_secondRead
        type = 'p'
    else:
        if (len(normalSAMdicts_firstRead) >= 1):
            listInsertSAMdicts_filtered_paired += normalSAMdicts_firstRead
        if (len(normalSAMdicts_secondRead) >= 1):
            listInsertSAMdicts_filtered_paired += normalSAMdicts_secondRead
        type = 's'

    mapped = False
    mappedUnique = False

    #process listInsertSAMdicts_filtered_firstRead_besthit, which includes both forward and reverse reads in case of inserts (both reads) mapping properly to one reference
    if (len(listInsertSAMdicts_filtered_paired) >= 1): # if is unique
        setReferences = set()
        for dictSAMline in listInsertSAMdicts_filtered_paired:
            currRefname = dictSAMline["refname"]
            for (geneName, overlapCount) in dictSAMline["overlapCount"]:
                if (overlapCount > 0):
                    setReferences.add(currRefname)
                    dictGene2basecount_insert_paired[geneName] += overlapCount
                    mapped = True
        if (len(setReferences) == 1):
            mappedUnique = True
        if (len(setReferences) > 1):
            mappedUnique = False
        if (mapped):
            insertInfo = (mapped, mappedUnique, dictGene2basecount_insert_paired, listInsertSAMdicts_filtered_paired, type)
            listInsertInfos.append(insertInfo)

    return(listInsertInfos)


#please note that samLines is a handle/iteratable
#function that runs the bwa output parsing and calls all the filtering steps
#checked
def parseBWA_SAMoutput(samLines, dictGene2counts, dictGene2basecount, dictReference2geneLocation, listMultipleMappers, minAlignLength, minClippedAlignLength, minPercID, dictGene2mOTUs, profile_mode,verbose,version_information_map_read):

    #listMultipleMappers = []
    listInsertSAMdicts = []

    intMappedInserts = 0

    prevInsert = ''

    count0 = 0 # total
    count1 = 0 # insert
    count2 = 0 # filtered insert
    count3 = 0 # unique
    count4 = 0 # multiple
    count5 = 0 # not used
    count6 = 0 # not used

    for strSamline in samLines:
        if not (profile_mode):
            strSamline = strSamline.decode('ascii')
        if (strSamline.startswith('@SQ') or strSamline.startswith('@HD') or strSamline.startswith('@RG') or strSamline.startswith('@CO')):
            continue
        if (strSamline.startswith('@PG')):
            version_information_map_read[0] = strSamline[14:].rstrip()
            continue

        strSamline = strSamline.strip()

        count0 += 1

        dictSAMline = parseSamLine(strSamline)
        if dictSAMline is None:
            sys.stderr.write("[calc_mgc] Warning. Skip line: "+strSamline+"\n")
            continue
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
            count1 += 1

            #filterInserts, if there are no references hit by both reads (forward and reverse), split the Insert into 2 independent reads/fake inserts
            #listInsertInfos = filterInsert(listInsertSAMdicts)
            listInsertInfos = filterInsert(listInsertSAMdicts, dictGene2mOTUs,verbose)

            for insertInfo in listInsertInfos:
                (mapped, mappedUnique, dictGene2basecount_insert, listInsertSAMdicts_filtered, type) = insertInfo

                count2 += 1
                #process uniquely mapping inserts (easy and straightforward)
                if (mappedUnique):
                    #geneName = dictGene2basecount_insert.keys()[0]
                    intMappedInserts += 1
                    numberOfGenesHit = len(dictGene2basecount_insert.keys())
                    for geneName in dictGene2basecount_insert:
                        dictGene2counts[geneName] += (1/float(numberOfGenesHit))
                        dictGene2basecount[geneName] += dictGene2basecount_insert[geneName]
                    for curr_dictSamline in listInsertSAMdicts_filtered:
                        strSamline = curr_dictSamline["samline"]
                    count3 += 1

                #process inserts that map to multiple references
                elif(mapped):
                    listMultipleMappers_insert = []
                    for geneName in dictGene2basecount_insert:
                        tupleGeneBaseCount = (geneName, dictGene2basecount_insert[geneName])
                        listMultipleMappers_insert.append(tupleGeneBaseCount)
                    listMultipleMappers.append((listMultipleMappers_insert, type))
                    intMappedInserts += 1

                    count4 += 1
                else:
                    count5 += 1

                if(mapped or mappedUnique):
                    for curr_dictSamline in listInsertSAMdicts_filtered:
                        strSamline = curr_dictSamline["samline"]

            listInsertSAMdicts = []

        if (not dictSAMflag["unmapped"]):
            keepAlignment = filterAlignment(dictSAMline, minClippedAlignLength, minAlignLength, minPercID)

            if (keepAlignment):
                calculateOverlap(dictSAMline, dictReference2geneLocation)
                listInsertSAMdicts.append(dictSAMline)

            prevInsert = currInsert
                #count6 += 1

    #last Insert
    #listInsertInfos = filterInsert(listInsertSAMdicts)
    listInsertInfos = filterInsert(listInsertSAMdicts, dictGene2mOTUs,verbose)
    for insertInfo in listInsertInfos:
        (mapped, mappedUnique, dictGene2basecount_insert, listInsertSAMdicts_filtered, type) = insertInfo

        count2 += 1
        if (mappedUnique):
            numberOfGenesHit = len(dictGene2basecount_insert.keys())
            for geneName in dictGene2basecount_insert:
                dictGene2counts[geneName] += (1/float(numberOfGenesHit))
                dictGene2basecount[geneName] += dictGene2basecount_insert[geneName]
                intMappedInserts += 1
            for dictSamline in listInsertSAMdicts_filtered:
                strSamline = dictSamline["samline"]
            count3 += 1

        elif(mapped):
            listMultipleMappers_insert = []
            for geneName in dictGene2basecount_insert:
                tupleGeneBaseCount = (geneName, dictGene2basecount_insert[geneName])
                listMultipleMappers_insert.append(tupleGeneBaseCount)
            listMultipleMappers.append((listMultipleMappers_insert, type))
            intMappedInserts += 1

            count4 += 1
        else:
            count5 += 1

        for dictSamline in listInsertSAMdicts_filtered:
            strSamline = dictSamline["samline"]


    #sys.stderr.write("Total SAM lines:\t" + str(count0)+'\n')
    #sys.stderr.write("Inserts:\t" + str(count1)+'\n')
    #sys.stderr.write("Filtered Inserts:\t" + str(count2)+'\n')
    #sys.stderr.write("Inserts mapped unique:\t" + str(count3)+'\n')
    #sys.stderr.write("Inserts mapped multiple:\t" + str(count4)+'\n')
    #print("Inserts unmapped :\t" + str(count5))
    #print("SAMlines filtered: " + str(count6))


def get_mOTU_abundances(dictUniqueInsertCounts, dictUniqueBaseCounts, listMultipleMappers, dictGene2mOTUs, dictGene2Lengths, nonUniqueMultThreshold, winnerThreshold, loserThreshold, sampleName, output, type_output,profile_mode,return_dictionary,verbose):

    uniqueInsertCount = 0
    totalInsertCount = 0

    ignoreMMInsertCount = 0

    dictmOTU_insert_rawCounts = defaultdict(float)
    dictmOTU_insert_coverage = defaultdict(float)
    dictmOTU_bases_rawCounts = defaultdict(float)
    dictmOTU_bases_coverage = defaultdict(float)

    dictmOTU_insert_rawCounts_uniq = defaultdict(float)
    dictmOTU_insert_coverage_uniq = defaultdict(float)
    dictmOTU_bases_rawCounts_uniq = defaultdict(float)
    dictmOTU_bases_coverage_uniq = defaultdict(float)

    totalCoverage_insert_scaled = defaultdict(float)
    totalCoverage_bases_scaled = defaultdict(float)

    mmCount_s = 0
    mmCount_p = 0
    mmCount_sp = 0

    #winnerThreshold = 0.95
    #loserThreshold = 0.01
    #nonUniqueMultThreshold = 3

    #get unique counts
    #this is straight forward
    for geneName in dictUniqueInsertCounts:
        geneLength = dictGene2Lengths[geneName]

        mOTUname = "no_mOTU"
        if (geneName in dictGene2mOTUs):
            mOTUname = dictGene2mOTUs[geneName]
        else:
            sys.stderr.write("[W::calc_mgc] Warning: "+geneName+' not in the map gens->mOTUs\n')

        currInsertCount = dictUniqueInsertCounts[geneName]
        currInsertCoverage = float(currInsertCount)/float(geneLength)
        totalInsertCoverage = currInsertCoverage # this variable is not used?
        dictmOTU_insert_rawCounts[mOTUname] += currInsertCount
        dictmOTU_insert_coverage[mOTUname] += currInsertCoverage
        dictmOTU_insert_rawCounts_uniq[mOTUname] += currInsertCount
        dictmOTU_insert_coverage_uniq[mOTUname] += currInsertCoverage

        currBaseCount = dictUniqueBaseCounts[geneName]
        currBaseCoverage = float(currBaseCount)/float(geneLength)
        totalBaseCoverage = currBaseCoverage # this variable is not used?
        dictmOTU_bases_rawCounts[mOTUname] += currBaseCount
        dictmOTU_bases_coverage[mOTUname] += currBaseCoverage
        dictmOTU_bases_rawCounts_uniq[mOTUname] += currBaseCount
        dictmOTU_bases_coverage_uniq[mOTUname] += currBaseCoverage

        uniqueInsertCount += currInsertCount
        totalInsertCount += currInsertCount

    #split multiple mappers one insert at a time
    for (listMultipleMappers_insert, type) in listMultipleMappers:
        setmOTUnames = set()
        listmOTUs = []
        listUnique_mOTU_inserts_rawCounts = []
        listUnique_mOTU_inserts_coverage = []
        listUnique_mOTU_bases_rawCounts = []
        listUnique_mOTU_bases_coverage = []

        dict_curr_mOTU_insert_count = defaultdict(int)
        dict_curr_mOTU_insert_coverage = defaultdict(int)
        dict_curr_mOTU_base_count = defaultdict(int)
        dict_curr_mOTU_base_coverage = defaultdict(int)

        for (geneName, baseCount) in listMultipleMappers_insert:
            geneLength = dictGene2Lengths[geneName]
            mOTUname = "no_mOTU"
            if (geneName in dictGene2mOTUs):
                mOTUname = dictGene2mOTUs[geneName]
            else:
                sys.stderr.write("[W::calc_mgc] Warning: "+geneName+' not in the map genes -> mOTU\n')
            setmOTUnames.add(mOTUname)

            dict_curr_mOTU_insert_count[mOTUname] += 1
            dict_curr_mOTU_insert_coverage[mOTUname] += 1.0/float(geneLength)
            dict_curr_mOTU_base_count[mOTUname] += baseCount
            dict_curr_mOTU_base_coverage[mOTUname] += float(baseCount)/float(geneLength)

        listmOTUs = list(setmOTUnames)
        #maps uniquely to one mOTU
        if (len(listmOTUs) == 1):
            dictmOTU_insert_rawCounts[mOTUname] += 1
            dictmOTU_insert_coverage[mOTUname] += dict_curr_mOTU_insert_coverage[mOTUname]/float(dict_curr_mOTU_insert_count[mOTUname])

            currBaseCount = dict_curr_mOTU_base_count[mOTUname]/float(dict_curr_mOTU_insert_count[mOTUname])
            dictmOTU_bases_rawCounts[mOTUname] += currBaseCount
            dictmOTU_bases_coverage[mOTUname] += dict_curr_mOTU_base_coverage[mOTUname]/float(dict_curr_mOTU_insert_count[mOTUname])
            uniqueInsertCount += 1
            totalInsertCount += 1

        elif(len(listmOTUs) > 1):
            totalInsertCount += 1

            #get fractions to from unique mappers
            for mOTUname in listmOTUs:
                listUnique_mOTU_bases_coverage.append(dictmOTU_bases_coverage_uniq[mOTUname])

            sumUnique_coverage_bases = sum(listUnique_mOTU_bases_coverage)

            listmOTUs2 = []
            listmOTUs_winner = []
            dictUnique_mOTUCounts_coverage_bases_proportions = defaultdict(float)
            if (sumUnique_coverage_bases > 0):

                for mOTUname in listmOTUs:
                    proportion = float(dictmOTU_bases_coverage_uniq[mOTUname])/float(sumUnique_coverage_bases)
                    if (proportion >= winnerThreshold):
                        listmOTUs_winner.append(mOTUname)
                    if (proportion >= loserThreshold):
                        listmOTUs2.append(mOTUname)

                if (len(listmOTUs_winner) > 0):
                    listmOTUs = listmOTUs_winner
                elif(len(listmOTUs2) > 0):
                    listmOTUs = listmOTUs2

            #if non of the hit references have a unique mapper and more this insert should be shared by more than 3 references --> ignore it
            elif (sumUnique_coverage_bases == 0):
                if (len(listmOTUs) > nonUniqueMultThreshold):
                    listmOTUs = []
                    ignoreMMInsertCount += 1

            if (len(listmOTUs) == 1):
                uniqueInsertCount += 1
            elif(len(listmOTUs) > 1):
                if (type == 's'):
                    mmCount_s += 1
                if (type == 'p'):
                    mmCount_p += 1
                if (type == 'sp'):
                    mmCount_sp += 1

            for mOTUname in listmOTUs:
                listUnique_mOTU_inserts_rawCounts.append(dictmOTU_insert_rawCounts_uniq[mOTUname])
                listUnique_mOTU_inserts_coverage.append(dictmOTU_insert_coverage_uniq[mOTUname])
                listUnique_mOTU_bases_rawCounts.append(dictmOTU_bases_rawCounts_uniq[mOTUname])
                listUnique_mOTU_bases_coverage.append(dictmOTU_bases_coverage_uniq[mOTUname])

            sumUnique_rawCounts_inserts = sum(listUnique_mOTU_inserts_rawCounts)
            sumUnique_coverage_inserts = sum(listUnique_mOTU_inserts_coverage)
            sumUnique_rawCounts_bases = sum(listUnique_mOTU_bases_rawCounts)
            sumUnique_coverage_bases = sum(listUnique_mOTU_bases_coverage)

            #using base coverage the winner and losers are determined (to be consistent)
            #dictmOTU_bases_coverage_uniq

            dictUnique_mOTUCounts_rawCounts_inserts_proportions = defaultdict(float)
            #how can div by zero happen here. check it out -> not one uniq mapper
            if (sumUnique_rawCounts_inserts > 0):
                proportionSum = 0
                for mOTUname in listmOTUs:
                    proportion = float(dictmOTU_insert_rawCounts_uniq[mOTUname])/float(sumUnique_rawCounts_inserts)
                    dictUnique_mOTUCounts_rawCounts_inserts_proportions[mOTUname] = proportion
                    proportionSum += proportion
            else:
                for mOTUname in listmOTUs:
                    dictUnique_mOTUCounts_rawCounts_inserts_proportions[mOTUname] = 1.0/float(len(listmOTUs))

            dictUnique_mOTUCounts_coverage_inserts_proportions = defaultdict(float)
            if (sumUnique_coverage_inserts > 0):
                proportionSum = 0
                for mOTUname in listmOTUs:
                    proportion = float(dictmOTU_insert_coverage_uniq[mOTUname])/float(sumUnique_coverage_inserts)
                    dictUnique_mOTUCounts_coverage_inserts_proportions[mOTUname] = proportion
                    proportionSum += proportion
            else:
                for mOTUname in listmOTUs:
                    dictUnique_mOTUCounts_coverage_inserts_proportions[mOTUname] = 1.0/float(len(listmOTUs))


            dictUnique_mOTUCounts_rawCounts_bases_proportions = defaultdict(float)
            if (sumUnique_rawCounts_bases > 0):
                proportionSum = 0
                for mOTUname in listmOTUs:
                    proportion = float(dictmOTU_bases_rawCounts_uniq[mOTUname])/float(sumUnique_rawCounts_bases)
                    dictUnique_mOTUCounts_rawCounts_bases_proportions[mOTUname] = proportion
                    proportionSum += proportion
            else:
                for mOTUname in listmOTUs:
                    dictUnique_mOTUCounts_rawCounts_bases_proportions[mOTUname] = 1.0/float(len(listmOTUs))

            dictUnique_mOTUCounts_coverage_bases_proportions = defaultdict(float)
            if (sumUnique_coverage_bases > 0):
                proportionSum = 0
                for mOTUname in listmOTUs:
                    proportion = float(dictmOTU_bases_coverage_uniq[mOTUname])/float(sumUnique_coverage_bases)
                    dictUnique_mOTUCounts_coverage_bases_proportions[mOTUname] = proportion
                    proportionSum += proportion
            else:
                for mOTUname in listmOTUs:
                    dictUnique_mOTUCounts_coverage_bases_proportions[mOTUname] = 1.0/float(len(listmOTUs))


            #apply fractions
            for mOTUname in listmOTUs:

                rawCounts_prop_insert = dictUnique_mOTUCounts_rawCounts_inserts_proportions[mOTUname]
                dictmOTU_insert_rawCounts[mOTUname] += 1.0*rawCounts_prop_insert

                coverage_prop_insert = dictUnique_mOTUCounts_coverage_inserts_proportions[mOTUname]
                dictmOTU_insert_coverage[mOTUname] += (dict_curr_mOTU_insert_coverage[mOTUname]/float(dict_curr_mOTU_insert_count[mOTUname]))*coverage_prop_insert

                rawCounts_prop_bases = dictUnique_mOTUCounts_rawCounts_bases_proportions[mOTUname]
                dictmOTU_bases_rawCounts[mOTUname] += (dict_curr_mOTU_base_count[mOTUname]/float(dict_curr_mOTU_insert_count[mOTUname]))*rawCounts_prop_bases

                coverage_prop_bases = dictUnique_mOTUCounts_coverage_bases_proportions[mOTUname]
                dictmOTU_bases_coverage[mOTUname] += (dict_curr_mOTU_base_coverage[mOTUname]/float(dict_curr_mOTU_insert_count[mOTUname]))*coverage_prop_bases

        else:
            if verbose>=2: sys.stderr.write("[W::calc_mgc] Warning: Could not find any mOTU for multiple mapper: " + str(listMultipleMappers_insert) + '\n')

    totalCount_insert = 0
    totalCoverage_insert = 0.0
    totalCount_bases = 0
    totalCoverage_bases = 0.0

    for mOTUname in dictmOTU_insert_rawCounts:
        totalCount_insert += dictmOTU_insert_rawCounts[mOTUname]
        totalCoverage_insert += dictmOTU_insert_coverage[mOTUname]
        totalCount_bases += dictmOTU_bases_rawCounts[mOTUname]
        totalCoverage_bases += dictmOTU_bases_coverage[mOTUname]

    #scaling coverages into insert count space
    for mOTUname in dictmOTU_insert_coverage:

        totalCoverage_insert_scaled[mOTUname] = (float(dictmOTU_insert_coverage[mOTUname])/float(totalCoverage_insert))*totalCount_insert
        totalCoverage_bases_scaled[mOTUname] = (float(dictmOTU_bases_coverage[mOTUname])/float(totalCoverage_bases))*totalCount_insert

    #Todo: handle unassigned (for general application)
    header = sampleName
    if type_output == 'all':
        outfileName_insert_rawCounts = output+".insert.rawCounts.tsv"
        outfileName_insert_coverage = output+ ".insert.coverage.tsv"
        outfileName_bases_rawCounts = output+ ".bases.rawCounts.tsv"
        outfileName_bases_coverage = output + ".bases.coverage.tsv"
        outfileName_insert_scaled = output + ".insert.scaled.tsv"
        outfileName_bases_scaled = output + ".bases.scaled.tsv"

        printDictToFile(dictmOTU_insert_rawCounts, header, outfileName_insert_rawCounts,profile_mode)
        printDictToFile(dictmOTU_insert_coverage, header, outfileName_insert_coverage,profile_mode)
        printDictToFile(dictmOTU_bases_rawCounts, header, outfileName_bases_rawCounts,profile_mode)
        printDictToFile(dictmOTU_bases_coverage, header, outfileName_bases_coverage,profile_mode)
        printDictToFile(totalCoverage_insert_scaled, header, outfileName_insert_scaled,profile_mode)
        printDictToFile(totalCoverage_bases_scaled, header, outfileName_bases_scaled,profile_mode)

        outfileName_insert_rawCounts_uniq = output+ ".uniq.insert.rawCounts.tsv"
        outfileName_insert_coverage_uniq = output+ ".uniq.insert.coverage.tsv"
        outfileName_bases_rawCounts_uniq = output + ".uniq.bases.rawCounts.tsv"
        outfileName_bases_coverage_uniq = output+ ".uniq.bases.coverage.tsv"

        printDictToFile(dictmOTU_insert_rawCounts_uniq, header, outfileName_insert_rawCounts_uniq,profile_mode)
        printDictToFile(dictmOTU_insert_coverage_uniq, header, outfileName_insert_coverage_uniq,profile_mode)
        printDictToFile(dictmOTU_bases_rawCounts_uniq, header, outfileName_bases_rawCounts_uniq,profile_mode)
        printDictToFile(dictmOTU_bases_coverage_uniq, header, outfileName_bases_coverage_uniq,profile_mode)

        #sys.stderr.write("UniqueCounts: " + str(uniqueInsertCount) + '\n')
        #sys.stderr.write("TotalCounts: " + str(totalInsertCount) + '\n')
        #sys.stderr.write("Ignored multiple mapper without unique hit: " + str(ignoreMMInsertCount) + '\n')

        #sys.stderr.write("Multiple Mappers (both ends mapping to different refs):  " + str(mmCount_s) + '\n')
        #sys.stderr.write("Multiple Mappers (both ends mapping to same ref): " + str(mmCount_p) + '\n')
        #print("mmCount_sp: " + str(mmCount_sp))

    if type_output == 'insert.raw_counts':
        if return_dictionary:
            return printDictToFile(dictmOTU_insert_rawCounts, header, output,return_dictionary)
        else:
            printDictToFile(dictmOTU_insert_rawCounts, header, output,return_dictionary)
    if type_output == 'insert.coverage':
        if return_dictionary:
            return printDictToFile(dictmOTU_insert_coverage, header, output,return_dictionary)
        else:
            printDictToFile(dictmOTU_insert_coverage, header, output,return_dictionary)
    if type_output == 'insert.scaled_counts':
        if return_dictionary:
            return printDictToFile(totalCoverage_insert_scaled, header, output,return_dictionary)
        else:
            printDictToFile(totalCoverage_insert_scaled, header, output,return_dictionary)
    if type_output == 'base.coverage':
        if return_dictionary:
            return printDictToFile(dictmOTU_bases_coverage, header, output,return_dictionary)
        else:
            printDictToFile(dictmOTU_bases_coverage, header, output,return_dictionary)
    if type_output == 'bases.raw_counts':
        if return_dictionary:
            return printDictToFile(dictmOTU_bases_rawCounts, header,output,return_dictionary)
        else:
            printDictToFile(dictmOTU_bases_rawCounts, header,output,return_dictionary)
    if type_output == 'bases.scaled':
        if return_dictionary:
            return printDictToFile(totalCoverage_bases_scaled, header,output,return_dictionary)
        else:
            printDictToFile(totalCoverage_bases_scaled, header,output,return_dictionary)
    if type_output == 'uniq.bases.coverage':
        if return_dictionary:
            return printDictToFile(dictmOTU_bases_coverage_uniq, header,output,return_dictionary)
        else:
            printDictToFile(dictmOTU_bases_coverage_uniq, header,output,return_dictionary)
    if type_output == 'uniq.bases.raw_counts':
        if return_dictionary:
            return printDictToFile(dictmOTU_bases_rawCounts_uniq, header, output,return_dictionary)
        else:
            printDictToFile(dictmOTU_bases_rawCounts_uniq, header, output,return_dictionary)
    if type_output == 'uniq.insert.coverage':
        if return_dictionary:
            return printDictToFile(dictmOTU_insert_coverage_uniq, header, output,return_dictionary)
        else:
            printDictToFile(dictmOTU_insert_coverage_uniq, header, output,return_dictionary)
    if type_output == 'uniq.insert.raw_counts':
        if return_dictionary:
            return printDictToFile(dictmOTU_insert_rawCounts_uniq, header, output,return_dictionary)
        else:
            printDictToFile(dictmOTU_insert_rawCounts_uniq, header, output,return_dictionary)



#functions still needed:
#get_and_check_mOTU_db

def run_mOTUs_v2_mapping(listInputFiles, databaseDir, databasePrefix, sampleName, nonUniqueMultThreshold, winnerThreshold, loserThreshold, minClippedAlignLength, output, msam_script, type_output, verbose, profile_mode, input_sam_file_for_profile, return_dictionary,min_perc_id,min_len_align,min_perc_align):
    if verbose>2: start_time = time.time()
    #constants
    minAlignLength = min_len_align
    minPercID = min_perc_id
    msamOverlap = min_perc_align

    if verbose >= 5: sys.stderr.write("Filter in calc_mgc: MIN_PERC_ID:"+str(min_perc_id)+" MIN_LENGTH_ALIGN: "+str(min_len_align)+" MIN_PERC_COVER: "+str(min_perc_align)+" \n")


    ## check that the input files exists
    if profile_mode == False:
        cont = 0
        for inputFile in listInputFiles:
            cont = cont + 1
            if not os.path.isfile(inputFile):
                if (cont == 1): sys.stderr.write("[E::calc_mgc] Error. "+inputFile+': No such file.\n')
                if (cont > 1): sys.stderr.write("[E::calc_mgc] Error. Input file number "+str(cont)+" ("+inputFile+'): No such file.\n')
                sys.exit(1)
            else:
                if not (inputFile.endswith(".bam") or inputFile.endswith(".sam")):
                    if verbose>=2: sys.stderr.write("[W::calc_mgc] Warning. File "+inputFile+" does not have extension \".bam\" or \".sam\" \n")


    geneLocationFileName = os.path.sep.join([databaseDir,  databasePrefix + ".padded.coords"])
    dictReference2geneLocation = getReferenceDict(geneLocationFileName)

    geneLengthFileName = os.path.sep.join([databaseDir, databasePrefix + ".genes.len"])
    dictGene2Lengths = parse2columnFile_int(geneLengthFileName)

    gene2mOTUfileName = os.path.sep.join([databaseDir, databasePrefix + ".map"])
    dictGene2mOTUs = getGene2mOTUdict(gene2mOTUfileName)

    listMultipleMappers = []
    dictGene2counts = defaultdict(int)
    dictGene2basecount = defaultdict(int)

    # version used in the 1st script
    version_information_map_read = {0:"no_info"}

    if profile_mode is False:
        cont = 0
        all_infor_map_reads = dict()
        for inputFile in listInputFiles:
            samLines = readSAMfile(inputFile, minPercID, minAlignLength, msam_script, msamOverlap)
            parseBWA_SAMoutput(samLines, dictGene2counts, dictGene2basecount, dictReference2geneLocation, listMultipleMappers, minAlignLength, minClippedAlignLength, minPercID, dictGene2mOTUs, profile_mode,verbose,version_information_map_read)
            all_infor_map_reads[cont] = version_information_map_read[0]
            cont = cont + 1
        # check that the first script is consistent

        first_info = all_infor_map_reads[0].split(" | ")
        error_flag_not_correct_header = False
        if len(first_info) != 3:
            error_flag_not_correct_header = True
        try:
            for k in all_infor_map_reads:
                k_info = all_infor_map_reads[k].split(" | ")
                if len(k_info) != 3:
                    error_flag_not_correct_header = True
                if first_info[0] != k_info[0]:
                    if verbose>1: sys.stderr.write(" [W::calc_mgc] Warning: reads mapped with different version of the tool ("+first_info[0]+","+k_info[0]+")\n")
                if first_info[1] != k_info[1]:
                    if verbose>1: sys.stderr.write(" [W::calc_mgc] Warning: reads mapped with different version of the database ("+first_info[1]+","+k_info[1]+")\n")
        except:
            error_flag_not_correct_header = True

        if len(version_information_map_read[0].split(" | ")) != 3:
            error_flag_not_correct_header = True

        if error_flag_not_correct_header:
            if verbose>1: sys.stderr.write(" [W::calc_mgc] Warning: couldn't find any header in the sam/bam file(s)\n")
            version_information_map_read[0] = "map_tax unknown | gene database: unknown | 100"
        else:
            version_information_map_read[0] = all_infor_map_reads[0]
    else: # if profile mode is true
        samLines = input_sam_file_for_profile
        parseBWA_SAMoutput(samLines, dictGene2counts, dictGene2basecount, dictReference2geneLocation, listMultipleMappers, minAlignLength, minClippedAlignLength, minPercID, dictGene2mOTUs, profile_mode,verbose,version_information_map_read)

    if verbose>2:
        if len(listInputFiles) == 1: errstr = "1 sam/bam file"
        if len(listInputFiles) > 1: errstr = str(len(listInputFiles))+ " files"
        sys.stderr.write(" [calc_mgc](parse " + errstr + ") "+ str("{0:.2f}".format(time.time() - start_time))+" sec\n")

    if verbose>2: start_time = time.time()
    if not(return_dictionary):
        get_mOTU_abundances(dictGene2counts, dictGene2basecount, listMultipleMappers, dictGene2mOTUs, dictGene2Lengths, nonUniqueMultThreshold, winnerThreshold, loserThreshold, sampleName, output, type_output, profile_mode, return_dictionary,verbose)
    else:
        dict_temp = get_mOTU_abundances(dictGene2counts, dictGene2basecount, listMultipleMappers, dictGene2mOTUs, dictGene2Lengths, nonUniqueMultThreshold, winnerThreshold, loserThreshold, sampleName, output, type_output, profile_mode, return_dictionary,verbose)
    if verbose>2: sys.stderr.write(" [calc_mgc](get mgc abundances) " + str("{0:.2f}".format(time.time() - start_time))+" sec\n")

    if (return_dictionary): return version_information_map_read[0],dict_temp


def main(argv=None):

    parser = argparse.ArgumentParser(description='This program calculates mOTU abundances for one sample', add_help = True)
    parser.add_argument('--inputFile', '-i', action="store", dest='listInputFiles', default="", help='name of input file(s); sam or bam formatted files. If it is a list: insert all files separated by a comma')
    parser.add_argument('--databaseDir', '-d', action="store", dest='databaseDir', default=".", help='name of database directory')
    parser.add_argument('--databasePrefix', '-p', action="store", dest='databasePrefix', default="mOTU.v2b", help='name of database (prefix)')
    parser.add_argument('--output', '-o', action="store", dest='output', default="", help='name of output file, if not present is stdout')
    parser.add_argument('--type', '-y', action="store", dest='type_output', default='bases.coverage', help='type of output that you want to print',choices=['base.coverage','bases.raw_counts','bases.scaled','insert.coverage','insert.raw_counts','insert.scaled_counts','uniq.bases.coverage','uniq.bases.raw_counts','uniq.insert.coverage','uniq.insert.raw_counts','all'])
    parser.add_argument('--sampleName', '-sn', action="store", dest='sampleName', default="", help='sample name for the current mapping')
    parser.add_argument('--multThreshold', '-m',action="store", dest='multThreshold', type=int, default=3, help='Threshold that regulates when a read should be discarded if it maps to over x mOTUs that don\'t have any unique alignments.')
    parser.add_argument('--winnerThreshold', '-w',action="store", dest='winnerThreshold', type=float, default=0.95, help='Threshold that regulates when a multiple mapper should be completely assigned to the highest unique abundance (relative abundance among all references hit by this insert).')
    parser.add_argument('--loserThreshold', '-l',action="store", dest='loserThreshold', type=float, default=0.01, help='Threshold that regulates when a multiple mapper should be removed if the unique abundance is less than this (relative abundance among all references hit by this insert).')
    parser.add_argument('--minClippedAlignLength', '-a',action="store", dest='minClippedAlignLength', type=int, default=60, help='Minimum alignment length for clipped reads (local alignments).')
    parser.add_argument('-v', action='store', type=int, default=3, dest='verbose', help='Verbose levels')
    args = parser.parse_args()

    if (args.listInputFiles == ""):
        sys.stderr.write("[E::calc_mgc] Error: Please provide at least one input file.\n")
        sys.exit(1)

    # find the position of msamtools_python.py
    # we expect to find it in the same directory as runBWA.py
    path_mOTUs = os.path.realpath(__file__)
    path_array = path_mOTUs.split("/")
    relative_path = "/".join(path_array[0:-1])
    msam_script = relative_path+"/msamtools_python.py"

    # if sample name is not provided, then we select the name that was given to the input file
    if args.sampleName == '':
        args.sampleName = 'temp'


    # check that with all there is also output - for now you cannot select 'all'
    if args.type_output == 'all' and args.output == "":
        sys.stderr.write("[E::calc_mgc] Error: --output not present. \nWhen --type is equal to 'all', you must insert also --output\n")
        sys.exit(1)

    # convert --inputFile that is a string to a list
    listInputFiles = args.listInputFiles.split(",")

    profile_mode = False # when using motu profile, this is set to True
    input_sam_file_for_profile = "" # use for motu profile
    return_dictionary = False

    run_mOTUs_v2_mapping(listInputFiles, args.databaseDir, args.databasePrefix, args.sampleName, args.multThreshold, args.winnerThreshold, args.loserThreshold, args.minClippedAlignLength, args.output, msam_script, args.type_output,args.verbose,profile_mode,input_sam_file_for_profile,return_dictionary,97,45,45)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
