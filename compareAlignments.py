#!/zfs/datastore0/software/v_envs/py333/bin/python3.3
# Author: Ioannis Moustakas
# Title: Input a number of Read IDs and return their position on the genome. Use 

import pysam
import sys
import bisect
import re
import gc
from pprint import pprint 
import time


 
# Save bam File from argv 
genomeAlignFile =str(sys.argv[1])
transcAlignFile = str(sys.argv[2])
annotationFile =str(sys.argv[3])

# Store readIDs from STDIN in an array
#readIDs=[]
#for line in sys.stdin:
#    # chomp the last character - new line    
#    line = line[:-1]       
#    # If line not empty
#    if line: 
#        readIDs.append(line)


# Go object oriented
### A read class to store all indivisual Alignments
class Alignment:
    def __init__(self, readID):
        self.readID=readID
    
    def alignmentStart(self, alignmentStart):
        self.alignmentStart=alignmentStart
    
    def referenceID(self, referenceID):
        self.referenceID=referenceID
        
    def misMatches(self, misMatches):
        self.misMatches=misMatches
        
    def mapQual(self, mapQual):
        self.mapQual=mapQual    
        
    def readLength(self, readLength):
        self.readLength=readLength    

    def cigarTuple(self, cigarTuple):
        self.cigarTuple=cigarTuple 

# search into a list 
def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

# A function to read an alignment file and return a list of alignment (read) instances
def getAlignment(alignmentFile):
    if transcAlignFile.endswith(".sam"):
        samFile=pysam.AlignmentFile(alignmentFile, 'r')
    elif transcAlignFile.endswith(".bam"):
        samFile=pysam.AlignmentFile(alignmentFile, 'rb')
    else: 
        print("Unknown type of file. Expecting .bam or .sam extensions only")
        sys.exit()
        
    # readID to chromosome, position
    readInstList=[]
    for read in samFile.fetch(until_eof=True):
        if not read.is_unmapped:
            readID=read.query_name
            readIDSplit=readID.split('/')[0]            
            readInst=Alignment(readIDSplit)
            readInst.cigarTuple(read.cigartuples)
            readInst.referenceID(samFile.getrname(read.reference_id))
            # Not all alignment files have NM tag, catch the key error             
            try:            
                readInst.misMatches(read.get_tag('NM'))
            except KeyError:
                pass
            readInst.mapQual(read.mapping_quality)
            readInst.readLength(read.query_length)
            readInst.alignmentStart(read.reference_start)
            readInstList.append(readInst)            

        else:
            readID=read.query_name
            readIDSplit=readID.split('/')[0]
            readInst=Alignment(readIDSplit)
            
    return(readInstList)
    
geneReads=getAlignment(transcAlignFile)  
print("Transcriptome Alignment Loaded")
#read=readIDs[0]
#print(read.readID, "\n", read.referenceID, "\n", read.misMatches, "\n", read.mapQual, "\n",  read.readLength, "\n", read.cigarTuple)
#
#print(len(readIDs))
#geneReads=[]
#for read in readIDs:
#    if read.referenceID=="ENST00000604952":
#        geneReads.append(read)
    

# take chromosome and postition, anotation dictionary and last saved start for looping (j) and return geneID(s)        
def annotateLocation(chrom, pos, annotation, j):
    geneIDsList=[]
    chromosomeAnnot=annotation[chrom]
    first=True
    length=len(chromosomeAnnot)
    # iterate over all annotated features
    for i in range(j, length):
        annotElement = chromosomeAnnot[i]
        if annotElement[0]<=pos and annotElement[1]>=pos:
            geneIDsList.append(annotElement[2])
        # Start positions are ordered. If start > pos, break the loop
        if annotElement[0]>pos:
            break
        # get the index of the first stop position that is larger than pos and save it in j.
        # next iteration start iterating from j, to make the search space smaller.
        # This requires that current pos is always larger than the previous.
#        if first:    
        if annotElement[1]>=pos and first:   
                j=i
                first=False
    # if the lsit is empty, it is a non annotated region
    if len(geneIDsList)==0:
        geneIDsList=["NonAnnotatedRegion"]
    return(geneIDsList, j)    
    


""" A function to read a gff3 file and return a dictionary 
    chromosome => position => [transcriptID(s)]
""" 
def getAnnotation(annotationFile):
    fileHandle=open(annotationFile, 'r')
    annotation={}
    firstLine=True
    chromosomesLen={}
    for line in fileHandle:
        line=line[:-1] 
        #if line starts with #        
        if line[0]=="#":
            if firstLine:
                firstLine=False
                if line!="##gff-version   3":
                    print("Annotation not a GFF3 file")
                    #sys.exit
            # get the length for all chromosomes 
            if line[0:17]=="##sequence-region":
                splitLine=line.split("   ")
                chrom=splitLine[1].split(" ")[0]
                chromLen=int(splitLine[1].split(" ")[2])
                chromosomesLen[chrom]=chromLen
        else:    
            listOfAttr=line.split('\t')
            chrom=listOfAttr[0]
            typeOfFeature=listOfAttr[2]
            startPos=int(listOfAttr[3])
            stopPos=int(listOfAttr[4])
            attributes=listOfAttr[8]
            # if the type of feature is good, build a dictionary chromosome => [start, stop, geneID]
            if re.search("transcript", typeOfFeature):
                featureID=re.findall('ID=transcript:([A-Z0-9]+);', attributes)
                if len(featureID)>0:
                    if chrom not in annotation:
                        annotation[chrom]=[[startPos, stopPos, featureID[0] ]]   
                    else:
                        annotation[chrom].append([startPos, stopPos, featureID[0] ])
    return(annotation)
    
annotation=getAnnotation(annotationFile)
print("Annotation Loaded")

#print(annotation["1"][100])
#print(annotation["1"][13000])
#print(annotation["1"][44000])
#print(annotation["1"][52500])
#print(annotation["1"][55000])
#print(annotation["1"][111000])
#print(annotation["1"][90000])
#print(annotation["1"][145000])
#annotateLocation("1", 13000, annotation)
#annotateLocation("1", 44000, annotation)
#annotateLocation("1", 52500, annotation)
#annotateLocation("1", 55000, annotation)
#annotateLocation("1", 111000, annotation)
#annotateLocation("1", 90000, annotation)



# read the genomic alignment file, calculate and print some relevant statistics
def getReadLocation(genomeAlignFile, annotation, transcrReads):
    genomicAlignment=getAlignment(genomeAlignFile)
    
    
#    # readID to chromosome, position
#    genomicReadIDToAttr=[]
#    #readToRecord={}
#    for read in samFile.fetch(until_eof=True):
#        readID=read.query_name
#        readIDSplit=readID.split('/')[0]
#        if not read.is_unmapped:
#            alignmentStart=read.reference_start+1
#            cigarTuple=read.cigartuples
#            mappedOn=samFile.getrname(read.reference_id)
#            genomicReadIDToAttr.append([readIDSplit, mappedOn, alignmentStart, cigarTuple])
#            # Save the record in a dictionary             
#            #readToRecord[readIDSplit]=read            
#            
#        else:
#            genomicReadIDToAttr.append([readIDSplit, "*", "0", [()]])
            
    # get the attributes of readIDs 
    matching=0
    matchingNMTag=0
    matchingLength=0
    nonMatching=0  
    nonMatchingNMTag=0
    nonMatchingLength=0
    
    # sort the list on chromosome and start possition
    genomicReadIDToAttr = sorted(genomicAlignment, key = lambda x: (x.referenceID, x.alignmentStart))
    genomicReadIDs=[read.readID for read in genomicAlignment]
    # Sort on readID so I can use bisect on the resulting list
    transcrReadsSorted = sorted(transcrReads, key = lambda x: x.readID)
    transcrReadIDs=[read.readID for read in transcrReadsSorted] 
    cnt=0  
    j=0
    lastChrom=genomicReadIDToAttr[0].referenceID
    for read in genomicReadIDToAttr:
        genomicReadID=read.readID
        genomicReadChrom=read.referenceID
        genomicReadStartPos=read.alignmentStart
        # read list is shorted. Every time a new chromosme comes, set j=0
        if lastChrom!=genomicReadChrom:
            j=0
            lastChrom=genomicReadChrom
            print("Chromosome: ", lastChrom)
        
        cnt+=1       

        try:
            idx=index(transcrReadIDs, genomicReadID)
            transcrRead=transcrReadsSorted[idx]
            if genomicReadChrom!="*":
                genomicAlignmentGeneID, j = annotateLocation(genomicReadChrom, genomicReadStartPos, annotation, j)
#                print("refID: {0}, gene ID: {1}, chromosome: {2}, location: {3} ".format(read.referenceID, genomicAlignmentGeneID, chromosome, location ))                
            else:
                genomicAlignmentGeneID=[]
                
            if transcrRead.referenceID in genomicAlignmentGeneID: 
                matching+=1
                matchingNMTag+=transcrRead.misMatches
                matchingLength+=transcrRead.readLength
            else:
                nonMatching+=1
                nonMatchingNMTag+=transcrRead.misMatches
                nonMatchingLength+=transcrRead.readLength
            if cnt%100000==0:
                print("Processing read n. ", str(cnt))
        # index function might raise ValueError            
        except ValueError:
            print(genomicReadID + " not in transcriptome alignment file")

    print("reads matching", matching)
    print("reads non matching", nonMatching)
    print("nonMatching/Total: ", nonMatching/(matching+nonMatching) )
    print("reads Matching NM/readLength: ", matchingNMTag/matchingLength)
    print("reads nonMatching NM/readLength: ", nonMatchingNMTag/nonMatchingLength)
 
getReadLocation(genomeAlignFile, annotation, geneReads)
      
#    for key in chromToLocations:
#        fivePrime=min(chromToLocations[key])
#        threePrime=max(chromToLocations[key])
#        regionSpan=threePrime-fivePrime
#        numberOfReads=len(chromToLocations[key])      
#        
#        print("In region {0}:{1}-{2}, spanning {4}nts there are {3} Reads".format(key, fivePrime, threePrime, numberOfReads, regionSpan))
#    
#    print("\nMoreover, {0} read(s) were left unmapped:".format(len(unMappedReadIDs)))
#    print(unMappedReadIDs)
#    
#    # Flatten the gene IDs list
#    geneIDsList=[gene for sublist in geneIDs for gene in sublist]
#    countReadsInGenes={gene:geneIDsList.count(gene) for gene in geneIDsList}
#    
#    for gene in countReadsInGenes:
#        print("{1} Reads Mapping On: {0}".format(gene, countReadsInGenes[gene]))
#    for record in records:
#        matchedReads.write(record)    
#    samFile.close()
#    matchedReads.close()
#    pysam.sort("matchedReads.bam", "matchedReads.sorted")    
#    pysam.index("matchedReads.sorted.bam")
    



#            cigarTuples=readIDToAttr[read.readID][2]
            #records.append(readToRecord[read.readID])
            # iterate over tuples and print the skipped reference length            
#            for tple in cigarTuples:
#                if tple[0]==3:
#                    skipedRefLength=tple[1]
#                    #print(skipedRefLength)
    # sort to speed up(?) and build the dict chrom => pos => geneID(s)
#    annotatedGenome={}
#    gc.disable()
#    for chromosome in annotation:
#        print("Annotating chromosome: ", chromosome)
#        annotation[chromosome].sort()
#        for line in annotation["1"]:
#            print(line)
#        chromLen=chromosomesLen[chromosome]
##        annotatedGenome[chromosome]={}
#        chromosomeAnnotTable=[]
#        j=0
#        for pos in range(chromLen): 
#            if pos%1000000==0:
#                print("base n:", pos, "on chromosome:", chromosome)
##                print("index: ", j)
#            featureID, j =annotateLocation(pos, annotation[chromosome], j)
#            chromosomeAnnotTable.append(featureID)
#        annotatedGenome[chromosome]=chromosomeAnnotTable
#    gc.enable()