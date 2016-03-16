#!/zfs/datastore0/software/v_envs/py333/bin/python3.3
# Author: Ioannis Moustakas
# Title: Input a number of Read IDs and return their position on the genome. Use 

import pysam
import sys
import bisect
 
 
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
    

# A function to read the transcriptome alignment file
def getTranscriptome(transcAlignFile):
    samFile=pysam.AlignmentFile(transcAlignFile, 'r')
    # readID to chromosome, position
    readInstList=[]
    for read in samFile.fetch(until_eof=True):
        if not read.is_unmapped:
            readID=read.query_name
            readIDSplit=readID.split('/')[0]            
            readInst=Alignment(readIDSplit)
            readInst.cigarTuple(read.cigartuples)
            readInst.referenceID(samFile.getrname(read.reference_id))
            readInst.misMatches(read.get_tag('NM'))
            readInst.mapQual(read.mapping_quality)
            readInst.readLength(read.query_length)
            readInstList.append(readInst)            

        else:
            readID=read.query_name
            readIDSplit=readID.split('/')[0]
            readInst=Alignment(readIDSplit)
    return(readInstList)
    
geneReads=getTranscriptome(transcAlignFile)  
#read=readIDs[0]
#print(read.readID, "\n", read.referenceID, "\n", read.misMatches, "\n", read.mapQual, "\n",  read.readLength, "\n", read.cigarTuple)
#
#print(len(readIDs))
#geneReads=[]
#for read in readIDs:
#    if read.referenceID=="ENST00000604952":
#        geneReads.append(read)
    

""" A function to read the simplified annotation file and return a dictionary 
    chromosome => [start, stop, geneID] for all gene IDs
""" 
def getAnnotation(annotationFile):
    fileHandle=open(annotationFile, 'r')
    annotation={}
    for line in fileHandle:
        line=line[:-1] 
        listOfAtr=line.split('\t')
        chrom=listOfAtr[0]
        startPos=int(listOfAtr[1])
        stopPos=int(listOfAtr[2])
        geneID=listOfAtr[3]
        
        # build a dictionary chromosome => [start, stop, geneID]
        if chrom not in annotation:
            annotation[chrom]=[[startPos, stopPos, geneID]]   
        else:
            annotation[chrom].append([startPos, stopPos, geneID])     
    # sort for later use with bisect        
    for chromosome in annotation:
        annotation[chromosome].sort()
    return(annotation)


# take chromosome and postition and anotation dictionary and return geneID(s)        
# use bisect to make it faster, given the huge amount of comparissons
def annotateLocation(chrom, pos, annotation):

    geneIDsList=[]
    chromosomeAnnot=annotation[chrom]
    for i in range(len(chromosomeAnnot)):
        annotElement = chromosomeAnnot[i]
        if int(annotElement[0])<=pos and int(annotElement[1])>=pos:
            geneIDsList.append(annotElement[2])
    if len(geneIDsList)==0:
        geneIDsList=["NonAnnotatedRegion"]
    return(geneIDsList)    
    


def getReadLocation(genomeAlignFile, annotation, readIDs):
    samFile=pysam.AlignmentFile(genomeAlignFile, 'rb')
    
    # readID to chromosome, position
    readIDToAttr=[]
    #readToRecord={}
    for read in samFile.fetch(until_eof=True):
        readID=read.query_name
        readIDSplit=readID.split('/')[0]
        if not read.is_unmapped:
            alignmentStart=read.reference_start+1
            cigarTuple=read.cigartuples
            mappedOn=samFile.getrname(read.reference_id)
            readIDToAttr.append([readIDSplit, mappedOn, alignmentStart, cigarTuple])
            # Save the record in a dictionary             
            #readToRecord[readIDSplit]=read            
            
        else:
            readIDToAttr.append([readIDSplit, "*", "0", [()]])
            
    # get the attributes of readIDs 
    matching=0
    matchingNMTag=0
    matchingLength=0
    nonMatching=0  
    nonMatchingNMTag=0
    nonMatchingLength=0
    
    # sort the list. Should use the 1st element of the nested lists
    readIDToAttr.sort()
    readIDsOnly=[read[0] for read in readIDToAttr]
    cnt=0    
    for read in readIDs:
        cnt+=1
        try: 
            idx=index(readIDsOnly, read.readID)
            chromosome=readIDToAttr[idx][1]
            location=readIDToAttr[idx][2]
#            cigarTuples=readIDToAttr[read.readID][2]
            #records.append(readToRecord[read.readID])
            # iterate over tuples and print the skipped reference length            
#            for tple in cigarTuples:
#                if tple[0]==3:
#                    skipedRefLength=tple[1]
#                    #print(skipedRefLength)
                    
            if chromosome!="*":
                genomicAlignmentGeneID=annotateLocation(chromosome, location, annotation)
                #print(genomicAlignmentGeneID)
#                print("refID: {0}, gene ID: {1}, chromosome: {2}, location: {3} ".format(read.referenceID, genomicAlignmentGeneID, chromosome, location ))                
                
                if read.referenceID in genomicAlignmentGeneID: 
                    matching+=1
                    matchingNMTag+=read.misMatches
                    matchingLength+=read.readLength
                else:
                    nonMatching+=1
                    nonMatchingNMTag+=read.misMatches
                    nonMatchingLength+=read.readLength
            if cnt%10000==0:
                print("Processing read n. ", str(cnt))
                sys.exit()
        # index function might raise ValueError            
        except ValueError:
            print(read.readID + " not in alignment file")

    print("reads matching", matching)
    print("reads non matching", nonMatching)
    print("nonMatching/Total: ", nonMatching/(matching+nonMatching) )
    print("reads Matching NM/readLength: ", matchingNMTag/matchingLength)
    print("reads nonMatching NM/readLength: ", nonMatchingNMTag/nonMatchingLength)
        
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
    
annotation=getAnnotation(annotationFile)
getReadLocation(genomeAlignFile, annotation, geneReads)