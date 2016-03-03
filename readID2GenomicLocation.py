#!/zfs/datastore0/software/v_envs/py333/bin/python3.3
# Author: Ioannis Moustakas
# Title: Input a number of Read IDs and return their position on the genome. Use 

import pysam
import sys
 
# Save bam File from argv 
genomeAlignFile =str(sys.argv[1])
annotationFile =str(sys.argv[2])

# Store readIDs from STDIN in an array
readIDs=[]
for line in sys.stdin:
    # chomp the last character - new line    
    line = line[:-1]       
    # If line not empty
    if line: 
        readIDs.append(line)


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
        startPos=listOfAtr[1]
        stopPos=listOfAtr[2]
        geneID=listOfAtr[3]
        
        # build a dictionary chromosome => [start, stop, geneID]
        if chrom not in annotation:
            annotation[chrom]=[[startPos, stopPos, geneID]]   
        else:
            annotation[chrom].append([startPos, stopPos, geneID])     
    return(annotation)


# take chromosome and postition and anotation dictionary and return geneID(s)        
def annotateLocation(chrom, pos, annotation):
    geneIDsList=[]
    chromosomeAnnot=annotation[chrom]
    for i in range(len(chromosomeAnnot)):
        annotElement = chromosomeAnnot[i]
        if int(annotElement[0])<pos and int(annotElement[1])>pos:
            geneIDsList.append(annotElement[2])
    if len(geneIDsList)==0:
        geneIDsList=["NonAnnotatedRegion"]
    return(geneIDsList)


def getReadLocation(genomeAlignFile, annotation):
    samFile=pysam.AlignmentFile(genomeAlignFile, 'rb')
    # object to write a bam file
    matchedReads = pysam.AlignmentFile("matchedReads.bam", "wb", template=samFile)
    # readID to chromosome, position
    readIDToAttr={}
    readToRecord={}
    for read in samFile.fetch(until_eof=True):
        if not read.is_unmapped:
            alignmentStart=read.reference_start+1
            cigarTuple=read.cigartuples
            readID=read.query_name
            readIDSplit=readID.split('/')[0]
            mappedOn=samFile.getrname(read.reference_id)
            readIDToAttr[readIDSplit]=[mappedOn, alignmentStart, cigarTuple]
            # Save the record in a dictionary             
            readToRecord[readIDSplit]=read            
            
        else:
            readID=read.query_name
            readIDSplit=readID.split('/')[0]
            readIDToAttr[readIDSplit]=["*", "0", [()]]            
            
    # get the attributes of readIDs 
    chromToLocations={}   
    unMappedReadIDs=[]
    geneIDs=[]
    records=[]     
    for readID in readIDs:
        try: 
            chromosome=readIDToAttr[readID][0]
            location=readIDToAttr[readID][1]
            cigarTuples=readIDToAttr[readID][2]
            records.append(readToRecord[readID])
            # iterate over tuples and print the skipped reference length            
            for tple in cigarTuples:
                if tple[0]==3:
                    skipedRefLength=tple[1]
                    #print(skipedRefLength)
                    
            if chromosome!="*":
                geneID=annotateLocation(chromosome, location, annotation)
                geneIDs.append(geneID)
                if chromosome not in chromToLocations:
                    #print(chromosome, location)
                    chromToLocations[chromosome]=[location]   
                    #print(chromToLocations)
                else:
                    chromToLocations[chromosome].append(location)                   
            else:
                unMappedReadIDs.append(readID)
        except KeyError:
            print(readID + " not in alignment file" )
    
    for key in chromToLocations:
        fivePrime=min(chromToLocations[key])
        threePrime=max(chromToLocations[key])
        regionSpan=threePrime-fivePrime
        numberOfReads=len(chromToLocations[key])      
        
        print("In region {0}:{1}-{2}, spanning {4}nts there are {3} Reads".format(key, fivePrime, threePrime, numberOfReads, regionSpan))
    
    print("\nMoreover, {0} read(s) were left unmapped:".format(len(unMappedReadIDs)))
    print(unMappedReadIDs)
    
    # Flatten the gene IDs list
    geneIDsList=[gene for sublist in geneIDs for gene in sublist]
    countReadsInGenes={gene:geneIDsList.count(gene) for gene in geneIDsList}
    
    for gene in countReadsInGenes:
        print("{1} Reads Mapping On: {0}".format(gene, countReadsInGenes[gene]))
    for record in records:
        matchedReads.write(record)    
    samFile.close()
    matchedReads.close()
    pysam.sort("matchedReads.bam", "matchedReads.sorted")    
    pysam.index("matchedReads.sorted.bam")
    
annotation=getAnnotation(annotationFile)    
getReadLocation(genomeAlignFile, annotation)