#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 13:29:36 2021

@author: laurinerolland
"""


import argparse
import os
from os.path import basename
import re

# Creation of the regular expressions
exprGeneLine = r'^(?P<seqname>\w+)\s(?P<source>\w+)\s(?P<feature>\w+)\s(?P<start>\d+)\s(?P<end>\d+).+\s(?P<strand>[+-])\s(?P<frame>[\d+|.])\sgene_id\s\"(?P<gene_id>[^\"]+)\";\sgene_type\s\"(?P<gene_type>[^\"]+)\";\sgene_name\s\"(?P<gene_name>[^\"]+)\";\slevel\s(?P<level>\d);.*$'
geneLine = re.compile(exprGeneLine)

exprTranscriptLine = r'^(?P<seqname>\w+)\s(?P<source>\w+)\s(?P<feature>\w+)\s(?P<start>\d+)\s(?P<end>\d+).+\s(?P<strand>[+-])\s(?P<frame>[\d+|.])\sgene_id\s\"(?P<gene_id>[^\"]+)\";\stranscript_id\s\"(?P<transcript_id>[^\"]+)\";\sgene_type\s\"(?P<gene_type>[^\"]+)\";\sgene_name\s\"(?P<gene_name>[^\"]+)\";\stranscript_type\s"(?P<transcript_type>[^\"]+)\";\stranscript_name\s"(?P<transcript_name>[^\"]+)\";\slevel\s(?P<level>\d);\stranscript_support_level\s\"(?P<transcript_support_level>[^\"]+)\";.*$'
transcriptLine = re.compile(exprTranscriptLine)

exprTranscriptCodingLine = r'^(?P<seqname>\w+)\s(?P<source>\w+)\s(?P<feature>\w+)\s(?P<start>\d+)\s(?P<end>\d+).+\s(?P<strand>[+-])\s(?P<frame>[\d+|.])\sgene_id\s\"(?P<gene_id>[^\"]+)\";\stranscript_id\s\"(?P<transcript_id>[^\"]+)\";\sgene_type\s\"(?P<gene_type>[^\"]+)\";\sgene_name\s\"(?P<gene_name>[^\"]+)\";\stranscript_type\s"(?P<transcript_type>[^\"]+)\";\stranscript_name\s"(?P<transcript_name>[^\"]+)\";\slevel\s(?P<level>\d);\sprotein_id\s\"(?P<protein_id>[^\"]+)\";\stranscript_support_level\s\"(?P<transcript_support_level>[^\"]+)\";.*$'
transcriptCodingLine = re.compile(exprTranscriptCodingLine)

exprExonLine = r'^(?P<seqname>\w+)\s(?P<source>\w+)\s(?P<feature>\w+)\s(?P<start>\d+)\s(?P<end>\d+).+\s(?P<strand>[+-])\s(?P<frame>[\d+|.])\sgene_id\s\"(?P<gene_id>[^\"]+)\";\stranscript_id\s\"(?P<transcript_id>[^\"]+)\";\sgene_type\s\"(?P<gene_type>[^\"]+)\";\sgene_name\s\"(?P<gene_name>[^\"]+)\";\stranscript_type\s"(?P<transcript_type>[^\"]+)\";\stranscript_name\s"(?P<transcript_name>[^\"]+)\";\sexon_number\s(?P<exon_number>\d);\sexon_id\s\"(?P<exon_id>[^\"]+)\";\slevel\s(?P<level>\d);\stranscript_support_level\s\"(?P<transcript_support_level>[^\"]+)\";.*$'
exonLine = re.compile(exprExonLine)

exprOtherLine1 = r'^(?P<seqname>\w+)\s(?P<source>\w+)\s(?P<feature>\w+)\s(?P<start>\d+)\s(?P<end>\d+).+\s(?P<strand>[+-])\s(?P<frame>[\d+|.])\sgene_id\s\"(?P<gene_id>[^\"]+)\";\stranscript_id\s\"(?P<transcript_id>[^\"]+)\";\sgene_type\s\"(?P<gene_type>[^\"]+)\";\sgene_name\s\"(?P<gene_name>[^\"]+)\";\stranscript_type\s"(?P<transcript_type>[^\"]+)\";\stranscript_name\s"(?P<transcript_name>[^\"]+)\";\sexon_number\s(?P<exon_number>\d+);\sexon_id\s\"(?P<exon_id>[^\"]+)\";\slevel\s(?P<level>\d);\sprotein_id\s\"(?P<protein_id>[^\"]+)\";.*$'
otherLine1 = re.compile(exprOtherLine1)

exprOtherLine2 = r'^(?P<seqname>\w+)\s(?P<source>\w+)\s(?P<feature>\w+)\s(?P<start>\d+)\s(?P<end>\d+).+\s(?P<strand>[+-])\s(?P<frame>[\d+|.])\sgene_id\s\"(?P<gene_id>[^\"]+)\";\stranscript_id\s\"(?P<transcript_id>[^\"]+)\";\sgene_type\s\"(?P<gene_type>[^\"]+)\";\sgene_name\s\"(?P<gene_name>[^\"]+)\";\stranscript_type\s"(?P<transcript_type>[^\"]+)\";\stranscript_name\s"(?P<transcript_name>[^\"]+)\";\sexon_number\s(?P<exon_number>\d+);\sexon_id\s\"(?P<exon_id>[^\"]+)\";\slevel\s(?P<level>\d);.*$'
otherLine2 = re.compile(exprOtherLine2)

def name_file(annotated_file): #extraction of the file name
    file_name_ext = basename(annotated_file) #to get the file name
    file_name = os.path.splitext(file_name_ext)[0]
    return (file_name)

def info_table(annotation_file):
    try:
        file=open(annotation_file,"r") #the annotated file is opened for reading
    except:
        print ("The file {} cannot be found".format(annotation_file))
    else:
        line=file.readline()#the file is read line by line
        list_table = [["seqname","source","feature","start","end","length","strand","frame","gene_id","transcript_id","gene_type","gene_name","transcript_type","transcript_name","exon_number","exon_id","level","protein_id","TSL"]]
        i=0
        while line != "":
            gene_line = geneLine.search(line)
            transcript_line = transcriptLine.search(line)
            transcript_coding_line = transcriptCodingLine.search(line)
            exon_line = exonLine.search(line)
            other_line1 = otherLine1.search(line)
            other_line2 = otherLine2.search(line)
            info=[]
            if gene_line:
                seqname1 = gene_line.group("seqname")
                source1 = gene_line.group("source")
                feature1 = gene_line.group("feature")
                start1 = int(gene_line.group("start"))
                end1 = int(gene_line.group("end"))
                strand1 = gene_line.group("strand")
                frame1 = gene_line.group("frame")
                gene_id1 = gene_line.group ("gene_id")
                gene_type1 = gene_line.group ("gene_type")
                gene_name1 = gene_line.group ("gene_name")
                level1 = gene_line.group ("level")
                length1= end1 - start1
                info = [seqname1, source1, feature1, start1, end1, length1 , strand1, frame1 , gene_id1 , "NA" , gene_type1, gene_name1,"NA","NA","NA","NA",level1,"NA","NA"]
            elif transcript_line:
                seqname2 = transcript_line.group("seqname")
                source2 = transcript_line.group("source")
                feature2 = transcript_line.group("feature")
                start2 = int(transcript_line.group("start"))
                end2 = int(transcript_line.group("end"))
                strand2 = transcript_line.group("strand")
                frame2 = transcript_line.group("frame")
                gene_id2 = transcript_line.group ("gene_id")
                transcript_id2 = transcript_line.group ("transcript_id")
                gene_type2 = transcript_line.group ("gene_type")
                gene_name2 = transcript_line.group ("gene_name")
                transcript_type2 = transcript_line.group ("transcript_type")
                transcript_name2 = transcript_line.group ("transcript_name")
                level2 = transcript_line.group ("level")
                level2 = transcript_line.group ("level")
                TSL2 = transcript_line.group ("transcript_support_level")
                length2 = end2 - start2
                info = [seqname2, source2, feature2, start2, end2, length2, strand2, frame2, gene_id2, transcript_id2, gene_type2, gene_name2, transcript_type2, transcript_name2, "NA", "NA", level2, "NA", TSL2]
            elif transcript_coding_line:
                seqname6 = transcript_coding_line.group("seqname")
                source6 = transcript_coding_line.group("source")
                feature6 = transcript_coding_line.group("feature")
                start6 = int(transcript_coding_line.group("start"))
                end6 = int(transcript_coding_line.group("end"))
                strand6 = transcript_coding_line.group("strand")
                frame6 = transcript_coding_line.group("frame")
                gene_id6 = transcript_coding_line.group ("gene_id")
                transcript_id6 = transcript_coding_line.group ("transcript_id")
                gene_type6 = transcript_coding_line.group ("gene_type")
                gene_name6 = transcript_coding_line.group ("gene_name")
                transcript_type6 = transcript_coding_line.group ("transcript_type")
                transcript_name6 = transcript_coding_line.group ("transcript_name")
                level6 = transcript_coding_line.group ("level")
                level6 = transcript_coding_line.group ("level")
                protein_id6 = transcript_coding_line.group ("protein_id")
                TSL6 = transcript_coding_line.group ("transcript_support_level")
                length6 = end6 - start6
                info = [seqname6, source6, feature6, start6, end6, length6, strand6, frame6, gene_id6, transcript_id6, gene_type6, gene_name6, transcript_type6, transcript_name6, "NA", "NA", level6, protein_id6, TSL6]
            elif exon_line:
                seqname3 = exon_line.group("seqname")
                source3 = exon_line.group("source")
                feature3 = exon_line.group("feature")
                start3 = int(exon_line.group("start"))
                end3 = int(exon_line.group("end"))
                strand3 = exon_line.group("strand")
                frame3 = exon_line.group("frame")
                gene_id3 = exon_line.group ("gene_id")
                transcript_id3 = exon_line.group ("transcript_id")
                gene_type3 = exon_line.group ("gene_type")
                gene_name3 = exon_line.group ("gene_name")
                transcript_type3 = exon_line.group ("transcript_type")
                transcript_name3 = exon_line.group ("transcript_name")
                exon_number3 = exon_line.group ("exon_number")
                exon_id3 = exon_line.group ("exon_id")
                level3 = exon_line.group ("level")
                TSL3 = exon_line.group ("transcript_support_level")
                length3 = end3 - start3
                info = [seqname3, source3, feature3, start3, end3, length3, strand3, frame3, gene_id3, transcript_id3, gene_type3, gene_name3, transcript_type3, transcript_name3, exon_number3, exon_id3, level3, "NA",TSL3]
            elif other_line1:
                seqname4 = other_line1.group("seqname")
                source4 = other_line1.group("source")
                feature4 = other_line1.group("feature")
                start4 = int(other_line1.group("start"))
                end4 = int(other_line1.group("end"))
                strand4 = other_line1.group("strand")
                frame4 = other_line1.group("frame")
                gene_id4 = other_line1.group ("gene_id")
                transcript_id4 = other_line1.group ("transcript_id")
                gene_type4 = other_line1.group ("gene_type")
                gene_name4 = other_line1.group ("gene_name")
                transcript_type4 = other_line1.group ("transcript_type")
                transcript_name4 = other_line1.group ("transcript_name")
                exon_number4 = other_line1.group ("exon_number")
                exon_id4 = other_line1.group ("exon_id")
                level4 = other_line1.group ("level")   
                protein_id4 = other_line1.group ("protein_id")
                length4 = end4- start4
                info = [seqname4, source4, feature4, start4, end4, length4, strand4, frame4, gene_id4, transcript_id4, gene_type4, gene_name4, transcript_type4, transcript_name4, exon_number4, exon_id4, level4, protein_id4,"NA"]
            elif other_line2:
                seqname5 = other_line2.group("seqname")
                source5 = other_line2.group("source")
                feature5 = other_line2.group("feature")
                start5 = int(other_line2.group("start"))
                end5 = int(other_line2.group("end"))
                strand5 = other_line2.group("strand")
                frame5 = other_line2.group("frame")
                gene_id5 = other_line2.group ("gene_id")
                transcript_id5 = other_line2.group ("transcript_id")
                gene_type5 = other_line2.group ("gene_type")
                gene_name5 = other_line2.group ("gene_name")
                transcript_type5 = other_line2.group ("transcript_type")
                transcript_name5 = other_line2.group ("transcript_name")
                exon_number5 = other_line2.group ("exon_number")
                exon_id5 = other_line2.group ("exon_id")
                level5 = other_line2.group ("level")   
                length5 = end5- start5
                info = [seqname5, source5, feature5, start5, end5, length5, strand5, frame5, gene_id5, transcript_id5, gene_type5, gene_name5, transcript_type5, transcript_name5, exon_number5, exon_id5, level5, "NA","NA"]
            list_table.append(info)
            line=file.readline()
            i+=1
            #print(list_table)
        return (list_table)

def bedfile(table_list, file):
    f=open(file,"w")
    for element in table_list:
        for info in element:
            f.write("{}\t".format(info))
        f.write("\n")  
        
def parse_arguments():
    '''creation of options 
    add -h to the command line to display the different options'''
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", help="bed file that will be used")
    parser.add_argument("--output", help="output file")
    return parser.parse_args()

def main(): 
    args=parse_arguments()
    if args.annotations:
        file_annotation = args.annotations #annotation file
        if args.output:
            output_file = args.output
        else:
            file_name_annotated = name_file(file_annotation)
            output_file = "{}.bed".format(file_name_annotated)#creation of the default file name from the name of the annotated file
        # creation de la liste
        table = info_table(file_annotation)
        print(table)
        # creation du fichier bed
        bedfile(table, output_file)
    print("The file have been created")
        
main()