#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 10:59:53 2021

@author: laurinerolland
"""

import argparse
import os
from os.path import basename
import re
import sys

# Creation of the regular expressions
exprTranscriptLine = r'^\w+\s\w+\stranscript\s(?P<start>\d+)\s(?P<stop>\d+).+\s[+|-]\s.+\s.*gene_id\s\"(?P<gene_id>[^\"]+)\";\s.*transcript_id\s\"(?P<transcript_id>[^\"]+)\";\s.*gene_type\s\"(?P<gene_type>[^\"]+)\";\s.*transcript_type\s"(?P<transcript_type>[^\"]+)\";\s.*level\s(?P<level>\d);\s.*transcript_support_level\s\"(?P<transcript_support_level>[^\"]+)\";.*$'
transcriptLine = re.compile(exprTranscriptLine)

exprExonLine = r"^\w+\s\w+\sexon\s(?P<start>\d+)\s(?P<stop>\d+).+\s[+|-]\s.+\s.*gene_id\s\"(?P<gene_id>[^\"]+)\";\s.*transcript_id\s\"(?P<transcript_id>[^\"]+)\";\s.*gene_type\s\"(?P<gene_type>[^\"]+)\";\s.*transcript_type\s\"(?P<transcript_type>[^\"]+)\";\s.*exon_id\s\"(?P<exon_id>[^\"]+)\";\s.*level\s(?P<level>\d);\s.*transcript_support_level\s\"(?P<transcript_support_level>[^\"]+)\";.*$"
exonLine = re.compile(exprExonLine)


def name_file (annotated_file): #extraction of the file name
    file_name_ext = basename(annotated_file) #to get the file name
    file_name = os.path.splitext(file_name_ext)[0]
    return (file_name)


def longest_coding_transcript(annotated_file, gene_type, exon, level_list, tsl_list):#exon is a boolean
    try:
        file=open(annotated_file,"r") #the annotated file is opened for reading
    except:
        print ("The file {} cannot be found".format(annotated_file))
    else:
        line=file.readline()#the file is read line by line
        dico_length_transcript={} #{gene_id:{transcript_id:length}}
        dico_transcript={} #{transcript_id:line}
        dico_transcript_exon={} #{transcript_id:[exon_id,exon_id]}
        dico_exon={} #{transcript_id:{exon_id:line}}
        while line != "": #We browse all the lines of the file  
            transcript_line = transcriptLine.search(line) #regular expression
            if transcript_line:
                geneType = transcript_line.group('gene_type')
                transcript_type = transcript_line.group('transcript_type')
                level = int(transcript_line.group('level')) 
                tsl_str = transcript_line.group('transcript_support_level')
                if tsl_str != "NA":
                    tsl = int(tsl_str)
                    if geneType == '{}'.format(gene_type) and transcript_type == '{}'.format(gene_type) and level in level_list and tsl in tsl_list:#if the gene type matches the one entered by the user 
                        gene_id = transcript_line.group('gene_id')
                        transcript_id = transcript_line.group('transcript_id')
                        if gene_id not in dico_length_transcript:#if the gene_id had not been encountered yet
                            dico_length_transcript[gene_id]={}#initialisation {gene_id:{}}
                        if transcript_id not in dico_transcript:#if the transcript_id had not been encountered yet
                            dico_transcript[transcript_id]=line #{transcript_id:line}
                            dico_length_transcript[gene_id][transcript_id]=0 #{gene_id:{transcript_id:0}}
                            if exon:
                                dico_transcript_exon[transcript_id]=[] #{transcript_id:[]}
                                dico_exon[transcript_id]={} #initialisation {transcript_id:{}}
            #print(dico_length_transcript)
            exon_line = exonLine.search(line) #regular expression
            if exon_line: 
                gene_type1 = exon_line.group('gene_type')
                transcript_type1 = exon_line.group('transcript_type')
                level1 = int(exon_line.group('level'))
                tsl1_str = exon_line.group('transcript_support_level')
                if tsl1_str != "NA":
                    tsl1 = int(tsl1_str)
                    if level_list: 
                        if gene_type1 == '{}'.format(gene_type) and transcript_type1 == '{}'.format(gene_type) and level1 in level_list and tsl1 in tsl_list:#if the gene type matches the one entered by the user 
                            gene_id1 = exon_line.group('gene_id')
                            transcript_id1 = exon_line.group('transcript_id')
                            exon_id = exon_line.group('exon_id')
                            start = exon_line.group('start')
                            stop = exon_line.group('stop')
                            exon_length = int(stop) - int(start)#length of the exon
                            dico_length_transcript[gene_id1][transcript_id1] += exon_length #the length of the exons corresponding to the transcript is added
                            if exon:
                                dico_exon[transcript_id1][exon_id]=line #{ transcript_id: {exon_id : line_exon} }
                                (dico_transcript_exon[transcript_id1]).append(exon_id) #{transcript_id:[exon_id,exon_id]}
            line=file.readline()
        return(dico_length_transcript, dico_transcript, dico_exon, dico_transcript_exon)


def filtered_file(filtered_file, dictionary_length, dictionary_transcript, exon_dictionary, dictionary_transcript_exon, exon, transcript):
    f=open(filtered_file,"w")#new file
    for k,v in dictionary_length.items(): #the dictionary {gene_id:{transcript_id:length}} is browsed k=gene_id and v=corresponding values {transcript_id:length}
        max_length=0#initialisation
        for k1,v1 in v.items(): #k1=transcript_id and v1=length
            if v1>max_length:
                max_length=v1 # max_length takes the length v1
                transcript_max=k1 #trasncrit_id of length max
        if transcript:
            f.write("{}".format(dictionary_transcript[transcript_max]))#the line corresponding to the longest coding transcript is written into the file
        if exon:#if the user has requested to display the exons
            exon_list=dictionary_transcript_exon[transcript_max] #list of exons of the longest coding transcript : {transcript_id:[exon_id,exon_id]} --> [exon_id,exon_id]
            for exonID in exon_list:
                f.write("{}".format(exon_dictionary[transcript_max][exonID]))#the line corresponding to the exons of the longest coding transcript is written into the file
    f.close()

def parse_arguments():
    '''creation of options 
    add -h to the command line to display the different options'''
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", help="bed file that will be used")
    parser.add_argument("--output", help="output file")
    parser.add_argument("--type", help="gene type, by default protein_coding")
    parser.add_argument("--gene", help="enter yes if you want to display gene lines")
    parser.add_argument("--transcript", help="enter yes if you want to display transcript lines")
    parser.add_argument("--exon", help="enter yes if you want to display exon lines")
    parser.add_argument("--level", help="level of annotation. Separated by ,")
    parser.add_argument("--TSL", help="transcript support level. Separated by ,")
    #parser.print_help() #to display the help
    return parser.parse_args()

def main():
    args=parse_arguments()
    if args.annotations:
        gene_annotations = args.annotations #annotation file
        #type : protein_coding, lncRNA, unprocessed_pseudogene
        if args.type:
            transcript_type= args.type
        else:
            transcript_type="protein_coding" #by default the type is protein_coding
        #output file name
        if args.output:
            output_file = args.output
        else:
            file_name_annotated = name_file(gene_annotations)
            output_file = "{}.filtered.gtf".format(file_name_annotated)#creation of the default file name from the name of the annotated file
        #display of transcripts
        if args.transcript=="yes":
            transcript=True
        else:
            transcript=False
        #display of exons
        if args.exon=="yes":
            exon = True
        else:
            exon  =False
        #level
        if args.level:
            level_list_str = (args.level).split(",")
            level_list = list(map(int, level_list_str)) # list of integer 
        else:
            level_list = [1] # by defaut level = 1
        #tsl
        if args.TSL:
            tsl_list_str = (args.TSL).split(",")
            tsl_list = list(map(int, tsl_list_str)) # list of integer 
        else:
            tsl_list = [1]
        #creation of dictionaries    
        (dictionary_transcript_length,dictionary_line_transcript, dictionary_exon, dictionary_exon_transcript)=longest_coding_transcript(gene_annotations, transcript_type, exon, level_list, tsl_list )
        #creation of the filtered file
        filtered_file(output_file, dictionary_transcript_length, dictionary_line_transcript, dictionary_exon, dictionary_exon_transcript, exon, transcript)
        print("The file have been created")
    else:
        print("The annotated file has not been entered as input")
        sys.exit()
        
main()