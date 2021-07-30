#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: laurinerolland
"""
import pybedtools
import argparse
import os
import math
from os.path import basename

def correspondence_matrix_interval(bedFile,matrixFile):
    ''' test the length of the intervals in the bed file and 
        return a booleen which indicates if the intervals are equal or not
        the length of the motif and the list of lines where the length of the intervals is different, 
        if there is no difference an empty list is returned
        We also check that the size of the matrix pattern is equal to the size of the intervals'''
    try:
        bed=open(bedFile,"r") #the bed file is opened for reading
    except:
        print ("The file {} cannot be found".format(bed))
    else:
        bed_name_ext = basename(bedFile) #to get the file name
        bed_name = os.path.splitext(bed_name_ext)[0]
        line=bed.readline().strip() #the file is read line by line
        i=0
        list_error=[]
        while line != "": #We browse all the lines of the file
            content_line=line.split() #the line of the file is transformed into a list, in this way elements 2 and 3 of the list will be the beginning and the end of the pattern position 
            start=content_line[1]
            end=content_line[2]
            if i==0: # if it is the first line
                length_motif=int(end)-int(start) # first value of length_motif
            else:
                if int(end)-int(start)!=length_motif: #the length of the line that is being read is compared with the previous line
                    list_error.append(i+1)
                    #print("Error in motif length at line {}".format(i+1))
                else:
                    length_motif=int(end)-int(start)
            line=bed.readline().strip()
            i+=1
        if list_error != []:
            test_length=False # peut être renommer length_interval
        else:
            test_length=True
    try:
        matrix=open(matrixFile,"r") #the matrix file is opened for reading
    except:
        print ("The file {} cannot be found".format(matrix))
    else:
        line=matrix.readline() #the file is read line by line
        number_line=0
        while line !="":
            number_line+=1
            line=matrix.readline()
        if int(number_line-1)!=length_motif:
            test_motif=False
        else:
            test_motif=True
        return (test_length, length_motif, list_error, test_motif, bed_name) #returns the length of the motif

def percent_gc(fastaFile): #calculate the base composition of a fasta file
    '''return the base composition of a fasta file in a dictionary format 
    {'A': 0.265124441336834, 'T': 0.2648310535578382, 'C': 0.2339301370383187, 'G': 0.23611436806700917}'''
    percent_background={"A":0, "T":0, "C":0, "G":0} #number of each nucleotide observed in the fasta file given as input to the program
    total_nucleotides=0
    try:
        fasta=open(fastaFile,"r")#the fasta file is opened for reading
    except:
        print ("The file {} cannot be found".format(fasta))
    else:
        line_fasta=fasta.readline().strip().upper() #the file is read line by line
        while line_fasta != "": #We browse all the lines of the file
            if line_fasta[0]!=">":
                for nucleotide in line_fasta: #we browse through the nucleotides in the line
                    if nucleotide != "N": 
                        percent_background[nucleotide] = percent_background[nucleotide] + 1
                        total_nucleotides+=1
            line_fasta=fasta.readline().strip().upper()
        for nucleotide in percent_background: #we look at the number of times each nucleotide has been observed and divide it by the total number of nucleotides
            percent_background[nucleotide]=((percent_background[nucleotide])/total_nucleotides)
        #dictionary of percentages for each nucleotide
        #exemple:{'A': 0.265124441336834, 'T': 0.2648310535578382, 'C': 0.2339301370383187, 'G': 0.23611436806700917}
    return (percent_background)

def dictionary_gc(gc): #make a dictionary from the GC percentage input by the user 
    ''' return the base composistion in a dictionary format from the GC composition input by the user '''
    dict_gc={}
    dict_gc["A"]=(1-float(gc))/2
    dict_gc["T"]=(1-float(gc))/2
    dict_gc["C"]=(float(gc))/2
    dict_gc["G"]=(float(gc))/2
    return(dict_gc)

def matrix_to_dictionary(matrixFile,length_motif):
    """converts the file containing the PWM into a dictionary"""
    dictionary={"A":{}, "C":{}, "G":{}, "T":{}}
    try:
        matrix=open(matrixFile,"r") #the bed file is opened for reading
    except:
        print ("The file {} cannot be found".format(matrix))
    else:
        first_line=matrix.readline() # ligne d'entête qu'on ne veut pas
        matrix_line=matrix.readline()
        content_line=matrix_line.split()
        i=1
        while i < length_motif+1:
            dictionary["A"][i]=content_line[0]
            dictionary["C"][i]=content_line[1]
            dictionary["G"][i]=content_line[2]
            dictionary["T"][i]=content_line[3]
            i+=1
            matrix_line=matrix.readline()
            content_line=matrix_line.split()
    return(dictionary)
    
def calcul_score(bedFile, fastaFile, newFile, dictionary_matrix, length, strand, seq_option, gc_freq):
    """score calculation and option to add the sequence column"""
    try:
        bed=open(bedFile,"r") #the bed file is opened for reading
    except:
        print ("The file {} cannot be found".format(bed))
    else:
        f=open(newFile,"w")
        nucleotides_list=list(dictionary_matrix.keys())
        #pybedtools
        a = pybedtools.BedTool(bedFile)
        if strand==True:
            a = a.sequence(fi=fastaFile,s=True)
        else:
            a = a.sequence(fi=fastaFile)
        sequence_fasta=open(a.seqfn) #this variable contains the file obtained with the previous command
        #form of the file:
        #>chr22:10541717-10541727(+)
        #GTAAAAATTA
        #>chr22:10541784-10541794(+)
        #ATAAAGAAAA
        #>chr22:10541810-10541820(-)
        #TATAAAGAAT 
        bed_line=bed.readline().strip()
        line=sequence_fasta.readline() #the file is read line by line
        #écounter_line=0 #number of sequence
        while bed_line != "": #We browse all the lines of the file
            if line[0] == ">": #lines beginning with > do not contain the sequence
                line=sequence_fasta.readline() 
            nucleotide_sequence_line=line.upper()#the sequence is uppercased because it is not uppercased on all lines
            i=0 #this counter is used to browse the sequence, it is initialized to 0
            score=0
            while i<(length):
                for nucleotide in nucleotides_list:
                    if nucleotide_sequence_line[i]==nucleotide:
                        score += float(math.log2((dictionary_matrix[nucleotide][i+1])/(gc_freq[nucleotide])))
                        #score += float(dictionary_matrix[nucleotide][i+1])
                i+=1
                    #counter_line+=1
            if seq_option == True:
                f.write("{} \t {} \t {} \n".format(bed_line,score,nucleotide_sequence_line))
            else:
                f.write("{} \t {} \n".format(bed_line,score))
            score=0
            line=sequence_fasta.readline()
            bed_line=bed.readline().strip()
        bed.close()
        f.close()
        print("The file has been generated")
    
def parse_arguments():
    '''creation of options 
    add -h to the command line to display the different options'''
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", help="bed file that will be used")
    parser.add_argument("--fasta", help="fasta file that will be used")
    parser.add_argument("--matrix", help="fasta file that will be used")
    parser.add_argument("--gc", help="percentage of GC, it can be displayed in different ways : enter -1 to calculate the GC percentage of the fasta file; enter a number between 0 and 1 corresponding to the percentage of GC; by default the percentage of GC is 0.5")
    parser.add_argument("--strand", help="type yes if you want to force strandedness. If the feature occupies the antisense, strand, the sequence will be reverse complemented.")
    parser.add_argument("--sequence", help="type yes if you want to display corresponding sequence")
    parser.add_argument("--output", help="output file")
    return parser.parse_args()

def main():
    args=parse_arguments()
    (length_consistancy, length, error, test_length_motif, bedName)= correspondence_matrix_interval(args.bed, args.matrix)
    if length_consistancy==False:
        error_str=[]#list that will contain the line numbers in str format
        for number in error:
            number=str(number)#conversion of numbers in int form to str 
            error_str.append(number)
        errors=",".join(error_str) #creation of a string containing the line numbers separated by commas
        print("Error in interval length at line(s) {} of the bed file".format(errors))
        quit()
    if test_length_motif==False:
        print("ERROR: there is a difference between the length of the motif and the length of the intervals in the bed file")
        quit()
    # GC content
    if args.gc == "-1": #the percentage of GC is calculated from the fasta file
        GC_frequency = percent_gc(args.fasta)
        GC_content = float(GC_frequency["G"]) + float(GC_frequency["C"])
    elif float(args.gc) >= 0 and float(args.gc) <= 1 : #the percentage is given by the user if it is between 0 and 1
        GC_frequency = dictionary_gc(args.gc)
        GC_content = args.gc
    else: #default percentage
        GC_frequency = dictionary_gc(0.5)
        #GC_frequency={"A":0.25,"T":0.25,"C":0.25,"G":0.25}
        GC_content = 0.5
    if args.sequence=="yes":
        sequence_option=True
    else:
        sequence_option=False
    matrix_dictionary=matrix_to_dictionary(args.matrix, length)
    if args.strand=="yes":
        strand=True
    else:
        strand=False
    if args.output:
        calcul_score(args.bed,args.fasta,args.output,matrix_dictionary,length,strand,sequence_option, GC_content)
    else:
        calcul_score(args.bed,args.fasta,"score_{}.bed".format(bedName),matrix_dictionary,length,strand,sequence_option, GC_content)

#main

main()
#help('score_PWM')