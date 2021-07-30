#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 16:04:06 2021

@author: laurinerolland
"""

import argparse
import re

expr = r'^>(?P<number>\w+)'
chrLine = re.compile(expr)

def positionN(fasta):
    try:
        file=open(fasta,"r") #the annotated file is opened for reading
    except:
        print ("The file {} cannot be found".format(fasta))
    else:
        line=file.readline().strip() #read the line and delete the \n
        dico = {}
        count=1
        past_letter = "" #initialisation
        chrNumber = "" # initialisation
        while line!="": # we browse the file as long as it is not empty
            chr_number = chrLine.search(line) #search if the line matches the regular expression
            if chr_number: #if the line matches the regular expression
                if past_letter =="N": #if the sequence of the previously read chromosome ends with N 
                    dico[chrNumber].append(count) # the end position is added to the dictionary
                chrNumber = chr_number.group("number") #chr1,chr2...
                dico[chrNumber]=[] #{chr1:[]}
                count=1 # initialisation
                past_letter = "" 
            else:
                for letter in line: #we browse the letters in the line
                    if letter == "N":
                        if past_letter != "N": #if the previous letter wasn't N
                            dico[chrNumber].append(count)
                    else:
                        if past_letter == "N": #  #if the previous letter was N
                            dico[chrNumber].append(count-1)
                    count+=1
                    past_letter = letter
            line=file.readline().strip()
        print(dico)
        return (dico)

def bedfile(dico, file):
    f=open(file,"w")
    for chromosome in dico:
        position = dico[chromosome]
        if position != [] and position != [0]:
            i=0
            while i < (len(position)-1):
                f.write("{}\t{}\t{}".format(chromosome,position[i],position[i+1]))
                f.write("\n")  
                i+=2
        
def parse_arguments():
    '''creation of options 
    add -h to the command line to display the different options'''
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", help="fasta file that will be used")
    parser.add_argument("--output", help="output file")
    return parser.parse_args()                    

def main(): 
    args=parse_arguments()
    if args.fasta:
        fasta_file = args.fasta
        dictionnary = positionN(fasta_file)
        output_file = args.output
        bedfile(dictionnary,output_file)
        print("The file have been created")
        
main()
            
            