#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: laurinerolland

"""
import pybedtools
import argparse
import os
from os.path import basename
import subprocess

def length_interval(bedFile): #test the length of the intervals
    ''' test the length of the intervals in the bed file and 
        return a booleen which indicates if the intervals are equal or not
        the length of the motif and the list of lines where the length of the intervals is different, 
        if there is no difference an empty list is returned'''
    try:
        bed=open(bedFile,"r") #the bed file is opened for reading
    except:
        print ("The file {} cannot be found".format(bed))
    else:
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
                else:
                    length_motif=int(end)-int(start)
            line=bed.readline().strip()
            i+=1
        if list_error != []:
            test_length=False
        else:
            test_length=True
        return (test_length, length_motif, list_error) #returns the length of the motif
    
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
    
def creation_PWM(bedFile,fastaFile, length, strand): #creation of the dictionary containing the PWM
        ''' return the name of the bed and fasta files used 
            and the values of the PWM in a dictionary format
            {'A': {postion_1: ..., position_n : ...}, 'T': {postion_1: ..., position_n : ...}, 'C': {postion_1: ..., position_n : ...}, 'G': {postion_1: ..., position_n : ...}}'''
        #recover file names
        bed_name_ext = basename(bedFile) #to get the file name
        bed_name = os.path.splitext(bed_name_ext)[0]
        fasta_name_ext = basename(fastaFile) #to get the file name
        fasta_name = os.path.splitext(fasta_name_ext)[0]
        ##study of the nucleotide sequence composition
        dictionary_nucleotides = {"A":{}, "T":{}, "C":{}, "G":{}} #this dictionary will contain the number of times each nucleotide has been found at each position and then the matrix
        # exemple : {'A': {1: 73, 2: 31, 3: 104, 4: 168, 5: 170, 6: 176, 7: 170, 8: 137, 9: 107, 10: 93}, 'T': {1: 106, 2: 139, 3: 72, 4: 7, 5: 3, 6: 6, 7: 9, 8: 45, 9: 74, 10: 89}, 'C': {1: 7, 2: 5, 3: 9, 4: 9, 5: 12, 6: 3, 7: 6, 8: 4, 9: 3, 10: 6}, 'G': {1: 4, 2: 15, 3: 5, 4: 6, 5: 5, 6: 5, 7: 5, 8: 4, 9: 6, 10: 2}}
        nucleotides_list=list(dictionary_nucleotides.keys())
        j=1 #the counter is initialized to 1 in order to make the first position 1
        while j<(length+1): #for each position
            for nucleotide in nucleotides_list: # for each nucleotide of the list
                dictionary_nucleotides[nucleotide][j]=0 #the dictionary is initialized to 0
            j+=1
        #From the bed file containing the pattern positions, find the nucleotide sequences with bedtools getfasta
        a = pybedtools.BedTool(bedFile)
        if strand==True:
            a = a.sequence(fi=fastaFile,s=True)#mettre en option
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
        line=sequence_fasta.readline() #the file is read line by line
        counter_line=0 #number of sequence
        while line != "": #We browse all the lines of the file
            if line[0] != ">": #lines beginning with > do not contain the sequence
                nucleotide_sequence_line=line.upper()#the sequence is uppercased because it is not uppercased on all lines
                i=0 #this counter is used to browse the sequence, it is initialized to 0
                while i<(length):
                    for nucleotide in nucleotides_list: # for each nucleotide of the list
                        if nucleotide_sequence_line[i]==nucleotide:#we look at which nucleotide corresponds each position of the sequence
                            dictionary_nucleotides[nucleotide][i+1]=dictionary_nucleotides[nucleotide][i+1] +1 # we add 1 to the counter i because the positions in the dictionary start from 1
                    i+=1
                    counter_line+=1
            line=sequence_fasta.readline()
        return (bed_name, fasta_name, dictionary_nucleotides)     
    
def position_to_matrix(dictionary_nucleotides):
    dictionary_matrix={"A":{}, "T":{}, "C":{}, "G":{}} 
    nucleotides_list=list(dictionary_nucleotides.keys())
    #frequency of nucleotides at each position
    positions=list(dictionary_nucleotides["A"].keys())
    for nucleotide in nucleotides_list:
        l=1
        while l<(len(positions)+1):
            dictionary_matrix[nucleotide][l]= (dictionary_nucleotides[nucleotide][l])/(dictionary_nucleotides["A"][l] + dictionary_nucleotides["T"][l] + dictionary_nucleotides["C"][l] + dictionary_nucleotides["G"][l])
            #dictionary_matrix[nucleotide][l]= (dictionary_nucleotides[nucleotide][l])/(counter_line/length)
            #if dictionary_matrix[nucleotide][l] !=0:
                #dictionary_matrix[nucleotide][l]= math.log2((dictionary_matrix[nucleotide][l])/(gc_freq[nucleotide]))
            #else:
                #dictionary_matrix[nucleotide][l]=0
            l+=1
        #W(b,i)=log2(p(b,i)/p(b))
        #p(b) = expected background frequency of letter b
        #p(b, i) is the frequency of letter b at position i
    return (dictionary_matrix)

def writePWM_homer(fichier,dictionary_matrix,length,bed_name, gc_content): #saves in a file a matrix in the format requested in the arguments (add parameters giving information about the titles --> they can be added in the header)
    '''create a file where the matrix is displayed in HOMER format'''
    f=open(fichier,"w")
    f.write("> motif L1 \t threshold \t {} \t expected gc : {}\n".format(bed_name,gc_content))
    i=1
    while i<(length+1):
        f.write("{:15} \t {:15} \t {:15} \t {:15} \n".format(dictionary_matrix["A"][i],dictionary_matrix["C"][i],dictionary_matrix["G"][i], dictionary_matrix["T"][i]))
        i+=1 
    f.close()
    print("The positon weight matrix has been generated")

def writePWM_jaspar(fichier,dictionary_matrix,length,bed_name, gc_content): #saves in a file a matrix in the format requested in the arguments (add parameters giving information about the titles --> they can be added in the header)
    '''create a file where the matrix is displayed in JASPAR format'''
    f=open(fichier,"w")
    f.write("> motif L1 \t threshold \t {} \t expected gc : {}\n".format(bed_name,gc_content))
    nucleotide_order=["A","C","G","T"]
    for nucleotide in nucleotide_order:
        f.write("{} [".format(nucleotide))
        i=1
        while i<length+1:
            f.write("{} ".format(dictionary_matrix[nucleotide][i]))
            i+=1 
        f.write("]\n")
    f.close()
    print("The positon weight matrix has been generated")
    
def writePWM_transfac(fichier, dictionary_matrix, length, matrix_ID, matrix_name, matrix_class): #saves in a file a matrix in the format requested in the arguments (add parameters giving information about the titles --> they can be added in the header)
    '''create a file where the matrix is displayed in TRANSFAC format'''
    f=open(fichier,"w")
    f.write("ID \t any_old_name_for_motif_1 \n")
    f.write("BF \t species_name_for_motif_1 \n")
    #f.write("DE \t {} \t {} \t {}\n".format(matrix_ID,matrix_name, matrix_class))
    f.write("P0 \t A \t C \t G \t T \n")
    i=1
    while i<(length+1):
        if i<11:
            f.write("0{} \t {} \t {} \t {} \t {} \n".format(i-1,dictionary_matrix["A"][i],dictionary_matrix["C"][i],dictionary_matrix["G"][i], dictionary_matrix["T"][i]))
        else:
            f.write("{} \t {} \t {} \t {} \t {} \n".format(i-1,dictionary_matrix["A"][i],dictionary_matrix["C"][i],dictionary_matrix["G"][i], dictionary_matrix["T"][i]))
        i+=1 
    f.write("XX \n")
    f.write("//")
    f.close()
    print("The positon weight matrix has been generated")
    
def writePWM_meme(fichier,dictionary_matrix,length, gc_freq): #saves in a file a matrix in the format requested in the arguments (add parameters giving information about the titles --> they can be added in the header)
    '''create a file where the matrix is displayed in MEME format'''
    f=open(fichier,"w")
    f.write("MEME version 4\n")
    f.write("\n")
    f.write("ALPHABET = ACGT\n")
    f.write("\n")
    f.write("strands: + - \n")
    f.write("\n")
    f.write("Background letter frequencies \n")
    f.write("\n")
    f.write("A {} C {} G {} T {}\n".format(gc_freq["A"], gc_freq["C"], gc_freq["G"], gc_freq["T"]))
    f.write("\n")
    f.write("MOTIF [identifier] L1[alternate name] \n")
    f.write("letter-probability matrix: alength= 4 w= {} nsites= [source site] E= [source E-value] \n".format(length))
    i=1
    while i<(length+1):
        f.write("{:15} \t {:15} \t {:15} \t {:15} \n".format(dictionary_matrix["A"][i],dictionary_matrix["C"][i],dictionary_matrix["G"][i], dictionary_matrix["T"][i]))
        i+=1 
    f.write("URL http://jaspar.genereg.net/matrix/MA0048.2")
    f.close()
    print("The positon weight matrix has been generated")

def writePWM_raw(fichier,dictionary_matrix,length,bed_name, gc_content): #saves in a file a matrix in the format requested in the arguments (add parameters giving information about the titles --> they can be added in the header)
    '''create a file where the matrix is displayed in RAW format'''
    f=open(fichier,"w")
    f.write("> motif L1 \t threshold \t {} \t expected gc : {}\n".format(bed_name,gc_content))
    nucleotide_order=["A","C","G","T"]
    for nucleotide in nucleotide_order:
        i=1
        while i<length+1:
            f.write("{} ".format(dictionary_matrix[nucleotide][i]))
            i+=1 #ordre A,C,G,T 
        f.write("\n")
    f.close()
    print("The positon weight matrix has been generated")
        
def parse_arguments():
    '''creation of options 
    add -h to the command line to display the different options'''
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", help="bed file that will be used")
    parser.add_argument("--fasta", help="fasta file that will be used")
    parser.add_argument("--gc", help="percentage of GC which will be displayed in the title of the generated matrix")
    parser.add_argument("--output", help="output file")
    parser.add_argument("--strand", help="type yes if you want to force strandedness. If the feature occupies the antisense, strand, the sequence will be reverse complemented.")
    parser.add_argument("--PWM", help="type yes if you want to generate the PWM")
    parser.add_argument("--format", choices=['HOMER', 'JASPAR', 'TRANSFAC', "MEME", "RAW"], help="format of the output file")
    parser.add_argument("--weblogo", help="type the name to give to the weblogo file. Format :eps, png, png_print, pdf, jpeg, svg, logodata")
    #parser.print_help() #to display the help
    return parser.parse_args()

def main():
    args=parse_arguments()
    # test if the bed file is entered
    if args.bed == None:
        print("ERROR: The bed file was not entered")
        quit()
    # test if the fasta file is entered
    if args.fasta == None:
        print("ERROR: The fasta file was not entered")
        quit()
    # test intervals
    (length_consistancy,length_motif, error)=length_interval(args.bed) #if the intervals are all equal we return true and the length of the motif
    if length_consistancy==False: #if the intervals are not all equal (booleen==False)
        error_str=[]#list that will contain the line numbers in str format
        for number in error:
            number=str(number)#conversion of numbers in int form to str 
            error_str.append(number)
        errors=",".join(error_str) #creation of a string containing the line numbers separated by commas
        print("Error in interval length at line(s) {} of the bed file".format(errors))
        quit()
    # creation of the PWM in a dictionary format
    if args.strand=="yes":
        strand=True
    else:
        strand=False
    (name_bedFile, name_fastaFile, counter_position)=creation_PWM(args.bed, args.fasta,length_motif,strand)
    # output the matrix into a file in the required format 
    if args.weblogo:
        writePWM_transfac("transfac.pwm",counter_position,length_motif, "matrix_ID", "matrix_name", "matrix_class")
        cmd = "weblogo -c classic -s large -i -3 < transfac.pwm > {}".format(args.weblogo) #cd ~/Desktop/STAGE/Workspace/Python |
        subprocess.call(cmd, shell=True)
        print("The WebLogo has been generated")
    if args.PWM=="yes":
        dict_matrix=position_to_matrix(counter_position)
        #print(dict_matrix)
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
        # Format matrix 
        if args.format=="HOMER":
            #name of the output file 
            if args.output: # if the name of the output file is given by the user 
                writePWM_homer(args.output,dict_matrix,length_motif, name_bedFile, GC_content)
            else:
                writePWM_homer("MotifL1_PWM_homer.motif", dict_matrix, length_motif, name_bedFile, GC_content) 
        elif args.format=="JASPAR":
            #name of the output file 
            if args.output: # if the name of the output file is given by the user 
                writePWM_jaspar(args.output,dict_matrix,length_motif, name_bedFile, GC_content)
            else:
                writePWM_jaspar("MotifL1_PWM_jaspar.motif", dict_matrix, length_motif, name_bedFile, GC_content) 
        elif args.format=="TRANSFAC":
            #name of the output file 
            if args.output: # if the name of the output file is given by the user 
                writePWM_transfac(args.output,dict_matrix,length_motif, "matrix_ID", "matrix_name", "matrix_class")
            else:
                writePWM_transfac("MotifL1_PWM_transfac.motif", dict_matrix, length_motif, "matrix_ID", "matrix_name", "matrix_class")
        elif args.format=="MEME":
            #name of the output file 
            if args.output: # if the name of the output file is given by the user 
                writePWM_meme(args.output,dict_matrix,length_motif, GC_frequency)
            else:
                writePWM_meme("MotifL1_PWM_meme.motif", dict_matrix, length_motif, GC_frequency)
        elif args.format=="RAW":
            #name of the output file 
            if args.output: # if the name of the output file is given by the user 
                writePWM_raw(args.output,dict_matrix,length_motif, name_bedFile, GC_content)
            else:
                writePWM_raw("MotifL1_PWM_raw.motif", dict_matrix, length_motif, name_bedFile, GC_content)
        else: # by default homer format
            if args.output: # if the name of the output file is given by the user 
                writePWM_homer(args.output,dict_matrix,length_motif, name_bedFile, GC_content)
            else:
                writePWM_homer("MotifL1_PWM_homer.motif", dict_matrix, length_motif, name_bedFile, GC_content)

    
#Main program

main()
#help('creation_PWM')