#!/usr/bin/env python

#######################################################
##    a program to demultiplex and quality screen    ##
##            illumina fastq output                  ##
####################################################### 

###MODULES
import fileinput
import re
import gzip
import argparse
import itertools

###ARGPARSE

def get_args():
    parser = argparse.ArgumentParser("a program to demultiplex and quality screen illumina fastq output")
    parser.add_argument("-R1", "--fastqR1", type=str, help="R1 fastq illumina output(biological read 1)", required=True)
    parser.add_argument("-R2", "--fastqR2", type=str, help="R2 fastq illumina output (index read 1)", required=True)
    parser.add_argument("-R3", "--fastqR3", type=str, help="R3 fastq illumina output (index read 2)", required=True)
    parser.add_argument("-R4", "--fastqR4", type=str, help="R4 fastq illumina output (biological read 2)", required=True)
    parser.add_argument("-b", "--barcodes", type=str, help="list of barcodes/indices", required=True)
    parser.add_argument("-q", "--cutoff", type=str, help="quality score cutoff", required=True)
    parser.add_argument("-o", "--output", type=str, help="naming prefix", required=True)
    return parser.parse_args()

args = get_args()

read1 = args.fastqR1
index1 = args.fastqR2
index2 = args.fastqR3
read2 = args.fastqR4
barcodes = args.barcodes
cutoff = args.cutoff
output = args.output


#########################FUNCTIONS###############################

def get_barcode_qscore(index_qscore):
    '''this function converts the letter/symbol in the qscore line to a qscore integer, 
    then calculates the average qscore for the string'''
    sumqscores = 0
    for letter in index_qscore:
        qscore=ord(letter)-33
        sumqscores += qscore
    avg_qscore = sumqscores/len(index_qscore)
    return avg_qscore 

def get_qscore(letter):
    '''This function will return the qscore for a given letter'''
    qscore=ord(str(letter))-33
    return qscore

#print(get_qscore('I'))

line_counter = 1
barcode_list = []
meta_data = {}
with open(barcodes, 'r') as fhbc:
    for line in fhbc:
        if line_counter > 1:
            barcode_split = line.strip().split("\t")
            barcode_seq = barcode_split[4]
            barcode_list.append(barcode_seq)
            smple = barcode_split[0]
            sample = smple.strip()
            grp = barcode_split[1]
            group = grp.strip()
            trtmnt = barcode_split[2]
            treatment = trtmnt.strip()
            ind_name = barcode_split[3]
            index_name = ind_name.strip()
            meta_data[barcode_seq] = [sample, group, treatment, index_name]
        line_counter += 1
permutations=itertools.product(barcode_list, repeat=2)

index_counter = {}
for perm in permutations:
    index_counter[perm] = 0


def get_reverse_complement(index2):
    '''this function will produce the reverse complement of a DNA strand'''
    output = ''
    for letter in index2:
        letter = letter.upper()
        if letter == 'A':
            output += 'T'
        elif letter == 'T':
            output += 'A'
        elif letter == 'G':
            output += 'C'
        elif letter == 'C':
            output += 'G'
        else:
            output += 'N'
    return(output[::-1])

def add_barcode_to_header(header,barcode1,barcode2):
    '''this function will add the barcode pair to the header line of a record'''
    str(header)
    edited_header = (str(header)+'_'+str(barcode1)+'-'+str(barcode2))
    return edited_header

#opening files for matched indices/barcodes and initialize counters in a loop
filenames_dict={}
filenames_list=[]
for key, value in meta_data.items():
    read1file = open(output+'_'+str(key)+'_'+str(value[3])+'_'+str(value[2])+'_'+str(value[1])+'_rd1.fastq', "w")
    read2file = open(output+'_'+str(key)+'_'+str(value[3])+'_'+str(value[2])+'_'+str(value[1])+'_rd2.fastq', "w")
    filenames_dict[key]=[read1file,read2file]
    filenames_list.append(read1file)
    filenames_list.append(read2file)

###initialize remaining counters
poorqual = 0
matched_unknown_barcodes = 0
index_hopped = 0
total_count = 0
total_matched = 0

############# OPENING THE REMAINING OUTPUT FILES ####################

index_hopped_read1 = open(output+"_hopped_rd1.fastq", "w")
index_hopped_read2 = open(output+"_hopped_rd2.fastq", "w")
poorqual_read1 = open(output+"_poorqual_rd1.fastq", "w")
poorqual_read2 = open(output+"_poorqual_rd2.fastq", "w")
summary_table = open(output+"_demultiplexing_summary", "w")

#################################### DEMULTIPLEXING ####################################
########################################################################################


with gzip.open(read1, 'rt') as rd1, gzip.open(index1, 'rt') as in1, gzip.open(index2, 'rt') as in2, gzip.open(read2, 'rt') as rd2:
    while True:
        header_rd1 = rd1.readline().rstrip()
        if header_rd1 == '':
            break
        seq_rd1 = rd1.readline().rstrip()
        comment_rd1 = rd1.readline().rstrip()
        qscore_rd1 = rd1.readline().rstrip()

        header_in1 = in1.readline().rstrip()
        seq_in1 = in1.readline().rstrip()
        comment_in1 = in1.readline().rstrip()
        qscore_in1 = in1.readline().rstrip()

        header_in2 = in2.readline().rstrip()
        seq_in2 = get_reverse_complement(in2.readline().rstrip())
        comment_ind2 = in2.readline().rstrip()
        qscore_in2 = in2.readline().rstrip()

        header_rd2 = rd2.readline().rstrip()
        seq_rd2 = rd2.readline().rstrip()
        comment_rd2 = rd2.readline().rstrip()
        qscore_rd2 = rd2.readline().rstrip()
        
        lowq = False
        for letter in qscore_in1 + qscore_in2:
            if get_qscore(letter) < int(cutoff):
                lowq = True

        if lowq:
            poorqual_read1.write(str(add_barcode_to_header(header_rd1,seq_in1,seq_in2))+'\n'+str(seq_rd1)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
            poorqual_read2.write(str(add_barcode_to_header(header_rd2,seq_in1,seq_in2))+'\n'+str(seq_rd2)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
            poorqual += 1
            total_count += 1

        elif seq_in1 != seq_in2:
            if seq_in1 not in barcode_list:
                poorqual_read1.write(str(add_barcode_to_header(header_rd1,seq_in1,seq_in2))+'\n'+str(seq_rd1)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                poorqual_read2.write(str(add_barcode_to_header(header_rd2,seq_in1,seq_in2))+'\n'+str(seq_rd2)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                poorqual += 1
                total_count += 1

            elif seq_in2 not in barcode_list:
                poorqual_read1.write(str(add_barcode_to_header(header_rd1,seq_in1,seq_in2))+'\n'+str(seq_rd1)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                poorqual_read2.write(str(add_barcode_to_header(header_rd2,seq_in1,seq_in2))+'\n'+str(seq_rd2)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                poorqual += 1
                total_count += 1

            else:
                index_hopped_read1.write(str(add_barcode_to_header(header_rd1,seq_in1,seq_in2))+'\n'+str(seq_rd1)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                index_hopped_read2.write(str(add_barcode_to_header(header_rd2,seq_in1,seq_in2))+'\n'+str(seq_rd2)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                index_hopped += 1
                index_counter[seq_in1,seq_in2] += 1
                total_count += 1

        elif seq_in1 == seq_in2:
            if seq_in1 not in barcode_list:
                poorqual_read1.write(str(add_barcode_to_header(header_rd1,seq_in1,seq_in2))+'\n'+str(seq_rd1)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                poorqual_read2.write(str(add_barcode_to_header(header_rd2,seq_in1,seq_in2))+'\n'+str(seq_rd2)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                matched_unknown_barcodes += 1
                total_count += 1
            elif seq_in2 not in barcode_list:
                poorqual_read1.write(str(add_barcode_to_header(header_rd1,seq_in1,seq_in2))+'\n'+str(seq_rd1)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                poorqual_read2.write(str(add_barcode_to_header(header_rd2,seq_in1,seq_in2))+'\n'+str(seq_rd2)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                matched_unknown_barcodes += 1
                total_count += 1
            else:
                filenames_dict[seq_in1][0].write(str(add_barcode_to_header(header_rd1,seq_in1,seq_in2))+'\n'+str(seq_rd1)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                filenames_dict[seq_in1][1].write(str(add_barcode_to_header(header_rd2,seq_in1,seq_in2))+'\n'+str(seq_rd2)+'\n'+str(comment_in1)+'\n'+str(qscore_rd1)+'\n')
                index_counter[seq_in1,seq_in2] += 1
                total_count += 1
                total_matched += 1



summary_table.write('category\tnumber of records\tpercent of total records\n')
summary_table.write('total matched\t'+str(total_matched)+'\t'+str((total_matched/total_count)*100)+'%\n')
summary_table.write('total index_hopped\t'+str(index_hopped)+'\t'+str((index_hopped/total_count)*100)+'%\n')
summary_table.write('matched unknown indices\t'+str(matched_unknown_barcodes)+'\t'+str((matched_unknown_barcodes/total_count)*100)+'%\n')
summary_table.write('poor quality: qscore cutoff of'+str(cutoff)+'\t'+str(poorqual)+'\t'+str((poorqual/total_count)*100)+'%\n')
summary_table.write('permutations summary\n')
for key, value in index_counter.items():
    summary_table.write('un/matched index:'+str(key)+'\t'+str(value)+'\t'+str((value/total_count)*100)+'%\n')


####### CLOSING MATCHED FILES ######
for fastqfile in filenames_list:
    fastqfile.close()
    #print(file.closed)
####################################

############# CLOSING REMAINING OUTPUT FILES ################
index_hopped_read1.close()
index_hopped_read2.close()
poorqual_read1.close()
poorqual_read2.close()
summary_table.close()
#############################################################