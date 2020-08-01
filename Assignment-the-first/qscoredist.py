#!/usr/bin/env python

import argparse
import matplotlib
import gzip
from matplotlib import pyplot as plt

import argparse
def get_args():
    parser = argparse.ArgumentParser("a program to produce average quality score histograms for each nucleotide position")
    parser.add_argument("-f", "--file", type=str, help="fastq file to pass through program", required=True)
    parser.add_argument("-l", "--length", type=str, help="length of qscore reads", required=True)
    parser.add_argument("-o", "--output", type=str, help="name of output file", required=True)
    return parser.parse_args()
args = get_args()

###assign argument for file to variable inside of program
file1 = args.file
lenqscore = args.length
output1 = args.output

def init_list(array, value=0.0):
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    init_list = [0.0]*int(lenqscore)
    return init_list

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter)-33

def populate_list(file):
    """This function opens a fastq file and uses the quality score line to
    calculate the sum of the quality scores at each position in the sequence"""
    mean_scores = init_list([])
    line_count = 0
    with gzip.open (file, 'rt') as fh:
        for line in fh:
            phred_counter = 0
            line_count +=1
            if line_count %4 == 0:
                qscore = line.strip()
                for character in qscore:
                    mean_scores[phred_counter] = mean_scores[phred_counter] + convert_phred(character)
                    phred_counter += 1
        return mean_scores, line_count

mean_scores, linect = populate_list(file1)
#print(linect)
#print(mean_scores)

def avg_qualityscores(mean_scores):
    '''This function takes quality score sequences and returns the avg qscore for each position in the sequence'''
    x=[]
    y=[]
    line_count = 0
    #print("# Base Pair\t"+"Mean Quality Score")
    for sum1 in mean_scores:
        mean = sum1 / (linect/4)
        #print(str(line_count)+"\t"+str(mean))
        mean_scores[line_count] = mean
        x.append(line_count)
        y.append(mean)
        line_count+=1
    return x, y

x, y = avg_qualityscores(mean_scores)

###plotting
x = x
y = y
plt.bar(x,y,align='center',alpha=0.75,color='green')
plt.xlabel('Nucleotide Position')
plt.ylabel('Quality Score')
plt.title(str(output1)+'\nMean Quality Scores at each Position in Sequence')
plt.rcParams['figure.figsize'] = (15,11)
plt.savefig('qscore_hist'+str(output1)+'.png')