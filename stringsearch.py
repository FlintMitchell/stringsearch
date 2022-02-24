#!/usr/bin/env python3

import subprocess, argparse, os, gzip
import matplotlib.pyplot as plt
import numpy as np

# Parses the arguments that are passed in with the executable, which include
# a sequence string to search for, the output prefix used to name result files,
# and the input path for the fastq file. It returns these arguments.
def getArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', metavar='searchstring', help="sequence to search through.")
    parser.add_argument('-o', metavar='output_prefix', help="The prefix to give result files.")
    parser.add_argument('-i', metavar='input_path', help="The fastq file, including its path.")
    parser.add_argument('-n', metavar='num_bases_after_string', help="The number of bases to analyze in the entries with a string match.")
    args = parser.parse_args()
    return args

# Takes in the output prefix, file path, sequence to search for, and length of the sequence
# and runs the bbduk.sh file to search for matches and save those matches in
# .fastq format in the current directory that this python command is used.
# returns the number of matches
def findMatches(outname, in_path, searchstring, sequencelen,dependencyPATH):
    subprocess.run([dependencyPATH+"/bbmap/bbduk.sh","in="+in_path, "outm="+outname+"_results.fastq", "literal="+searchstring, "k="+sequencelen, "copyundefined", "mm=f", "rcomp=f"])
    # subprocess.run([dependencyPATH+"/bbmap/bbduk.sh","in="+in_path, "outm="+outname+"_NOMATCH_results.fastq", "literal="+searchstring, "k="+sequencelen, "copyundefined", "rcomp=f"])

# Takes in a fastq file and finds/returns the total number of entries.
def findNumMatches(filename):
    print(filename)
    with open(filename, 'r') as f:
        numEntries = int(len(f.readlines()))
    return numEntries/4

# Passes in the newly generated fastq filename from findmatches,
# the number of matches of the searchstring, the searchstring, and its length
# and returns a list of the sequnces (of length 4) that follow the search_string
# in each match. For example in the Sequence 'AAACCCTTTGGG', if the searchstring
# was 'CCC', this would return a ['TTTG']

# edits: 1) Instead of returning a list, create a fasta file with the sequences
# 2) make it flexible for however many bases afterwards, not just 4
def getAfterBases(filename, nMatches, searchstring, sequencelen, nBasesAfterSeq):
    sequences=[]
    afterBases=[]
    with open(filename, 'r') as f:
        content = f.read().splitlines() #
    for i in range(0,int(nMatches)):
        sequences.append(content[i*4+1]) # get only the sequences from the fastq file and put into a list of sequences
        buff = sequences[i].find(searchstring) # find the location where the search string starts in the current sequence
        afterBases.append("")
        for a in range(0,nBasesAfterSeq):
            afterBases[i]=afterBases[i] + sequences[i][buff+int(sequencelen)+a]
    return afterBases

# Using the list of n-base sequences from getAfterBases(), for each base,
# calculate the percentage of A, C, G, or T at that position.
def percentBases(seqs, numNucs):
    # create list of all the nucleotides in each position
    position = []
    percentages = []
    for nucPos in range(0,numNucs):
        position.append("")
        for seq in seqs:
             position[nucPos] = position[nucPos] + seq[nucPos]

    # within each list, positions vs nucleotides will be as follows:
    # 0 - A; 1 - C; 2 - G; 3 - # T
    # create a list of lists, where each sublist pertains to a sequence and
    # contains the count of each of the 4 nucleotides in that sequence
    for nucPos in position:
        seqLength = len(nucPos)
        buff = [0,0,0,0]
        for i in range(0,seqLength):
            if nucPos[i] == "A":
                buff[0]=buff[0]+1
            elif nucPos[i] == "C":
                buff[1]=buff[1]+1
            elif nucPos[i] == "G":
                buff[2]=buff[2]+1
            elif nucPos[i] == "T":
                buff[3]=buff[3]+1
        percentages.append(buff)

    # Turn counts into percentages
    for i in range(0,len(percentages)):
        for p in range(0,len(percentages[i])):
            percentages[i][p] = round((percentages[i][p] / seqLength) * 100, 2)
    return percentages

# Create a _report.txt file that has the total count of the sequence being searched for
# and the base-position percentages
def createReport(outname,in_path,searchstring,numMatches,sequences,percentages,numBasesAfterString):
    totalInputEntries = int(findNumMatches(in_path))
    with open(outname+"_report.txt", "w+") as f:
        f.write("Results while searching for the search string:\n\n" + str(searchstring) +"\n\n")
        f.write("in:\n\n" + in_path + "\n\n")
        f.write("Total number of entries in the input fastq file: " + str(totalInputEntries) + ".\n")
        f.write("The string was found " + str(int(numMatches)) + " times.\n")
        f.write("Percentage of entries in the input fastq file with the string being searched for: " + str(round((numMatches/totalInputEntries*100),2)) + "%.\n\n")
        f.write("For the " + str(numBasesAfterString) +" nucleotides AFTER this sequence, they contain the following percentages:\n\n")
        f.write("        [  A  ,  C  ,  G  ,  T  ]\n")
        for i in range(0,len(percentages)):
            f.write("base " + str(i+1) + ": " + str(percentages[i]) + "\n")
    return

def main():
    args = getArguments()
    outname, in_path, searchstring, sequencelen, numBasesAfterString = args.o, args.i, args.s, str(len(args.s)), int(args.n)
    dependencyPATH=os.path.abspath(os.path.dirname(__file__))
    findMatches(outname, in_path, searchstring, sequencelen, dependencyPATH)
    numMatches = findNumMatches(outname+"_results.fastq") # Finds/return the total number of matches that were found by bbduk.
    sequences = getAfterBases(outname+"_results.fastq", numMatches, searchstring, sequencelen, numBasesAfterString)
    percentages = percentBases(sequences, numBasesAfterString)
    createReport(outname,in_path,searchstring,numMatches,sequences,percentages,numBasesAfterString)

    # print out report
    with open(outname+"_report.txt", 'r') as f:
        print(f.read())

if __name__ == "__main__":
    main()
