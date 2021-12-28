#!/usr/bin/python

import subprocess, argparse

# Parses the arguments that are passed in with the executable, which include
# a sequence string to search for, the output prefix used to name result files,
# and the input path for the fastq file. It returns these arguments.
def getArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', metavar='searchstring', help="sequence to search through")
    parser.add_argument('-o', metavar='output_prefix', help="The prefix to give result files")
    parser.add_argument('-i', metavar='input_path', help="The fastq file, including its path")
    args = parser.parse_args()
    return args

# Takes in the output prefix, file path, sequence to search for, and length of the sequence
# and runs the bbduk.sh file to search for matches and save those matches in
# .fastq format in the current directory that this python command is used.
# returns the number of matches
def findMatches(outname, in_path, searchstring, sequencelen):
    subprocess.run(["bbmap/bbduk.sh","in="+in_path, "outm="+outname+"_results.fastq", "literal="+searchstring, "k="+sequencelen, "copyundefined", "mm=f", "rcomp=f"])

# Takes in the newly generated fastq file name from findmatches
# and finds/returns the total number of matches that were found by bbduk.
def findNumMatches(filename):
    with open(filename, 'r') as f:
        numEntries = int(len(f.readlines()))
    return numEntries/4

# Passes in the newly generated fastq filename from findmatches,
# the number of matches of the searchstring, the searchstring, and its length
# and returns a list of the sequnces (of length 4) that follow the search_string
# in each match. For example in the Sequence 'AAACCCTTTGGG', if the searchstring
# was 'CCC', this would return a ['TTTG']
def getAfterBases(filename, nMatches, searchstring, sequencelen):
    sequences=[]
    afterBases=[]
    with open(filename, 'r') as f:
        content = f.read().splitlines()
    for i in range(0,int(nMatches)):
        sequences.append(content[i*4+1])
        buff = sequences[i].find(searchstring)
        afterBases.append("")
        for a in range(0,4):
            afterBases[i]=afterBases[i] + sequences[i][buff+a]
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
            percentages[i][p] = (percentages[i][p] / seqLength) * 100

    return percentages

# Create a _report.txt file that has the total count of the sequence being searched for
# and the base-position percentages
def createReport(outname,in_path,searchstring,numMatches,sequences,percentages):
    with open(outname+"_report.txt", "w+") as f:
        f.write("Results while searching for the search string:\n\n" + searchstring +"\n\n")
        f.write("in:\n\n" + in_path + "\n\n")
        f.write("The string was found " + str(int(numMatches)) + " times.\n\n")
        f.write("For the four nucleotides AFTER this sequence, they contain the following percentages:\n\n")
        f.write("XNNN\nA: " + str(percentages[0][0]) + " C: " + str(percentages[0][1]) + " G: " + str(percentages[0][2]) + " T: " + str(percentages[0][3]) + "\n")
        f.write("NXNN\nA: " + str(percentages[1][0]) + " C: " + str(percentages[1][1]) + " G: " + str(percentages[1][2]) + " T: " + str(percentages[1][3]) + "\n")
        f.write("NNXN\nA: " + str(percentages[2][0]) + " C: " + str(percentages[2][1]) + " G: " + str(percentages[2][2]) + " T: " + str(percentages[2][3]) + "\n")
        f.write("NNNX\nA: " + str(percentages[3][0]) + " C: " + str(percentages[3][1]) + " G: " + str(percentages[3][2]) + " T: " + str(percentages[3][3]) + "\n\n")

def main():
    args = getArguments()
    outname, in_path, searchstring, sequencelen = args.o, args.i, args.s, str(len(args.s))
    findMatches(outname, in_path, searchstring, sequencelen)
    numMatches = findNumMatches(outname+"_results.fastq")
    sequences = getAfterBases(outname+"_results.fastq", numMatches, searchstring, sequencelen)
    percentages = percentBases(sequences, 4)
    createReport(outname,in_path,searchstring,numMatches,sequences,percentages)

if __name__ == "__main__":
    main()
