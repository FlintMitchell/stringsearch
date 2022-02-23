# stringsearch
Searching for and analyzing substrings within fastq files 

Requires Python3, wget, git, bbmap, numpy, and matplotlib (this is installed via the instructions below)

## Installation.

Copy and paste the following instructions:

```
pip3 install matplotlib numpy
git clone https://github.com/FlintMitchell/stringsearch.git
cd stringsearch
wget 'https://sourceforge.net/projects/bbmap/files/BBMap_38.94.tar.gz'
tar -zxvf BBMap_38.94.tar.gz
rm BBMap_38.94.tar.gz
chmod +x stringsearch.py

```

After this, find the path into the stringsearch directory:
```
pwd
```

Then go to your desktop and edit your .bashrc or .bash_profile file by using:
```
vim .bash_profile
```
scroll to the bottom, press `i` to edit, add a new line and then the following:
```
export PATH=$PATH:[PATH/TO/STRINGSEARCH]
```
for example, this might look like:
```
export PATH=$PATH:/Users/flintmitchell/Desktop/stringsearch
```

Now, exit your terminal, reopen it, and the stringsearch.py program should be available
to use anywhere on your computer from the terminal.

## Usage and Output

This command uses three options:
```
-i PATH/TO/INPUT/file.fastq
-o desired_output_prefix
-s string to search for
```
example:
```
python3 stringsearch.py -i mydata/file001.fastq -o file001 -s ACGTACGTATGACTACGTCAG
```

This will output two files:
```
file001_report.txt      Contains the total number of matches of your string in the fastq file and the base-composition of the 4 bases following those matches
file001_result.fastq    Contains the fastq entries from your input fastq file that had matches 
```
