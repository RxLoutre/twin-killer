# From a fasta file that contains several sequences, remove one/sev sequences based on its/their names.
# Input				- fasta file
#					- a file with a list of ID to remove
# Output			- directly the cleaned fasta file.

import os
import sys

#Where the sequences to blast are put
file = sys.argv[1]
#The list file with ID to remove
toRemove = sys.argv[2]

#Read the file with the ID list
f = open(toRemove,'r')
# list of seq names to remove
lnames=list()
for line in f:
	lnames.append(line.rstrip())
f.close()

try:
	infile=open(file, 'r')
except IOError, e:
	print "Unknown file: ",file
	sys.exit()

OK=True
for line in infile:
	if line[0] == ">": # we enter a new sequence.
		name=line.split()[0].split(">")[1]
		if name in lnames: # is it a bad one that needs to be removed
			OK=False
		else: # We have entered a good sequence.
			OK=True
			print line.strip()
	else:# we are not within a sequence header. 
		if OK: #within an OK sequence. To keep
			print line.strip()
infile.close()
