from blast import parse
import csv
import argparse

EVALUE_THREESHOLD = 0.0001
IDENTITY_THREESHOLD = 80

parser = argparse.ArgumentParser()
parser.add_argument("input", help="The input blast result in outfmt 6")
parser.add_argument("output", help="The filtered blast file")
args = parser.parse_args()

fh = open(args.input,'r')
fout = open(args.output,'w')
good_hsp = 0
all_hsp = 0
print("Filtering good HSPs...")
for blast_record in parse(fh):
    #print('query id: {}'.format(blast_record.qid))
    for hit in blast_record.hits:
        for hsp in hit:
            #print('****Alignment****')
            #print('sequence:', hsp.sid)
            #print('length:', hsp.length)
            #print('e value:', hsp.evalue)
            #If we are under the evalue threeshold, keep the alignment
            all_hsp += 1
            if(hsp.evalue < EVALUE_THREESHOLD and hsp.pident > IDENTITY_THREESHOLD):
				good_hsp += 1
				fout.write(str(hsp)+"\n")
				#print('****Alignment****')
				#print('sequence:', hsp.sid)
				#print('length:', hsp.length)
				#print('e value:', hsp.evalue)
				
print("Analyzed HSP : "+str(all_hsp))
print("Retained HSP : "+str(good_hsp))
fh.close()
