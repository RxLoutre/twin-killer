"""
Differents class allowing to store data coming from a BUSCO output
"""

import csv
import argparse

class Contig:
	"""
	A contig is defined only by it's ID, it's length
	a gstart and gend, representing the genomic start/end on the whole genome
	and astart and aend representing the alignement start and end 
	"""
	def __init__(self,contig_id,gstart,gend,length,astart,aend):
		self.contig_id = contig_id
		self.gstart = gstart
		self.gend = gend
		self.length = length
		self.astart = astart
		self.aend = aend
		
	def __str__(self):
		s = self.contig_id + " " + str(self.length) + " " + str(self.astart) + " " + str(self.aend)
		return s
		
def make_contig(contig_id,start,end,length,astart,aend):
	"""
	A function that initialize a contig
	"""
	contig = Contig(contig_id,start,end,length,astart,aend)
	return contig

def read_contigs(index_file):
	"""
	Read a .fai file in order to create a primary contig_dictionary that provide informations
	about the contigs, like it's genomic start and end, and also it's length
	the astart and a end are set to 0 because it's not the final dictionnary
	The astart and aend information is onmy used when creating a list of contigs within a
	duplication dictionnary (see below)
	"""
	contig_dic={}
	with open(index_file,'r') as f:
		reader = csv.reader(f,delimiter='\t')
		for contig_id,length,start,fa1,fa2 in reader:
			c = make_contig(contig_id,int(start),int(start)+int(length)-1,int(length),0,0)
			contig_dic[contig_id]=c
			
	return contig_dic

	
def read_duplication_file(full_table_file,contig_dic):
	"""
	Read a output result of BUSCO, (full_table), but only with the DUPLICATED lines
	This function use some informations contained if the contig dictionnary, created with
	the .fai index files, and contained in the busco output file in order de create a dictionnary.
	This dictionnary has Busco ID as key, and as element, a list of contigs.
	The contigs are created one by one, partially copied with the correspondant contig from the
	contig_dictionnary, and with the astart and aend informations added.
	"""
	duplication_dic = {}
	with open(full_table_file,'r') as f:
		reader = csv.reader(f,delimiter='\t')
		for busco_id,tag,contig_id,start,end,score1,score2 in reader:
			contig_gstart = contig_dic[contig_id].gstart
			contig_gend = contig_dic[contig_id].gend
			contig_length = contig_dic[contig_id].length
			c = make_contig(contig_id,contig_gstart,contig_gend,contig_length,start,end[0:-1])
			if busco_id in duplication_dic:				
				duplication_dic[busco_id].append(c)
			else:
				duplication_dic[busco_id]=[]
				duplication_dic[busco_id].append(c)
			
	return duplication_dic
	
def print_dic_duplicated(duplication_dic):
	"""
	Does exactly what it says
	"""
	for key in duplication_dic:
		print key
		for elem in duplication_dic[key]:
			print elem
			
def print_elem_dic_duplicated(duplication_dic,key):
	"""
	Does exactly what it says
	"""
	for elem in duplication_dic[key]:
		print elem
		
def duplication_analysis(duplication_dic):
	ratios = []
	s = ""
	for key in duplication_dic:
		max_len = 0
		max_contig = ""
		for contig in duplication_dic[key]:
			if(contig.length > max_len):
				max_contig_index = str(duplication_dic[key].index(contig))
				max_len = contig.length
		max_contig_id = duplication_dic[key][int(max_contig_index)].contig_id
		#print "Le plus grand contig est : " + duplication_dic[key][int(max_contig_index)].contig_id
		for contig in duplication_dic[key]:
			if contig.contig_id != max_contig_id:
				#ratios.append(float(contig.length) / float(max_len))
				s +=  max_contig_id + " " + contig.contig_id + "\n"
	return s
	
parser = argparse.ArgumentParser()
parser.add_argument("full_table", help="The full_table Busco generated")
parser.add_argument("index", help="The index file of your assembly (must be created with IGV)")
parser.add_argument("output", help="The contig couple file. (")
args = parser.parse_args()
duplicated_file = args.full_table
index_file = args.index
contig_dico = read_contigs(index_file)
duplication_dico = read_duplication_file(duplicated_file,contig_dico)
#print_dic_duplicated(duplication_dico)
string = duplication_analysis(duplication_dico)
output = open(args.output,'w')
output.write(str(string))
output.close()

			
		


	
