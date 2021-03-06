#############################################################
#### Motif searching tool for FASTA sequences   #############
####    by Elena P Sorokin   December 16, 2014  #############
#############################################################

'''Search for user-defined motif within a FASTA file.
Return # motif occurrences and motif location(s).'''



# 1. IMPORT MODULES, MEME RESULT FILE, AND FASTA TO SEARCH

# modules

from Bio.Seq import Seq

from Bio.Alphabet import IUPAC

from Bio import SeqIO

from Bio import motifs

import csv

# loaad meme result

filename1=raw_input('Enter MEME txt file (default = meme.txt): ')

if filename1=='':

	filename1='meme.txt'

handle=open(filename1) # Title of MEME file

record=motifs.parse(handle,'meme')

handle.close()

num = raw_input('Enter MEME motif number (default = 1)')

if num == '':

	num= 0
else: 
	num = int(num) - 1 # 1st motif=0
	
motif = record[num] 

print "\n.....Searching for motif:", motif.degenerate_consensus, "\n"

# load and parse fasta 

filename2= raw_input('Enter Fasta file name (default= test.fa): ') 

if filename2=='':

	filename2='test.fa'


fastafile = SeqIO.parse(open(filename2),'fasta')

# Can alternatively parse as Dict

# fastafile2=SeqIO.to_dict(SeqIO.parse('test.fa','fasta'))

# peak annotations

# annot_file='annotations.csv'

# annot_filehandle=csv.DictReader(open(annot_file, 'rb'), delimiter=',', quotechar='"')

# read in annotations as dict

# annot={}

# for line in annot_filehandle:

	# annot.update({line['Peak-ID']: line['Gene-Name']})

	
# prefix name for output file

prefix=raw_input('Enter output prefix (hit RETURN for default): ')

if prefix=='':
	
	prefix=  filename1.split('.')[0] + '_In' + filename2.split('.')[0]
	
	
# 2. PREPARE OUTPUT LISTS

motif_occurrences = []

header= list(["Gene ID", "Site Found", "Strand", "Location within 5' intergenic"]) 

motif_occurrences.append(header)

head_to_head_motifs = []

header2= list(["Gene ID", "Total No Sites", "Distance between half-sites"])

head_to_head_motifs.append(header2)


# 3. SEARCH FOR MEME MOTIF

for fasta in fastafile:

	# uppercase
	
	fasta.seq = fasta.seq.upper()
	
	# this is dna
	
	fasta.seq.alphabet=IUPAC.unambiguous_dna
	
	
	# check we got into the loop
	
	# print ".....Searching for motif in", fasta.id
	
	# search for fwd sites
	
	fwd_sites=0
	
	for pos1, seq1 in motif.instances.search(fasta.seq):
		
		fwd_sites +=1

		print ("%i %s" % (pos1, seq1)), "\t", fasta.id
		
		motif_entry=[fasta.id, str(seq1), '+', pos1]

		# print ".....Trying to add this motif entry: ", motif_entry
		
		motif_occurrences.append(motif_entry)
		
		# print ".....Trying to save this:", motif_occurrences
	# search for rev sites
	
	rc=motif.reverse_complement()

	rev_sites=0
	
	for pos2, seq2 in rc.instances.search(fasta.seq):
	
		rev_sites +=1
		
		print ("%i %s" % (pos2, seq2)), "\t", fasta.id
		
		motif_entry=[fasta.id, str(seq2), '-', pos2]
		
		# print ".....Trying to add this motif entry: ", motif_entry
		
		motif_occurrences.append(motif_entry)
		
# tally total sites and report distance between half sites

	total_sites = fwd_sites + rev_sites
	
	if fwd_sites > 0 and rev_sites> 0:
		
		head_to_head_entry = [fasta.id, total_sites, int(pos1 - pos2)]
		
		head_to_head_motifs.append(head_to_head_entry)
	
	
# 4.  WRITE TO OUTPUT FILE

print ".....Writing output files. Done."

with open('motif' + str(int(num+1)) + "_" + prefix + '_byPeak.csv', 'wb') as csvfile1:

	spamwriter= csv.writer(csvfile1, delimiter = ',')
	
	spamwriter.writerows(motif_occurrences)
	
with open('motif' + str(int(num+1)) + '_' + prefix + '_head2head.csv', 'wb') as csvfile:

	spamwriter= csv.writer(csvfile, delimiter = ',')
	
	spamwriter.writerows(head_to_head_motifs)