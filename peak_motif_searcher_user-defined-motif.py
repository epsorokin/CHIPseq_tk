#############################################################
#### Motif searching tool for FASTA sequences   #############
####    by Elena P Sorokin   December 16, 2014  #############
#############################################################

'''Search for user-defined motif within a FASTA file.
Return # motif occurrences and motif location(s).'''



# 1. IMPORT MODULES, CHOOSE MOTIF, FASTA TO SEARCH

# modules

from Bio.Seq import Seq

from Bio.Alphabet import IUPAC

from Bio import SeqIO

from Bio import motifs

import csv

# determine motif

instances= [Seq("CATGAGAA"), Seq("CGTGGGAA"), Seq("CATGGGAA"), Seq("CGTGAGAA"),
Seq("TATGAGAA"), Seq("TGTGGGAA"), Seq("TATGGGAA"), Seq("TGTGAGAA")]

m= motifs.create(instances)

print m.counts

print m.consensus

print m.degenerate_consensus

print "\n.....Searching for motif:", m.degenerate_consensus, "\n"

# load and parse fasta 

filename2= raw_input('Enter Fasta file name (default= test.fa): ') 

if filename2=='':

	filename2='test.fa'


fastafile = SeqIO.parse(open(filename2),'fasta')

	
# prefix name for output file

prefix=raw_input('Enter output prefix (hit RETURN for default): ')

if prefix=='':
	
	prefix=  'LBS_in_' + filename2.split('.')[0]
	
	
# 2. PREPARE OUTPUT LISTS

motif_occurrences = []

header= list(["Gene ID", "Gene Name", "Site Found", "Strand", "Location within 5' intergenic"]) 

motif_occurrences.append(header)


motif_summary = []

header2= list(["Gene ID", "Gene Name", "Total No Sites", "Locations"])

motif_summary.append(header2)


# 3. SEARCH FOR MEME MOTIF

for fasta in fastafile:

	# uppercase
	
	fasta.seq = fasta.seq.upper()
	
	# this is dna
	
	fasta.seq.alphabet=IUPAC.unambiguous_dna
	
	
	# search for fwd sites
	
	motif_locations = []
	
	fwd_sites=0
	
	for pos1, seq1 in m.instances.search(fasta.seq):
		
		fwd_sites +=1
	
		
		print ("%i %s" % (pos1, seq1)), "\t", fasta.id
		
		motif_entry=[fasta.id.split("|")[0], fasta.id.split("|")[1], str(seq1), '+', pos1]

		# print ".....Trying to add this motif entry: ", motif_entry
		
		motif_occurrences.append(motif_entry)
		
		motif_locations.append(pos1)
		
	# search for rev sites
	
	rc=m.reverse_complement()

	rev_sites=0
	
	for pos2, seq2 in rc.instances.search(fasta.seq):
	
		rev_sites +=1
		
		print ("%i %s" % (pos2, seq2)), "\t", fasta.id
		
		motif_entry=[fasta.id.split("|")[0], fasta.id.split("|")[1], str(seq2), '-', pos2]
				
		motif_occurrences.append(motif_entry)
		
		motif_locations.append(pos2)
		
# tally total sites and report their locations

	total_sites = fwd_sites + rev_sites
	
	if total_sites > 0:
		
		entry = [fasta.id.split("|")[0], fasta.id.split("|")[1], total_sites, motif_locations]
		
		motif_summary.append(entry)
	
	# tally total sites and report distance between half sites

	total_sites = fwd_sites + rev_sites
	
	if fwd_sites > 0 and rev_sites> 0:
		
		head_to_head_entry = [fasta.id, total_sites, int(pos1 - pos2)]
		
		head_to_head_motifs.append(head_to_head_entry)
# # tally total sites and report distance between half sites

	# total_sites = fwd_sites + rev_sites
	
	# if fwd_sites > 0 and rev_sites> 0:
		
		# head_to_head_entry = [fasta.id, total_sites, int(pos1 - pos2)]
		
		# head_to_head_motifs.append(head_to_head_entry)
		
# 4.  WRITE TO OUTPUT FILE

print ".....Writing output files. Done."

with open( prefix + '_MotifOccurrences.csv', 'wb') as csvfile1:

	spamwriter= csv.writer(csvfile1, delimiter = ',')
	
	spamwriter.writerows(motif_occurrences)
	
with open(prefix + '_MotifSummary.csv', 'wb') as csvfile:

	spamwriter= csv.writer(csvfile, delimiter = ',')
	
	spamwriter.writerows(motif_summary)