#############################################################
#### Locate motifs by position relative to peak center   ####
####    by Elena P Sorokin   March 4, 2015               ####
#############################################################

'''
Search for  motif within a fasta file.
Find motif locations, relative to the peak center.
Return csv file that can be directly used to make a histogram.
'''

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import motifs
import csv


# 2. CHOOSE MOTIF 

# a, dreme-based LBS motif, SGTGGAA

a= [Seq("CGTGGGAA"), Seq("GGTGGGAA")]
amotif=motifs.create(a)

# b, dreme-based unknown motif TACNGTA

b= [Seq("TACAGTA"), Seq("TACGGTA"), Seq("TACCGTA")]
bmotif=motifs.create(b)

# c, LBS from Christensen and Kimble 1996
c=[Seq("ATGGGAA"), Seq("GTGGGAA")]
cmotif=motifs.create(c)

# d, LBS from Yoo and Greenwald 

d=[Seq("CATGAGAA"), Seq("CGTGGGAA"), Seq("CATGGGAA"), Seq("CGTGAGAA"), Seq("TATGAGAA"), Seq("TGTGGGAA"), Seq("TATGGGAA"), Seq("TGTGAGAA")]
dmotif=motifs.create(d)
	
# 3. PARSE INPUT FASTA AND PEAK ANNOTATIONS

filename2= raw_input('Enter peak sequences file (default= partB_lag_peaks+1000bp.fa): ') 

if filename2=="":

	filename2='partB_lag_peaks+1000bp.fa'

fastafile = SeqIO.parse(open(filename2),'fasta')

# Peak annotations -- only necessary if they're not in the original fasta file

annot_file='annotations.csv'

annot_filehandle=csv.DictReader(open(annot_file, 'rb'), delimiter=',', quotechar='"')

# Read in annotations as hash table

annot={}

for line in annot_filehandle:

	annot.update({line['Peak-ID']: line['Gene-Name']})
	
# prefix name for output file

prefix=raw_input('Enter output prefix (hit RETURN for default): ')

if prefix=='':
	
	prefix=  'hist_' + filename2.split('.')[0]
	
group = raw_input("Enter group name, e.g. LAG.1:")	
if group=="":
	group="LAG.1"
	
# 4. PREPARE OUTPUT LISTS

motif_occurrences = []

header= list(["Peak or Gene ID", "Gene Name", "Motif", "Group","String Found", "Strand", "Distance from peak center"]) 

motif_occurrences.append(header)

# 5. SEARCH FOR MOTIF

for fasta in fastafile:

	# Uppercase
	fasta.seq = fasta.seq.upper()
	
	# This is DNA
	fasta.seq.alphabet=IUPAC.unambiguous_dna
	
	# Search for motif instances
	fwd_sites=0
	
	for pos1, seq1 in amotif.instances.search(fasta.seq):
				
		center = len(fasta.seq)/2
		
		fwd_sites +=1
		
		fasta.name=annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos1-center, seq1))  
		
		motif_entry=[fasta.id, fasta.name, amotif.degenerate_consensus, group,str(seq1), '+', pos1- center]
	
		motif_occurrences.append(motif_entry)
		
	# SEARCH FOR REVERSE COMPLEMENT
	
	rc=amotif.reverse_complement()

	rev_sites=0
	
	for pos2, seq2 in rc.instances.search(fasta.seq):
	
		rev_sites +=1
		
		center = len(fasta.seq)/2
		
		fasta.name= annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos2, seq2))
		
		motif_entry=[fasta.id, fasta.name, amotif.degenerate_consensus, group,str(seq2), '-', pos2-center]
		
		motif_occurrences.append(motif_entry)

# 6. LOOP B - reparse FASTA 		

fastafile = SeqIO.parse(open(filename2),'fasta')

for fasta in fastafile:

	# Uppercase
	fasta.seq = fasta.seq.upper()
	
	# This is DNA
	fasta.seq.alphabet=IUPAC.unambiguous_dna
	
	# Search for motif instances
	fwd_sites=0
	
	for pos1, seq1 in bmotif.instances.search(fasta.seq):
				
		center = len(fasta.seq)/2
		
		fwd_sites +=1
		
		fasta.name=annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos1-center, seq1))  
		
		motif_entry=[fasta.id, fasta.name, bmotif.degenerate_consensus, group,str(seq1), '+', pos1- center]
	
		motif_occurrences.append(motif_entry)
		
	# SEARCH FOR REVERSE COMPLEMENT
	
	rc=bmotif.reverse_complement()

	rev_sites=0
	
	for pos2, seq2 in rc.instances.search(fasta.seq):
	
		rev_sites +=1
		
		center = len(fasta.seq)/2
		
		fasta.name= annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos2, seq2))
		
		motif_entry=[fasta.id, fasta.name, bmotif.degenerate_consensus, group,str(seq2), '-', pos2-center]
		
		motif_occurrences.append(motif_entry)
			
# 7. LOOP C - reparse FASTA 		

fastafile = SeqIO.parse(open(filename2),'fasta')

for fasta in fastafile:

	# Uppercase
	fasta.seq = fasta.seq.upper()
	
	# This is DNA
	fasta.seq.alphabet=IUPAC.unambiguous_dna
	
	# Search for motif instances
	fwd_sites=0
	
	for pos1, seq1 in cmotif.instances.search(fasta.seq):
				
		center = len(fasta.seq)/2
		
		fwd_sites +=1
		
		fasta.name=annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos1-center, seq1))  
		
		motif_entry=[fasta.id, fasta.name, cmotif.degenerate_consensus, group,str(seq1), '+', pos1- center]
	
		motif_occurrences.append(motif_entry)
		
	# SEARCH FOR REVERSE COMPLEMENT
	
	rc=cmotif.reverse_complement()

	rev_sites=0
	
	for pos2, seq2 in rc.instances.search(fasta.seq):
	
		rev_sites +=1
		
		center = len(fasta.seq)/2
		
		fasta.name= annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos2, seq2))
		
		motif_entry=[fasta.id, fasta.name, cmotif.degenerate_consensus, group,str(seq2), '-', pos2-center]
		
		motif_occurrences.append(motif_entry)			
			
# 8. LOOP D - reparse FASTA 		

fastafile = SeqIO.parse(open(filename2),'fasta')

for fasta in fastafile:

	# Uppercase
	fasta.seq = fasta.seq.upper()
	
	# This is DNA
	fasta.seq.alphabet=IUPAC.unambiguous_dna
	
	# Search for motif instances
	fwd_sites=0
	
	for pos1, seq1 in dmotif.instances.search(fasta.seq):
				
		center = len(fasta.seq)/2
		
		fwd_sites +=1
		
		fasta.name=annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos1-center, seq1))  
		
		motif_entry=[fasta.id, fasta.name, dmotif.degenerate_consensus, group,str(seq1), '+', pos1- center]
	
		motif_occurrences.append(motif_entry)
		
	# SEARCH FOR REVERSE COMPLEMENT
	
	rc=dmotif.reverse_complement()

	rev_sites=0
	
	for pos2, seq2 in rc.instances.search(fasta.seq):
	
		rev_sites +=1
		
		center = len(fasta.seq)/2
		
		fasta.name= annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos2, seq2))
		
		motif_entry=[fasta.id, fasta.name, dmotif.degenerate_consensus, group,str(seq2), '-', pos2-center]
		
		motif_occurrences.append(motif_entry)			
									
# 6.  WRITE OUTPUT

with open( prefix + '_sites.csv', 'wb') as csvfile:
	spamwriter= csv.writer(csvfile, delimiter = ',')
	spamwriter.writerows(motif_occurrences)
	
print "\n.....Writing output file. Done."