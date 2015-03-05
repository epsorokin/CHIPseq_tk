#############################################################
#### Motif searching tool for FASTA sequences   #############
####    by Elena P Sorokin   December 16, 2014  #############
#############################################################

'''
Search for user-defined motif within a peak fastafile.
Get # motif occurrences and motif location(s), both by peak 
and as a summary file.
'''

# 1. IMPORT MODULES 

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import motifs
import csv

# 2. CHOOSE MOTIF 

# m, Yoo and Greenwald motif

yoo= [Seq("CATGAGAA"), Seq("CGTGGGAA"), Seq("CATGGGAA"), Seq("CGTGAGAA"), Seq("TATGAGAA"), Seq("TGTGGGAA"), Seq("TATGGGAA"), Seq("TGTGAGAA")]

y= motifs.create(yoo)

print "\n.....Searching for Yoo and Greenwald motif:", y.degenerate_consensus, "\n"
print y.counts

# k, Christiansen and Kimble motif

kimble= [Seq("ATGGGAA"), Seq("GTGGGAA")]

k= motifs.create(kimble)

# 3. PARSE INPUT FASTA AND PEAK ANNOTATIONS

filename2= raw_input('Enter peak Fasta file (default= test.fa): ') 

if filename2=='':

	filename2='test.fa'

fastafile = SeqIO.parse(open(filename2),'fasta')

# peak annotations

annot_file='annotations.csv'

annot_filehandle=csv.DictReader(open(annot_file, 'rb'), delimiter=',', quotechar='"')

# read in annotations as dict

'''
This next section  is specifically for peak fasta files.
'''

annot={}

for line in annot_filehandle:

	annot.update({line['Peak-ID']: line['Gene-Name']})
	
# prefix name for output file

prefix=raw_input('Enter output prefix (hit RETURN for default): ')

if prefix=='':
	
	prefix=  'NEW_lbs_in_' + filename2.split('.')[0]
	
# 4. PREPARE OUTPUT LISTS

motif_occurrences = []

header= list(["Peak or Gene ID", "Gene Name", "Site Found", "Strand", "Location within 5' intergenic"]) 

motif_occurrences.append(header)

summary = []

summaryHeader= list(["Peak or Gene ID",  "Gene Name", "Total No Sites"])

summary.append(summaryHeader)

# 5. SEARCH FOR YOO MOTIF

runningTotal = lbs_1 = lbs_2 = lbs_3= lbs_4= lbs_5orMore = 0

for fasta in fastafile:

	# uppercase
	fasta.seq = fasta.seq.upper()
	
	# this is dna
	fasta.seq.alphabet=IUPAC.unambiguous_dna
	
	# search for fwd sites
	fwd_sites=0
	
	for pos1, seq1 in y.instances.search(fasta.seq):
		
		fwd_sites +=1

		fasta.name=annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos1, seq1))  
		
		motif_entry=[fasta.id, fasta.name, str(seq1), '+', pos1]
	
		motif_occurrences.append(motif_entry)
		

	# search for rev sites
	rc=y.reverse_complement()

	rev_sites=0
	
	for pos2, seq2 in rc.instances.search(fasta.seq):
	
		rev_sites +=1
		
		fasta.name= annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t", ("%i %s" % (pos2, seq2))
		
		motif_entry=[fasta.id, fasta.name, str(seq2), '-', pos2]
		
		motif_occurrences.append(motif_entry)
		
# Still in loop, tally total sites and report distance between half sites
	total_sites = fwd_sites + rev_sites
	
	entry = [fasta.id, fasta.name, total_sites ]
		
	summary.append(entry)
	
	# update global totals
	if total_sites > 0:
		runningTotal +=1
		
		# Then, update individual totals
		print "updating individual totals"
		if total_sites == 1:
			lbs_1 += 1
		elif total_sites == 2:
			lbs_2 += 1
		elif total_sites ==3:
			lbs_3 += 1	
			print "\n*****Found 3 sites in ", fasta.name, "\n"
		elif total_sites ==4:
			lbs_4 +=1
			print "\n*****Found 4 sites in ", fasta.name, "\n"
		elif total_sites>= 5:
			lbs_5orMore += 1
			print "\n*****Found five or more sites in ", fasta.name, "\n"
		else:
			print "Error! Counting tally failed", fasta.name

# 6. PRINT YOO SUMMARY

print "\n\n.....Found the following YRTGRGAA (Yoo site) in", filename2, "\n"
print "TOTAL SITES:\t", runningTotal
print "One:\t\t", lbs_1
print "Two:\t\t", lbs_2
print "Three:\t\t", lbs_3
print "Four:\t\t ", lbs_4
print "5 or more: ", lbs_5orMore, "\n"

# 7.  WRITE TO OUTPUT FILES

print "\n.....Writing output files. Done."

# Not writing this file (for now)
#with open( prefix + 'Yoo_byPeak.csv', 'wb') as csvfile1:
	#spamwriter= csv.writer(csvfile1, delimiter = ',')
	#spamwriter.writerows(motif_occurrences)
	
with open(prefix + '_Yoo_summary.csv', 'wb') as csvfile:

	spamwriter= csv.writer(csvfile, delimiter = ',')
	
	spamwriter.writerows(summary)
	
# 8. RE-PARSE FASTA FILE 

fastafile = SeqIO.parse(open(filename2),'fasta')

# 9. SEARCH FOR KIMBLE MOTIF IN FASTA

print "\n.....Searching for Christiansen and Kimble motif: ", k.degenerate_consensus,"\n"

print k.counts

runningKimbleTotal = k1 = k2 = k3= k4= k5orMore = 0

# 10. PREPARE SECOND SET OF OUTPUT LSITS USING HEADERS ABOVE

motif_occurrences2 = []

motif_occurrences2.append(header) #header from Yoo 

summary2 = []

summary2.append(summaryHeader) #header from earlier

# 11. LOOP THROUGH FASTA AGAIN 

for fasta in fastafile:

	# uppercase
	fasta.seq = fasta.seq.upper()
	
	# this is dna	
	fasta.seq.alphabet=IUPAC.unambiguous_dna

	# search for fwd sites
	fwd_sites=0
	
	for pos1, seq1 in k.instances.search(fasta.seq):
		
		fwd_sites +=1

		fasta.name=annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t",("%i %s" % (pos1, seq1))
		
		motif_entry=[fasta.id, fasta.name, str(seq1), '+', pos1]
	
		motif_occurrences2.append(motif_entry)

	# search for rev sites
	
	rc=k.reverse_complement()

	rev_sites=0
	
	for pos2, seq2 in rc.instances.search(fasta.seq):
	
		rev_sites +=1
		
		fasta.name= annot[fasta.id]
		
		print fasta.id, "\t", fasta.name, "\t",("%i %s" % (pos2, seq2))
		
		motif_entry=[fasta.id, fasta.name, str(seq2), '-', pos2]
		
		motif_occurrences2.append(motif_entry)
		
# tally total sites and report distance between half sites
	total_sites = fwd_sites + rev_sites
	
	entry = [fasta.id, fasta.name, total_sites ]

	summary2.append(entry)
	
	# update global totals
	
	if total_sites > 0:
		runningKimbleTotal +=1
		
		# update individual totals
		if total_sites == 1:
			k1 += 1
		
		elif total_sites == 2:
			k2 += 1
		
		elif total_sites ==3:
			k3 += 1
			print "\n*****Found 3 sites in ", fasta.name, "\n"
	
		elif total_sites ==4:	
			k4 +=1
			print "\n*****Found 4 sites in ", fasta.name, "\n"
			
		elif total_sites>= 5:
			k5orMore += 1
			print "\n*****Found five or more sites in ", fasta.name, "\n"
			
		else:
			print "Error! Counting tally failed", fasta.name

# 12. PRINT YOO SUMMARY AGAIN

print "\n\n.....Found the following", y.degenerate_consensus, " (Yoo site) in", filename2, "\n"
print "TOTAL SITES:\t", runningTotal
print "One:\t\t", lbs_1
print "Two:\t\t", lbs_2
print "Three:\t\t", lbs_3
print "Four:\t\t ", lbs_4
print "5 or more: ", lbs_5orMore
			
# 13. PRINT KIMBLE SUMMARY

print "\n\n.....Found the following", k.degenerate_consensus,"(Kimble site) in", filename2, "\n"
print "TOTAL SITES:\t", runningKimbleTotal
print "One:\t\t", k1
print "Two:\t\t", k2
print "Three:\t\t", k3
print "Four:\t\t ", k4
print "5 or more: ", k5orMore


# 14.  WRITE TO OUTPUT FILES

# Not writing this file (for now)

#with open( prefix + '_byPeak.csv', 'wb') as csvfile1:

	#spamwriter= csv.writer(csvfile1, delimiter = ',')
	
	#spamwriter.writerows(motif_occurrences)
	
with open(prefix + '_Kimble_summary.csv', 'wb') as csvfile:

	spamwriter= csv.writer(csvfile, delimiter = ',')
	
	spamwriter.writerows(summary2)