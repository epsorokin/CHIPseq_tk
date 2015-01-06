==========
CHIPseq_tk
==========

Toolkit for ChIP-seq data analysis

**Annotation** 

(1) Annotating ChIP-seq peaks by nearest promoters: Identify any transcription start sites 
that occur within one kilobase of the midpoint of a ChIP peak in either direction. 
Report output with one annotation per line. 

* Usage: python annotate_bidirectional_promoters_v2.py 

(2) A second script performs the same task but spits out one line per peak, thus 
presenting the annotation information in a more condensed manner. 

* Usage: python annotate_bidirectional_promoters_one_row_per_peak.py

**Motif analysis**

(3) Motif searching within FASTA files using output from MEME as a starting point.

* Usage: 

python motif_searcher_v3.py # Search from within peak sequences

python global_motif_searcher.py # Search for within 5' intergenic regions or promoters

(4) Motif searching within FASTA files beginning with a user-defined motif

* Usage: 

python global_motif_searcher-user-defined-motif.py # Search from within peak sequences

python peak_motif_searcher-user-defined-motif.py # Search from within genomic sequences 
