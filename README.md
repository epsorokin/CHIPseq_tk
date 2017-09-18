CHIPseq_tk
=========

An all-purpose toolkit for downstream analyses of peaks derived from a ChIP-Seq study.

* Annotation of ChIP peaks 
* Motif analysis

Peak Annotation

Annotating ChIP-seq peaks by nearest promoters: Identify any transcription start sites 
that occur within one kilobase of the midpoint of a ChIP peak in either direction. 
Report output with one annotation per line. 

`Usage: python annotate_bidirectional_promoters_v2.py `

A second script performs the same task but splits out one line per peak, thus 
presenting the annotation information in a more condensed manner. 

`Usage: python annotate_bidirectional_promoters_one_row_per_peak.py`

Motif Analysis 

Motif searching within FASTA files using output from `MEME` as a starting point.

`Usage: `

`python motif_searcher_v3.py # Search from within peak sequences`

`python global_motif_searcher.py # Search for within 5' intergenic regions or promoters`

Motif searching within FASTA files beginning with a user-defined motif. The first script searches from within peak sequences, the second script searches from across all genomic sequences for any pre-specified genome.

`Usage: `

`python global_motif_searcher-user-defined-motif.py`

`python peak_motif_searcher-user-defined-motif.py`

