# 2019_TRiC_Hsp70_cotranslationalFolding

This repository contains the code and annotation files used to examine the interactions of TRiC and the Hsp70 Ssb with the nascent proteome as published in:

Stein KC, Kriel A, Frydman J. Nascent Polypeptide Domain Topology and Elongation Rate Direct the Cotranslational Hierarchy of Hsp70 and TRiC/CCT. Mol Cell. 2019 Sep 19;75(6):1117-1130.e5. doi: 10.1016/j.molcel.2019.06.036. Epub 2019 Aug 7. [PMID: 31400849](https://pubmed.ncbi.nlm.nih.gov/31400849/)

Raw data and pre-processed codon counts tables can be downloaded from GEO: [GSE114882](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSE114882)


## Post_Processing directory:

**PostProcessing.R:** Import codon counts tables and calculate relevant metrics by position and gene.

**FishersPeaks.R:** Identify positions of cotranslational chaperone interactions using Fisher's exact test and isolate peak position with five ore more consecutive residues with significant chaperone enrichment.

**LiteratureComparison.R:** Parse chaperone binding data from previous work.


## Domain_Analysis directory:

**DomainAnalysis.R:** Analyze the proximity of chaperone binding sites relative to protein domains.

**WDanalysis.R:** Prepare table for analysis around WD repeat elements.


## Gene-level analysis code (Gene_Enrichment directory):

**DESeq.R:** Uses DESeq2 to determine significantly enriched genes from selective ribosome profiling data. 

**ProteomeProperties.R:** Analyze enriched protein properties of chaperone substrates.


## Codon-level analysis code around chaperone binding sites (Peptide_Analysis directory):

**PeptideProperties.R:**
Characterize the secondary structure, aggregation propensity, hydrophobicity, charge, and other miscellaneous biophysical propeties of the nascent polypeptide.

**TranslationRate.R:**
Analyze ribosome occupancy as a proxy local translation rates around chaperone binding sites.

**SequenceAnalysis.R:** 
Analyze the amino acid composition and codon optimality of the nascent polypeptide.


## Other relevant code:

**Annotation directory:** Relevant annotation files, summary tables, and data used for secondary structure analysis (psipred) and amino acid composition (sequenceAnalysis).

**Figures directory:** Specific visualization code for generating plots included in publication.

**MovingAverage.R:** function to average ribosome occupancy over user-defined window.

