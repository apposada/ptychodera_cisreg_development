# ptychodera_cisreg_development
**Analyses of the Ptychodera flava RNASeq+ATACSeq project**

![logo](graphics/pfla.png?raw=true)

## About

This repository contains the code and most of the data for the analyses of the manuscript **Perez-Posada A, Lin CY, Lin ChingY, Chen YC, Gómez-Skarmeta JL, Yu JK, Su YH, Tena JJ: Insights into deuterostome evolution from the biphasic transcriptional programmes of hemichordates. bioRxiv 2022: https://doi.org/10.1101/2022.06.10.495707**.


## Organisation

The repository is mainly organised in:
 - code: the necessary code for the analyses. Importantly here:
	- `code/code_markdowns` contains the markdown files recapitulating all the tools and software used from mapping to generation of networks and orthology datasets.
	- `r_code/r_markdowns` contains the R markdown files for all the different analyses, therefore also including the code for generating all main and supplementary figures. All analyses were done in these R markdowns and the remaining R code is either legacy or used as supporting functions (see `code/r_code/r_functions`)
 - `data` contains the raw-mapped counts of RNA-seq, ATAC-seq. It also contains the full dataset of called consensus regions of open chromatin (OCRs).
 - `outputs` contains most of the outputs of the manuscript:
	- `outputs/rda` contains most of the .rda files from the R analysis that are generated and loaded by the different R markdowns. **Note: See details below**
	- `outputs/ananse` contains the input and output data for ANANSE. Most of these are omitted due to file size (fully available upon request)
	- `outputs/functional_annotation` contains a number of lookup tables with attributes deriving from transdecoder, eggNOG, Blast2GO and BLAST outputs to help identify TF genes and the identity of the genes of Ptychodera.
	- `outputs/comparative` contains most of the data from the reanalysis of other species, from OMA and Orthofinder, and from the Gene Age analyses. A number of these files have been omitted due to file size (fully available upon request)
	- `outputs/homer` contains the final output of the homer motif enrichment analyses and de novo motif discovery.
	- `outputs/reports_*` contain reports on the RNA-seq and ATAC-seq mapping.

## Notes

This repository is accompanied by a FigShare repository where to download other additional information: https://figshare.com/projects/Perez-Posada_et_al_Supplementary_Data/192362 . Please refer to the FigShare for supplementary files and R objects; see below.

This repository hosts a feature-freeze copy of the comparABle wrapper found at https://github.com/apposada/comparABle . Please refer to that repository for updates.

This repository also hosts a collection of functions and wrappers for graph analysis in both directed and non-directed graphs, such as those from ANANSE and WGCNA. These will likely become a separate repository in the future as well.

Most of these files and the necessary input are available in this repo, but one-click or full-code reproducibility might not be possible: some files have been omitted due to several kind of constraints such as file size limit, intermediate files, or because they are originally from a different publication. The most prominent example are the .rda and .RDS files that are needed to run the markdowns, and where all of the analyses were generated. These files are available in the figshare of this project.

If needed, they are fully available upon request. In that case, please contact Alberto Perez-Posada and/or Juan Tena, who can provide further information.

Any missing information is fully available upon request at `ap.posada1[AT]gmail[DOT]com`.
