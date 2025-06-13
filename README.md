# PEDIA-BRAIN: A single nuclei multiomic encyclopedia of the human pons provides a resource for normal development and disease vulnerability

This is a repository for the manuscript entitled "PEDIA-BRAIN: A single nuclei multiomic encyclopedia of the human pons provides a resource for normal development and disease vulnerability". Available here (link here).

Authors: Tianli Ding†, Gary Schweickart†, Adithe Rivaldi, Cole Harrington,  Aaron Varghese, Ke Qin, Ben Kelly, Grant Lammi, Benjamin D. Sunkel, Kathryn L. Stahl, Alex H. Wagner, Jeffrey R. Leonard, Albert M. Isaacs, Katherine E. Miller, Elaine R. Mardis, Michelle A. Wedemeyer*

†These authors contributed equally to this work

\* Lead contact: 700 Children’s Dr, Columbus, OH 43205

## About
Advancements in scalable and reliable single cell sequencing technologies have significantly enhanced our understanding of the human brain's diversity in terms of gene expression and chromatin accessibility. Numerous research groups have mapped the cellular diversity of the fetal brain with unprecedented detail. However, most studies on the transcriptional and chromatin signatures of the developing human brain have primarily concentrated on the telencephalon or limited developmental timeframes.

The PEDIA-BRAIN (Pediatric Encyclopedia of Development from Infancy to Adulthood - Brain Regions and Associated Integrative Neuroscience) is a public resource encompassing the complete spectrum of human brain development. It incorporates droplet-based single nuclei RNA sequencing and open chromatin profiling (snRNA-seq + snATAC-seq) of normal tissues corresponding to critical stages of postnatal development including infant and toddler (0-2 years), child (4-8 years), adolescent (13-14 years), and early adulthood, at the completion of myelination (31-33 years).

For a user-friendly analysis tool, check out https://pediabrain.nchgenomics.org/.

## Structure of the Repository

* **analysis** contains all code used to create figures an tables 
* **code** contains functions relevant to the creation of the pons atlas
* **data** contains the raw data files (not included see GEO)
* **output** contains figures, tables, and RDS files created during the analysis steps

## Creating Figures/Tables

* Download Pediatric and Adult files from GEO (Link)
* Download Fetal pons data by following the instructions here: https://github.com/linnarsson-lab/fetal_brain_multiomics?tab=readme-ov-file 
* Follow MultiomeWorkflow.Rmd, in the **code** folder, for each project
* Merge Projects following Merge_Projects.Rmd, in the **code** folder.
* Use created atlas for subsequent analysis files in the **analysis** folder

Note: You can skip the creation of the atlas by downloading the full or downsampled atlas here: https://pediabrain.nchgenomics.org/
