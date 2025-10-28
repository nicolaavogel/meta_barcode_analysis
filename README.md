# barcodeMiner 

This is a repository for the barcode analysis regarding the Lake Tjoernin data of the preprint "Environmental DNA Reveals Reykjavík’s Human and Ecological History"; doi: https://doi.org/10.1101/2025.10.08.681091

barcodeMiner aims to analyse the proportion and probability of shotgun sequenced aeDNA reads that would hit taxonomically informative regions of the barcode sequence used in metabarcoding analysis. We expand this experiment to the full trnL gene sequence and compare it with the full chloroplast analysis. This repository contains databases, simulated data, and results used for assessing the performance and limitations of the trnL barcode region in plant DNA metabarcoding.
The project explores two main aspects:

## Database
The full databases are not provided in this repository (due to size limitations from GitHub), but can be sent when requested. We list the NCBI accession numbers contained in each database in accession lists stored in the repository. Additionally, we created lists of missing accession numbers for all the combinations between the databases, with additional taxonomic information, stored here: ```DBs
/viridiplantae_chloroplast/```. 

1. Full chloroplast DB – 12,816 complete chloroplast genomes from Viridiplantae. (downloaded using the ```script/down_DB.py```)

2. trnL gene DB – 11,808 trnL gene sequences extracted from the chloroplast DB. (downloaded using the ```script/down_trnl.py```)

3. trnL barcode DB – 10,851 barcode regions defined by the standard primers (extracted using rCRUX https://github.com/CalCOFI/rCRUX ```script/rcrux.R```)
(forward: GGGCAATCCTGAGCCAA, reverse: CCATTGAGTCTCTGCACCTATC).

4. Mapping and LCA analysis were done using bowtie2 and metaDMG: ```script/bar_miner_pip.sh``` or ```script/bar_miner_pip_emp.sh```

## Data simulation

To ensure traceable results, we simulated data using gargammel and ART Illumina: ```scripts/run_sim.sh```. We simulated ancient damage probability and the fragment length distribution after Günther, Torsten, et al. "Ancient genomes link early farmers from Atapuerca in Spain to modern-day Basques." Proceedings of the National Academy of Sciences 112.38 (2015): 11917-11922; doi: https://doi.org/10.1073/pnas.1509851112. The deamination files can be found in the ```sim_data``` folder starting with ```dhigh```. 
All simulated fastq files used for the analysis in this paper can be found in the ```sim_data``` folder. The results of the simulations to the barcode and full trnL gene database can be found in the ```results``` folder. 

## Empirical data experiments
The raw empirical data will be uploaded to ENA (identification number will follow). As described in the publication, we used all reads assigned to Viridiplantae and remapped them to our specific databases. The results of that mapping to the full trnL gene and the barcode database can be found in the ```emp_results``` folder. The results files for the full chloroplast database could not be uploaded due to the GitHub size limitation, but can be send when requested. 

## Data visualisation 
All data was visualised using custom R scripts: ```scripts/plotting.R```. Plots from the paper can be found in the ```plots``` folder. 



