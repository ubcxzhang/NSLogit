> The code assumes that all files are placed under a same working directory on Compute Canada.  

### Directory Layout
    .
    ├── conquer_framework  
    │   ├── create_subfile.R      # Script for generating conquer subsets           
    │   ├── experiment.R          # Script for reproducing conquer experiments
    │   ├── merge_sub_result.R    # Script for merging results for pre_plot_sim123.R
    │   ├── pre_plot_sim123.R     # Combine NsLogit&LR results with conquers' and format them for plot
    │   └── plot_sim123.R         # FDR&TPR&F1 for each method on conquer data
    │
    ├── scDD_framework                
    │   ├── simulation_scDD.R     # Function to do simulation with scDD 
    │   ├── simulation_script.R   # Script to generate simulated data from scDD
    │   ├── simu_eval.R           # Evaluate DE methods using simulated data from scDD
    │   ├── simu_eval_comb.R      # Script to merge results for pre_scDD_simu_plot.R
    │   ├── pre_scDD_simu_plot.R  # Format scDD DE results files for plot
    │   └── scDD_simu_plot.R      # FDR&TPR&F1 for each method on scDD simulated data
    │
    ├── real_data
    │   ├── partition_fullsize.R  # Partition full size GSE135893 according to cell types
    │   ├── GSE135893_DE.R        # Evaluate DE methods using GSE135893
    │   ├── real_data_analysis.R  # Summarize real data results
    │   ├── Differentiating_Ciliated.R # Plot pairwise DE compare on Differentiating Ciliated cells
    │   └── Seurat_NB.R           # Details of the method used in GSE135893 published work (not in the pipeline)
    │
    ├── real_data_GSE136831
    │   ├── partition_fullsize_GSE136831.R  # Partition full size GSE136831 according to cell types
    │   ├── analysis_real_data_res.R        # R script for combining paper results and ns results for later anaylsis (for both 135 and 136)
    │   └── GSE136831_DE.R                  # Evaluate DE methods using GSE136831
    │
    ├── DE_methods                # Define all DE methods in this study
    ├── computecanada             # All shell scripts for job submission on Compute Canada
    ├── load_pkg.R                # Script for testing whether user can load required packages with correct version
    ├── utils.R                   # Define helper functions used in other R scripts
    └── README.md

### Folder structure on Compute Canada

    . /home/xsslnc/projects/def-ubcxzh/xsslnc/ # the variable `work_dir`
    ├── MasterProject             # Up-to-date GitHub Repository
    ├── D3E                       # D3E code (download from https://github.com/hemberg-lab/D3E)
    ├── 3.5                       # R 3.5.0 library with required dependencies
    ├── tmp                       # D3E python output
    ├── ENV                       # D3E python environment library
    ├── py_requirements.txt       # D3E python environment requirements
    ├── conquer_data              # Original data download from conquer database (unzipped from http://imlspenticton.uzh.ch/robinson_lab/conquer_de_comparison/Soneson_datasets.tar.gz)
    ├── conquer_sim123            # Conquer's published results file (unzipped from http://imlspenticton.uzh.ch/robinson_lab/conquer_de_comparison/Soneson_results.tar.gz)
    ├── sim123_filtered           # Results from reproducing conquer's filtered experiments on NsLogit and LR
    ├── sim123                    # Results from reproducing conquer's unfiltered experiments on NsLogit and LR
    ├── GSE135893_count           # Processed GSE135893 (counts slot of SCT for each cell type)
    ├── GSE135893_with_tpm        # Data file from GSE135893_count with tpm 
    ├── GSE135893_file
    │   ├── GSE135893_ILD_annotated_fullsize.rds # Full size GSE135893
    │   ├── NsLogit_count_result  # Results of applying NsLogit on GSE135893_count
    │   ├── Disease_vs_Control    # Excel results files from published paper (Table S5): https://advances.sciencemag.org/content/6/28/eaba1972/tab-figures-data
    │   └── ...
    ├── GSE136831_file
    │   ├── NsLogit_count_result  # Results of applying NsLogit on GSE136831_count
    │   ├── paper_result          # media-8.txt lists DE genes of IPF vs Control, the folder contains all txt results files from published paper: https://www.biorxiv.org/content/10.1101/759902v3.supplementary-material    
    │   └── ...                   # Data download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831
    ├── Figure                    # Where plots will be placed
    ├── real_data_comb_res        # Folder containing combined paper & ns results for each cell types
    ├── simu_data                 # Data simulated by scDD
    └── scDD_simu_result          # Where results on `simu_data` will be stored

Users need to check whether their working directory (`work_dir`) matched the above structure before job submission.

---

- **IMPORTANT**: following script will ask user for the full path containing all files (e.g. /home/xsslnc/projects/def-ubcxzh/xsslnc for me)
- The input should not end with `/`

---

### Reproduce conquer experiments using new method

- Generate conquer subfile

As known, conquer subsets each data 16 times so that they run 16 experiments with one simulated data. Modified from [their script](https://github.com/csoneson/conquer_comparison/blob/master/scripts/generate_subsets.R), `create_subfile.R` read data from `conquer_data` and create subfile for each, and save subfiles in `conquer_data` (overwrite existing ones).

```R
Rscript create_subfile.R "/home/xsslnc/projects/def-ubcxzh/xsslnc"
```

- Generate results for NsLogit and LR on conquer data

Use following script to generate results (filtered data & unfiltered data) for **NsLogit** and **LR**. Each script take input data from `conquer_data`.

```shell
./run_GSE74596sim123_sub.sh
./run_GSE60749-GPL13112sim123_sub.sh
./run_GSE45719sim123_sub.sh
```

When all results are available in `sim123_filtered` and `sim123`, run the merge script to combine them for plot:

```shell
./run_merge_sub_result.sh
```

- Merge and format results together with conquer's

```shell
./run_pre_plot_sim123.sh
```
It will create `sim123_res.rds` under `sim123`, and `sim123_res_TPM_1_25p.rds` under `sim123_filtered`, which are the formatted full results (conquer's, NsLogit and LR).

- Plot

When `sim123_res.rds` and `sim123_res_TPM_1_25p.rds` are ready, by running `plot_sim123.R`, `no_filt.png` and `filt.png` will be available under `/home/xsslnc/projects/def-ubcxzh/xsslnc/Figure/`.

---

### Simulation study with scDD

- Check R dependencies

By sroucing `load_pkg.R` (add working directory as input parameter), user can check required packages' version and whether they can be loaded.

```R
Rscript load_pkg.R "/home/xsslnc/projects/def-ubcxzh/xsslnc"
```
The script will print `sessionInfo`, which should be the same as following:

```
other attached packages:
 [1] scDD_1.6.1                  metagenomeSeq_1.24.1
 [3] RColorBrewer_1.1-2          glmnet_2.0-18
 [5] edgeR_3.24.3                limma_3.38.3
 [7] DESeq2_1.14.1               samr_2.0
 [9] impute_1.48.0               MAST_1.8.2
[11] SingleCellExperiment_1.4.1  SummarizedExperiment_1.12.0
[13] DelayedArray_0.8.0          BiocParallel_1.16.6
[15] matrixStats_0.55.0          GenomicRanges_1.34.0
[17] GenomeInfoDb_1.18.2         IRanges_2.16.0
[19] S4Vectors_0.20.1            BPSC_0.99.1
[21] doParallel_1.0.15           iterators_1.0.12
[23] foreach_1.4.7               statmod_1.4.32
[25] DEsingle_1.2.1              scde_1.99.1
[27] flexmix_2.3-13              lattice_0.20-38
[29] Seurat_2.2.1                cowplot_1.0.0
[31] monocle_2.10.1              DDRTree_0.1.4
[33] irlba_2.3.3                 VGAM_1.1-2
[35] ggplot2_3.2.1               Biobase_2.42.0
[37] BiocGenerics_0.28.0         Matrix_1.2-14
[39] ROTS_1.10.1

```

- Simulate scDD data

This script submits 54 jobs, each for simulating one set from scDD. The input data are in `conquer_data`, which are original conquer data.
```shell
./run_simulation_script.sh
```

- Apply DE methods on simulated data

When the simulation is completed, simulated data will be available in `simu_data`. Then user can run the script for applying DE methods on the simulated sets (Include methods to evaluate in the variable `mt`. `run_simu_eval.sh` lists all of them).

```shell
./run_simu_eval.sh
```

Before running `./run_simu_eval.sh`, user need to make sure that folders named after those methods to be run are available under `scDD_simu_result`. The DE results for each method will appear under corresponding folder (e.g. `scDD_simu_result/NsLogit`).

#### Time estimates for scDD experiments

The following table lists percentage of experiments that can be finished within a time limit for each method.

| Method   | 0 ~ 1 h | 1 ~ 3 h | 3 ~ 5 h | 5 h + | max running time (hours) |
| ---      | ---     | ---     | ---     | ---   | --    |
| DEsingle       | /       | 25.00%  | /       | 75.00%   | 59.70 |
| edgeRLRT       | 96.30%  | 3.70%   | /       | /        | 1.16  |
| edgeRLRTcensus | 85.19%  | 14.81%  | /       | /        | 1.26  |
| edgeRLRTdeconv | 94.44%  | 5.56%   | /       | /        | 1.13  |
| edgeRLRTrobust | 22.22%  | 22.22%  | 16.67%  | 38.89%   | 13.84 |
| edgeRQLF       | 98.15%  | 1.85%   | /       | /        | 1.03  |
| edgeRQLFDetRate| 81.48%  | 18.52%  | /       | /        | 1.42  |
| limmatrend     | 100%    | /       | /       | /        | 0.10  |
| MASTcpm        | 81.48%  | 18.52%  | /       | /        | 1.42  |
| MASTcpmDetRate | 75.93%  | 24.07%  | /       | /        | 1.75  |
| MASTtpmDetRate | 79.63%  | 20.37%  | /       | /        | 1.38  |
| MASTtpm        | 88.89%  | 11.11%  | /       | /        | 1.28  |
| metagenomeSeq  | 40.74%  | 18.52%  | 11.11%  | 29.63%   | 9.75  |
| monocle        | 33.33%  | 46.30%  | 20.37%  | /        | 4.08  |
| monoclecensus  | 100%    | /       | /       | /        | 0.99  |
| monoclecount   | 100%    | /       | /       | /        | 0.72  |
| ROTScpm        | 22.22%  | 11.11%  | 9.26%   | 57.41%   | 13.52 |
| ROTStpm        | 22.22%  | 11.11%  | 9.26%   | 57.41%   | 14.16 |
| ROTSvoom       | 36.84%  | 15.79%  | 15.79%  | 31.58%   | 9.25  |
| SAMseq         | 37.04%  | 42.59%  | 20.37%  | /        | 4.74  |
| scDD           | 100%    | /       | /       | /        | 0.12  |
| SeuratBimod    | 100%    | /       | /       | /        | 0.18  |
| SeuratBimodIsExpr2 | 100%| /       | /       | /        | 0.19  |
| SeuratBimodnofilt| 26.67%| 13.33%  | 13.33%  | 46.67%   | 16.30 |
| SeuratTobit    | /       | 33.96%  | 24.53%  | 41.51%   | 6.76  |
| SeuratTobitnofilt| /     | 50.00%  | 8.33%   | 41.67%   | 9.05  |
| voomlimma      | 100%    | /       | /       | /        | 0.16  |
| Wilcoxon       | 100%    | /       | /       | /        | 0.83  |
| ttest          | 100%    | /       | /       | /        | 0.59  |
| BPSC           | /       | 33.33%  | 62.96%  | 3.70%    | 6.87  |
| LR             | 100%    | /       | /       | /        | 0.46  |
| NsLogit        | 100%    | /       | /       | /        | 0.64  |
| Summary        | 67.55%  | 13.46%  | 6.70    | 12.30%   | 59.70 |

- Schedule via batch submission

In practice, we do not want to submit 1728 jobs at a time by running `./run_simu_eval.sh`. Jobs are divided into batches according to their resources requirements. And jobs in a batch are further defined by two shell scripts `run_simu_eval_batchi_123.sh` and `run_simu_eval_batchi_456.sh`, each submits experiments with data characterized by simulation seed 123 and 456. 

```shell
./run_simu_eval_batch1_123.sh
```

- Notes on applying **D3E**

Unlike other methods, **D3E** is a python library. To apply **D3E** in this section, users need to create empty folders named `tmp` (for python output from D3E) and `ENV` (for creating python virtual environments inside of job submitted to Compute Canada) in the working directory. And put **D3E** code (download from https://github.com/hemberg-lab/D3E) under a folder named `D3E` in the working directory. Like other methods, DE results of **D3E** will be stored in `rds` format under `scDD_simu_result/D3E`. File `py_requirements.txt` should also be placed under the working directory for environment requirements checking. The `py_requirements.txt` looks like following for my Compute Canada environment (with python 2.7.4):

```txt
-f /cvmfs/soft.computecanada.ca/custom/python/wheelhouse/avx2
-f /cvmfs/soft.computecanada.ca/custom/python/wheelhouse/generic
numpy==1.16.3
scipy==1.2.1
mpmath==1.1.0
```

- Combine DE results of each method

The `run_simu_eval_comb.sh` will merge any existing results (might not all 54) under `scDD_simu_result`, and overwrite combined result files under `scDD_simu_result/comb_results` for methods listed in `mt` variable in `run_simu_eval_comb.sh`.

```shell
./run_simu_eval_comb.sh
```

- Prepare for plot

The `pre_scDD_simu_plot.R` will format result files under `scDD_simu_result/comb_results`, which are generated by `run_simu_eval_comb.sh`. The formatted file `scDD_simu_res.rds` will be stored under `scDD_simu_result/comb_results`. It will be used as input for `scDD_simu_plot.R`.

```R
Rscript pre_scDD_simu_plot.R "/home/xsslnc/projects/def-ubcxzh/xsslnc"
```

- Plot

When `scDD_simu_res.rds` is ready, by running `scDD_simu_plot.R`, `scDD_result1.png` and `scDD_result2.png` will be available under `/home/xsslnc/projects/def-ubcxzh/xsslnc/Figure/`.

```R
Rscript scDD_simu_plot.R "/home/xsslnc/projects/def-ubcxzh/xsslnc"
```

- Runtime plot

When `scDD_simu_res.rds` is ready, by running `scDD_time.R`, `scDD_time.png` will be available under `/home/xsslnc/projects/def-ubcxzh/xsslnc/Figure/`.

```R
Rscript scDD_time.R "/home/xsslnc/projects/def-ubcxzh/xsslnc"
```

---

### Real data study

- Split data according to cell types

Firstly, we need to subset full size GSE135893 (~20G) into smaller files, where the `counts` slot from `SCT` assay of each cell type will be stored.

```shell
./run_partition_fullsize.sh
```
The script take `GSE135893_ILD_annotated_fullsize.rds` as input. When it finishes, folder `/home/xsslnc/projects/def-ubcxzh/xsslnc/GSE135893_count/` will have data files named after cell types. 

- DE analysis on GSE135893

`run_GSE135893_DE.sh` will submit jobs to evaluate one method on 31 subtypes of GSE135893. Users can change **NsLogit** to other method’s name listed in folder `DE_methods`. Make sure that a folder named `DE_methods_count_result/` exists under `/home/xsslnc/projects/def-ubcxzh/xsslnc/GSE135893_file/`, since results will be placed there.

```shell
./run_GSE135893_DE.sh
```

- DE results analysis

When `NsLogit_count_result` is ready, we can source script `real_data_analysis.R` to summarize 31 results of **NsLogit** and published results listed in `GSE135893_file/Disease_vs_Control`. The script takes working directory as input parameter.

```R
Rscript real_data_analysis.R "/home/xsslnc/projects/def-ubcxzh/xsslnc"
```

- Pairwise DE compare on Differentiating Ciliated cells

The following script submit jobs to apply counts or cpm based methods (other than NsLogit) on Differentiating Ciliated cells of GSE135893. Make sure that folders names `DE_methods_count_result/` exists under `/home/xsslnc/projects/def-ubcxzh/xsslnc/GSE135893_file/` for those methods.

```shell
./run_GSE135893_DE_Differentiating_Ciliated.sh
```

- Pairwise DE plot

The following script will save `Differentiating_Ciliated.png` under `Figure` folder.

```R
Rscript Differentiating_Ciliated.R "/home/xsslnc/projects/def-ubcxzh/xsslnc"
```

---
> The actual script that submits job to ComputeCanada when run `run_script_name.sh` is `script_name.sh`. User need to modify `script_name.sh` for corresponding ComputeCanada configurations.

- Add tpm

To apply all DE method, tpm should be calculated. The following script take each file from `GSE135893_count` as input, query gene length from `biomaRt` and `AnnotationDbi`, calculate tpm and store updated file in `GSE135893_with_tpm`.

```shell
./run_add_tpm.sh
``` 

I manually find length for following genes. They appeared in the published results, but could not be found by `biomaRt` or `AnnotationDbi`.

```
AC090498.1       74; CTD-3252C9.4       359; SERPINA3.1    6081
CH17-360D5.2 523886; PP7080            2455; AC058791.1    3197
RP11-707G18.1   571; CH17-373J23.1      346; MATR3.1      57576
RP11-147L13.11 2918; FO538757.2        1478; RP11-80H18.3   153
CH507-42P11.8 12605; LL22NC03-N95F10.1 9936; AC013461.1  192575
1-Mar         14624; 3-Mar             4283; 9-Sep       219097
15-Sep        53789; 7-Sep           114777; RP11-295M3.4  1150
RP11-670E13.6   348;  
```

AC013461.1 from http://may2012.archive.ensembl.org/Homo\_sapiens/Gene/Summary?g=ENSG00000091436;r=2:173940163-174132738

RP11-670E13.6 from https://doi.org/10.1111/php.12858

RP11-295M3.4  from http://c-it-loci.uni-frankfurt.de/region/CIL\_0009C3


- Apply all DE methods

```shell
./run_GSE135893_DE_all_method.sh
```

- Modify filter

Since we wish to find a reasonable filter for our method, this section describe steps to calculate pair-wise density difference and KL divergence for each subtype, and refine our DE method using them instead of the published filter. 

The R script `/real_data/fckl.R` will generate new profile for each subtype using data from `GSE135893_count` folder and also store results there.

```shell
./run_fckl.sh
```

For each type, in addition to `n` genes that pass the published filter with absolute logfc > 0.25, here we also include genes that are among the top `n` in terms of pair-wise density difference `kl_diff` in our final result. And the adjusted p-value is re-calculated. The following command will print the analysis table with an additional column of the percentage of `n` in total number of genes in corresponding subtypes. 

**NOTE** For each gene, denote `p1` and `p2` as the probability density of count for each group, and `kld(px, py) = sum(px*log px/log py)`. The column `kl_max` is `max(kld(p1, p2), kld(p2, p1))`, column `kl_mean` is `mean(kld(p1, p2), kld(p2, p1))`, and `kl_diff = sum(abs(p1-p2))`.

```R
Rscript real_data_analysis_new_filter.R "/home/xsslnc/projects/def-ubcxzh/xsslnc"
```

or

```shell
./run_real_data_analysis_new_filter.sh
```

### Confirm DE genes (from GSE135893) using GSE136831

- Split data according to cell types

Similarly, raw data files of GSE136831 which are placed under `GSE136831_file` are split according to cell types (Manuscript Identity) and sub-files are stored under `GSE136831_count`.

```shell
./run_partition_fullsize_GSE136831.sh
```

- DE analysis on GSE136831

Then we apply the method on each sub-files and store the results under `GSE136831_file/NsLogit_count_result`.

```shell
./run_GSE136831_DE.sh
```
 

