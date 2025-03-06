# **Single-cell RNA atlas of endometrial tissue to study endometriosis**

## **Introduction**

This repository is intended to support the analysis presented in the accompanying publication "Tracing Endometriosis: Coupling deeply phenotyped, single-cell based Endometrial Differences and AI for disease pathology and prediction".
It contains custom code for atlas integration, Sample-wise Principal Component Analysis, annotation transfer, time trajectory analysis, cluster-wise differential expression analysis, and ligand-receptor analysis, as well as the code to recreate figures and extended data figures.

## **Associated Publication**
### **Publication Title**

Tracing Endometriosis: Coupling deeply phenotyped, single-cell based Endometrial Differences and AI for disease pathology and prediction

### **Publication Abstract**

Endometriosis, affecting 1 in 9 women, presents treatment and diagnostic challenges. To address these issues, we generated the biggest single-cell atlas of endometrial tissue to date, comprising 466,371 cells from 35 endometriosis and 25 non-endometriosis patients without exogenous hormonal treatment. Detailed analysis reveals significant gene expression changes and altered receptor-ligand interactions present in the endometrium of endometriosis patients, including increased inflammation, adhesion, proliferation, cell survival, and angiogenesis in various cell types. These alterations may enhance endometriosis lesion formation and offer novel therapeutic targets. Using ScaiVision, we trained neural network models to predict endometriosis of varying disease severity (median AUC = 0.83) including one model based solely on a set of  11 genes confirmed as dysregulated in endometriosis patients through differential expression analysis. Our models, while not yet externally validated, can serve as a tool for hypothesis generation as well as a starting point for further clinical development. In conclusion, our findings illuminate numerous pathway and ligand-receptor changes in the endometrium of endometriosis patients, offering insights into pathophysiology, targets for novel treatments, and predictive models for enhanced outcomes in endometriosis management.

### **Link to Publication**
Will be integrated upon publication

## **Data availability**

All raw sequencing data and the processed Seurat object are available at NCBI’s Gene Expression Omnibus (series accession number: GSE111976; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE266265).


## **Installation**

To install Project Title, follow these steps:

1. Clone the repository: 
    ```{bash}
    cd {path to save code}
    git clone https://github.com/DuempelmannLea/endometriosis_endometrium_scRNA_atlas.git
    ```
2. All processed data is stored in the GSE111976. To reproduce the analysis, download the data and link it to the _Data folder within the endometriosis_endometrium_scRNA_atlas code repository:
    ```{bash}
    cd endometriosis_endometrium_scRNA_atlas/_Data
    ln -s /path/to/EndoAtlas.rds EndoAtlas.rds
    ln -s /path/to/64_1r_L1_R1_001_GiiMsim2jmMj.fastq.gz 64_1r_L1_R1_001_GiiMsim2jmMj.fastq.gz
    ...
    ```

## **Usage**

Follow these steps:

1. Open the project in your favorite code editor.
2. Modify the source code to fit your needs.
3. Use the project as desired.

## **License**

This repository is intended to support the analysis presented in the accompanying publication. Portions of the code may include adaptations from vignettes or external packages. 
While the code is available for reuse, it is primarily provided as a resource to complement the publication and does not include a formal license.

## **Authors and Acknowledgement**
### **Authors**
**[Lea Duempelmann](https://github.com/DuempelmannLea)**,
**[Shaoline Sheppard](https://github.com/Shaoline)**,
**[Angelo Duo](https://github.com/duoa)**,
Dennis Goehlsdorf,
Sukalp Muzumdar,
Ryan Lusby,
Sarah Carl

### **Acknowledgment**
We would like to acknowledge the additional authors of the accompanying publication for their invaluable contributions to this project:
Jitka Skrabalova, Brett McKinnon, Thomas Andrieu, Wiebke Solass, Cinzia Donato, Peter Nestorov, and Michael D. Mueller.

## **Contact**

If you have any questions, feedback, or issues, please open an issue on this repository, tagging @DuempelmannLea. 
