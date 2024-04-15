The impact of sex on DNA methylation across nine cancer types
=
### Abstract
Sex differences in cancer tissue genomics have been explored, yet focusing solely on intra-tumoral comparisons may overlook the complexities of sex effects on malignancies. Here, we investigated the influence of sex on DNA methylation (DNAm) across nine non-reproductive cancers and their matched normal adjacent tissues (NATs) from The Cancer Genome Atlas. We observed a striking reduction of sex differential methylation in cancers compared to NATs. Comparing the correlation and magnitude of cancer-related DNAm changes between sexes revealed highly correlated DNAm changes, but with differing magnitudes (termed “amplification”). Among cancer-related CpGs, we prioritized 3,361 female-amplifiers and 11,837 male-amplifiers across cancers, which were enriched in various transcription factor binding sites (e.g., SMAD2/4 and SOX2). Integrating sex-amplifiers with the downstream transcriptome highlighted consistent sex-specific differential expression (e.g., GATA2 and ZBP1). Genes regulated by sex-amplifiers were implicated in cancer development-related pathways, including vasculature development and sprouting angiogenesis. These findings reveal an extensive attenuation of sex differences in DNAm across cancers and emphasize the amplification effects as the primary driver of sex-related influences on malignancies. 

## Script summary

### Differential methylation analyses

- *Preprocessing and quality control (discovery): 01_A_Preprocessing_DNAm*<br />
- *Sex-related differential methylation: 02_Attenuation_DNAm - Attenuation_DNAm.R*<br />
- *AMP identification and annotation: 02_Attenuation_DNAm - Attenuation_DNAm.R and AMP_annotation_enrichment.R*<br />
- *Sex-stratified differential methylation: 03_Differential_Interaction_DNAm*<br />
- *Sex-by-cancer interaction analyses: 03_Differential_Interaction_DNAm*<br />
- *Permutation analyses & Downsampling analyses & Bootstrap test: 04_StatisticalEvaluation*<br />
- *Slope calculation: 04_StatisticalEvaluation - slope_significance.R*<br />
- *Replication for sex-stratified results: 01_B_Replication_DNAm*<br />

### Amplification effects estimation

- *Mash model: 05_AmplificationEffect_mash - mash_hnsc_nominalCpG.R and mash_hnsc_allCpG.R*<br />
- *Sex-biased amplification groups: 05_AmplificationEffect_mash - sexBias_by_amplification.R*<br />
- *Limma vs. Mash: 05_AmplificationEffect_mash - limma_vs_mash.R*<br />

### Sex amplifiers identification

- *Sex-heterogeneity effects: 06_SexAmplifier - sex-het.R*<br />
- *Sex-amplifiers and sex-shared effects: 06_SexAmplifier - sex-amplifier.R*<br />
- *Transcription factor binding sites enrichment: 06_SexAmplifier - TFBSenrich_preparation.R and EnrichmentAnalyses-sex-amplification-effects.R*<br />

### RNA-seq

- *Preprocessing and quality control: 07_RNA-seq.R*<br />
- *Sex-stratified differential expression: 07_RNA-seq - sex_stratified_DEG.R*<br />

### Methylome and Transcriptome correlation (eQTM)

- *Sex-stratified eQTM screening: 08_MethyExpression - eQTM.R*<br />
- *Integration of sex-amplifiers with downstream transcriptome using eQTM pairs: 08_MethyExpression - IntegrationAnalyses.R*<br />

## Data availability

All data, including the main text and supplementary materials, are available at xxx

## Citation

If you use anything in this repository please cite:<br />

Jiaqi Zhou, Miao Li, Yu Chen, Shangzi Wang, Danke Wang, and Xingdong Chen*. The impact of sex on DNA methylation across nine cancer types. 2024<br />

*Corresponding author<br />

DOI: xxx
