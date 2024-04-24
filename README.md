The impact of sex on DNA methylation across nine cancer types
=
### Abstract
Sex differences among cancer tissues have been investigated, yet focusing solely on intra-tumoral comparisons may overlook the complexities of sex effects on malignancies. Here, we investigated the influence of sex on DNA methylation (DNAm) across nine types of non-reproductive cancers and their matched normal adjacent tissues (NATs) from The Cancer Genome Atlas. We observed a striking reduction of sex differential methylation in cancers compared to NATs. Comparing the cancer-related aberrant DNAm between sexes revealed differences in the magnitudes of DNAm changes, termed “amplification”, contributing to the sex differences observed in cancers. We then prioritized 3,361 female-amplifiers and 11,837 male-amplifiers across cancers, representing the observed sex-related magnitude differences. These amplifiers were enriched in binding sites of various transcription factors, including P53, SOX2, and CTCF. Integrating sex-amplifiers with downstream gene expression highlighted consistent amplified differential expression, indicating sex-dependent regulation in cancers. Together, these findings reveal an extensive attenuation of sex differences in DNAm across cancers and emphasize that amplification effects can account for sex differences in malignancies. 

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

Jiaqi Zhou, Miao Li, Yu Chen, Shangzi Wang, Danke Wang, Chen Suo, and Xingdong Chen*. The impact of sex on DNA methylation across nine cancer types. 2024<br />

*Corresponding author<br />

DOI: xxx
