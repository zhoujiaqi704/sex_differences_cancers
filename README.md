The impact of sex on DNA methylation across nine cancer types
=
### Abstract
DNA methylation (DNAm) plays a critical role in both sex differences and cancer development, but the mechanisms underlying their interplay are incompletely understood. To decipher the effects of sex on DNAm in cancers, we characterized the sex differences within and between nine types of non-reproductive cancers and their paired normal adjacent tissues (NATs) from The Cancer Genome Atlas. Using paired samples, we observed a striking reduction of sex differences in DNAm within cancers compared to NATs. We identified 3,452 CpGs (Bonferroni-adjusted P < 0.05) associated with such reduction and showed that the linked genes were involved in X chromosome inactivation. Analyses of differential methylation between cancers and NATs in a sex-stratified design revealed significant variation in the number of differentially methylated positions (DMPs) across sexes. By quantifying heterogeneity in effect sizes among sexes spanning a wide range of correlations, we found that although aberrant DNAm events were extensively correlated, effect sizes varied greatly among sexes within cancers, which we termed “amplification”. Based on the weights of sex-biased amplification effects, we categorized these cancer types into female- and male-biased groups. We prioritized a series of cancer-related DMPs that exhibit these observed sex amplification events. These sex-amplified DMPs were enriched in binding sites of various transcription factors, including P53, SOX2, and CTCF. Integrating sex-amplified DMPs with downstream gene expression through eQTMs pinpointed consistent sex-biased differential expression in cancers. Overall, our findings delineate the epigenetic basis of sex differences in cancers, with implications for the development of sex-specific biomarkers and precision medicine.

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

All data, including the main text and supplementary materials, are available at xxx.<br />
The supplementary tables are available at [Zenodo](https://zenodo.org/records/11238927).<br />

## Citation

If you use anything in this repository please cite:<br />

Jiaqi Zhou, Miao Li, Yu Chen, Shangzi Wang, Danke Wang, Chen Suo, and Xingdong Chen*. The impact of sex on DNA methylation across nine cancer types. 2024<br />

*Corresponding author<br />

DOI: xxx
