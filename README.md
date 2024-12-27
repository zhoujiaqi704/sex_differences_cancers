Attenuated sex-related DNA methylation differences in cancer highlight the magnitude bias mediating existing disparities
=

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

Zhou, J., Li, M., Chen, Y. et al. Attenuated sex-related DNA methylation differences in cancer highlight the magnitude bias mediating existing disparities. Biol Sex Differ 15, 106 (2024). https://doi.org/10.1186/s13293-024-00682-4<br />
