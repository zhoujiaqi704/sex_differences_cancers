The impact of sex on DNA methylation across nine cancer types
=
### Abstract
Sex differences exist in numerous non-reproductive cancers, but whether the extent of these differences is reorganized during tumorigenesis remains uncharacterized. Here, we investigated the influence of sex on DNA methylation (DNAm) across nine types of cancers from The Cancer Genome Atlas. Comparative analysis with matched normal adjacent tissues revealed a predominant attenuation of sex differences in DNAm across all nine cancers. Comparing the correlation and magnitude relationships of cancer-related DNAm changes between sexes indicated that the DNAm changes were highly correlated, but the magnitude differed between sexes (termed “amplification”). Among cancer-related CpGs, we prioritized 3361 female-amplifiers and 11837 male-amplifiers across cancers. These sex-amplifiers were sex-specifically enriched in a wide range of transcription factors (e.g., SMAD2/4, SOX2, and CTCF). Integrating sex-amplifiers with the downstream transcriptome highlighted consistent sex-specific differential expression (e.g., GATA2 and ZBP1). These genes regulated by sex-amplifiers were enriched for cancer development-related pathways, including vasculature development and sprouting angiogenesis. Taken together, our findings reveal an extensive attenuation of sex differences in DNAm across cancers and emphasize the amplification effects as the primary driver of sex-related influences on malignancies.  

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

- *Sex-stratified eQTM screeaning: 08_MethyExpression - eQTM.R*<br />
- *Integration of sex-amplifiers with downstream transcriptome using eQTM pairs: 08_MethyExpression - IntegrationAnalyses.R*<br />

## Data availability

All data, including the main text and supplementary materials, are available at [insert location]

## Citation

If you use anything in this repository please cite: xxx
