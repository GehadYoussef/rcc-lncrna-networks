# Non-exonic read content reorganises lncRNA networks and masks tumour cell programmes in renal cancer

Code and results for the article by Gehad Youssef and Namshik Han (Milner
Therapeutics Institute, University of Cambridge). Correspondence: Namshik Han,
nh417@cam.ac.uk.

Scripts in `src/` write every number, table and figure in the article to
`results/` or `figures/`.

## Repository layout

```
src/
  R/                 bulk RNA-seq pipeline (numbered stages, 00_config.R, run_all.R)
  singlecell/        single-cell pipeline (numbered stages, lib_sc.R, run_all.R)
  config/            single-cell pipeline settings (singlecell.yml)
data/
  manifests/         single-cell download manifest and sample sheets
  raw/               downloads (GDC, GEO, Mendeley), created on first run
  derived/           caches and intermediate objects, created on first run
results/             result tables (tab-separated) and lock files
  singlecell/        single-cell results and candidate locks
figures/             main and supplementary figures (PNG and SVG)
manuscript/          article, Supplementary Information and figure files with legends (PDF)
```

## Requirements

Tested with R 4.4.2 on Windows. A full run downloads about 60 GB. `data/raw/` and `data/derived/` are not tracked by git.

Principal packages, with the versions used: WGCNA 1.73, survival 3.7-0, glmnet 4.1-10, data.table 1.16.4, clusterProfiler 4.12.6, estimate 1.0.13, edgeR 4.2.2, DESeq2 1.44.0, limma 3.60.6, GenomicRanges 1.56.2, metafor 4.8.0 and ruv 0.9.7.2.

Other packages:
- **Bioconductor:** org.Hs.eg.db, GO.db, AnnotationDbi, IRanges, S4Vectors, BiocParallel, SingleCellExperiment, SingleR, scuttle and rhdf5.
- **CRAN:** Matrix, matrixStats, car, timeROC, survminer, energy, SeuratObject, yaml, digest, jsonlite, httr, curl, xml2 and readxl.
- **Plotting:** ggplot2, ggrepel, patchwork, pheatmap, scales, svglite, ragg and systemfonts.
- **Parallel computing:** foreach and doParallel.

## Running the pipeline

From the repository root:

```bash
Rscript src/R/run_all.R
```

`run_all.R` runs the stages in order. `00_download_gdc.R` runs first and fetches or md5-verifies the TCGA-KIRC files. The single-cell pipeline (`src/singlecell/run_all.R`) runs between stages 35 and 36, and the tables and figures are written last.

A stage whose cache file exists in `data/derived/cache/` is skipped. To re-run it, delete that file (the header of `run_all.R` lists them). Thresholds, filters, seeds and model parameters are set in `src/R/00_config.R`, and single-cell settings in `src/config/singlecell.yml`.

Not run by default:
- `02b_lncrna_tuning.R`, the outcome-blind sweep that chose the lncRNA network parameters now fixed in `00_config.R`. Run it by hand.
- `04_dcor.R` and `04b_dcor_diagnostics.R`, a distance-correlation screen not used in the article. Set `RUN_DCOR` in `run_all.R` to include them.

To run the single-cell pipeline on its own:

```bash
Rscript src/singlecell/run_all.R [--no-download] [--force] [--new-lock]
```

## Data

No new data were generated. All inputs are open access, and the pipeline downloads them.

- **Bulk RNA-seq.** STAR-Counts files and harmonised clinical fields from the NCI Genomic Data Commons (https://portal.gdc.cancer.gov), Data Release 46.0, for TCGA-KIRC, TCGA-KIRP, TCGA-KICH, CPTAC-3 and all 33 TCGA projects. TCGA-KIRC grade comes from the Biospecimen Core Resource clinical XML supplements.
- **Outcomes.** Disease-specific survival, progression-free interval and pan-cancer overall survival from the TCGA Clinical Data Resource (Liu et al., Cell 2018).
- **Single-cell RNA-seq.** GEO series GSE159115, GSE222703 and GSE207493, and the dataset of Li et al. (Cancer Cell 2022) on Mendeley Data (g67bkbnhhg). `data/manifests/` lists every file and its source.
- **Annotation.** GENCODE v36.

`results/SupplementaryData1_case_list.tsv` and `results/38_pancancer_file_manifest.tsv` list every expression file used, with its GDC identifier, md5 checksum and retention status.

## Bulk pipeline stages (`src/R`)

| Script | Purpose |
|---|---|
| `00_config.R` | paths, seed, every filter and model parameter, shared helpers |
| `00_download_gdc.R` | TCGA-KIRC STAR-Counts files and clinical XML from the GDC, md5-verified |
| `01_build_data.R` | discovery cohort, sample inclusion, CONSORT counts |
| `02_wgcna.R` | adjustment for the alignment metrics, protein-coding and lncRNA co-expression networks, eigengenes |
| `02b_lncrna_tuning.R` | outcome-blind sweep of lncRNA network parameters (run by hand) |
| `03_survival.R` | module Cox models under the three specifications, proportional-hazards tests |
| `04_dcor.R`, `04b_dcor_diagnostics.R` | distance-correlation screen (not run by default) |
| `05_ml.R` | penalised Cox models with repeated cross-validation |
| `06_enrichment.R` | Gene Ontology over-representation of protein-coding modules |
| `07_purity_qc.R` | ESTIMATE scores, STAR alignment metrics, the leading lncRNA axis |
| `08_validation_data.R` | CPTAC-3 cohort through the same pipeline |
| `09_validate.R` | locked model and association replication in CPTAC-3 |
| `10_paper_analyses.R` | validation bootstrap, guilt-by-association enrichment, module composition |
| `11_subtype_specificity.R` | projection onto TCGA-KIRP and TCGA-KICH |
| `12_model_diagnostics.R` | calibration, Brier scores, decision curves, fold-wise module rebuild |
| `13_figures.R` | main figures |
| `14_tables.R` | main tables |
| `15_sensitivity_tissue_source_site.R` | tissue source site as a batch variable |
| `18_signsplit_enrichment.R` | sign-split guilt-by-association enrichment |
| `19_eigengene_correlation.R` | coupling of lncRNA eigengenes to the catabolic protein-coding module |
| `20_turquoise_units_and_locked_model.R` | drop-one-module test and unit-matched replication |
| `22_locked_model_export.R` | the locked prediction model written out in full |
| `23_technical_axis_extended.R` | RNA integrity, biospecimen and batch variables, the unadjusted lncRNA network |
| `24_module_preservation.R` | WGCNA module preservation in the external cohorts |
| `25_matched_normal_control.R` | matched normal kidney as a tumour-free control |
| `26_published_signatures.R` | published TCGA-KIRC lncRNA signatures re-scored |
| `27_metric_biology.R` | correlates of the non-feature read fraction |
| `28_mutations_and_lncRNA_classes.R` | driver mutations and positional lncRNA classes |
| `29_normalisation_check.R` | the axis under five normalisations |
| `30_endpoint_sensitivity.R` | disease-specific and progression-free endpoints |
| `31_supplementary_data_and_index.R` | Supplementary Data 1 and package versions |
| `32_supplementary_figures.R` | supplementary figures |
| `33_sex_stratified.R` | sex-disaggregated module effects |
| `34_proximal_tubule_reading.R` | proximal-tubule marker score and ccA/ccB |
| `35_axis_projection_subtypes.R` | discovery axis projected onto TCGA-KIRP and TCGA-KICH |
| `36_singlecell_module_localisation.R` | single-cell localisation of the 29 modules |
| `37_network_rewiring.R` | adjusted against unadjusted lncRNA network partitions |
| `38_pancancer_build.R` | 33 TCGA projects from the GDC under the discovery rules |
| `39_pancancer_axis.R` | lncRNA-specific axis rule and positional-class gradient per project |
| `40_pancancer_network.R` | network reorganisation rule per project |
| `41_pancancer_metric_survival.R` | hazard of the metric per project and random-effects meta-analysis |
| `42_singlecell_truth_correction.R` | single-cell ground truth for the correction |
| `43_pancancer_dose_response.R` | within-cohort spread of the metric against the rule quantities (exploratory) |
| `44_ruv3_prps_comparison.R` | comparison with RUV-III using pseudo-replicates of pseudo-samples |
| `pan_helpers.R` | shared pan-cancer functions, including the decision-rule check |
| `ora_hypergeometric.R` | hypergeometric over-representation test |
| `run_all.R` | runs every stage in order |

## Single-cell pipeline stages (`src/singlecell`)

| Script | Purpose |
|---|---|
| `00_config_sc.R`, `lib_sc.R` | paths, settings from `singlecell.yml`, shared functions |
| `01_validate_manifest.R` | validates the manifests and downloads inputs with checksums |
| `02_import_dataset.R` | imports each dataset and removes outcome fields |
| `03_harmonize_metadata.R` | harmonises metadata and cell compartment labels |
| `04_gene_annotation.R` | GENCODE v36 gene annotation and lncRNA coverage |
| `05_cell_qc_and_annotations.R` | cell quality control and malignant-cell labels |
| `12_render_report.R` | data audit and the quality gate that stage 06 requires |
| `06_build_pseudobulk.R` | donor-level pseudobulk |
| `07_within_dataset_de.R` | within-dataset differential expression (edgeR quasi-likelihood) |
| `08_meta_analysis.R` | random-effects meta-analysis across datasets and candidate tiers |
| `09_candidate_lock.R` | writes the candidate lock, which cannot be changed without `--new-lock` |
| `10_bulk_handoff.R` | availability and technical sensitivity of candidates in bulk cohorts |
| `11_bulk_prognostic.R` | prognostic model locked before CPTAC-3 outcomes are read |
| `run_all.R` | runs the stages in the order above |

## Where each figure and table comes from

| Item | Source tables in `results/` |
|---|---|
| Fig. 1 | `01_consort_flow.tsv`, `01_table1_cohort.tsv`, `08_validation_table1.tsv`, `11_subtype_cohort_summary.tsv`, `11_cohort_quality_metrics.tsv` |
| Fig. 2 | `27_axis_scores_per_sample.tsv`, `07_per_gene_quality_correlation.tsv`, `07_per_gene_quality_summary.tsv`, `23_axis_nested_extended.tsv`, `23_batch_variance.tsv`, `23_noFeature_determinants.tsv`, `27_metric_biology_correlations.tsv`, `25_normal_axis_projection.tsv`, `35_axis_projection_subtypes.tsv`, `07_axis_vs_purity_qc_correlations.tsv`, `09_axis_replication_cptac.tsv`, `11_axis_replication_subtypes.tsv` |
| Fig. 3 | `39_pancancer_axis.tsv`, `39_pancancer_class_gradient.tsv`, `40_pancancer_network.tsv`, `41_pancancer_metric_survival.tsv`, `41_pancancer_metric_meta.tsv`, `43_pancancer_dose_response_projects.tsv`, `43_pancancer_dose_response.tsv` |
| Fig. 4 | `12_technical_adjustment_sensitivity.tsv`, `02_lncRNA_module_trait.tsv`, `03_mRNA_module_survival.tsv`, `03_lncRNA_module_survival.tsv` |
| Fig. 5 | `09_module_replication.tsv`, `09_module_replication_parsimonious.tsv`, `11_subtype_module_effects.tsv`, `11_subtype_heterogeneity.tsv`, `24_module_preservation_lncRNA.tsv`, `24_module_preservation_mRNA.tsv` |
| Fig. 6 | `06_GO_enrichment_top10_per_module.tsv`, `03_mRNA_module_survival.tsv`, `03_lncRNA_module_survival.tsv`, `10_lncRNA_guilt_by_association_GO_signed.tsv`, `11_lnc_black_coupling_by_cohort.tsv`, `36_sc_module_localisation_per_dataset.tsv` |
| Fig. 7 | `09_delta_cindex_validation.tsv`, `05_delta_cindex.tsv`, `12_foldwise_per_assignment.tsv`, `09_validation_cindex.tsv`, `12_calibration_bins.tsv`, `12_calibration_summary.tsv`, `12_validation_decision_curve.tsv`, `12_brier_scores.tsv` |
| Table 1 | `01_table1_cohort.tsv`, `08_validation_table1.tsv`, `11_cohort_quality_metrics.tsv`, `11_subtype_cohort_summary.tsv` |
| Table 2 | `03_mRNA_module_survival.tsv`, `03_lncRNA_module_survival.tsv`, `09_module_replication.tsv`, `09_module_replication_parsimonious.tsv` |
| Table 3 | `09_validation_cindex.tsv`, `09_delta_cindex_validation.tsv`, `05_delta_cindex.tsv`, `12_foldwise_module_sensitivity.tsv`, `12_calibration_summary.tsv`, `12_brier_scores.tsv` |
| Supplementary Fig. 1 | `02_lncRNA_soft_threshold.tsv`, `02_mRNA_soft_threshold.tsv`, `02_network_summary.tsv`, `02_lncRNA_module_sizes.tsv`, `02_mRNA_module_sizes.tsv` |
| Supplementary Fig. 2 | `29_normalisation_pc1.tsv`, `29_normalisation_module_contrast.tsv` |
| Supplementary Fig. 3 | `25_tumour_vs_normal_distribution.tsv`, `25_paired_scores_per_sample.tsv`, `25_normal_axis_projection.tsv`, `25_normal_metric_survival.tsv` |
| Supplementary Fig. 4 | `23_biospecimen_kirc.tsv`, `23_axis_vs_biospecimen.tsv`, `27_axis_scores_per_sample.tsv` |
| Supplementary Fig. 5 | `27_metric_biology_correlations.tsv`, `28_noFeature_vs_mutations.tsv` |
| Supplementary Fig. 6 | `42_sc_truth_concordance.tsv`, `42_sc_truth_score_coherence.tsv`, `42_sc_truth_technical_candidates.tsv` |
| Supplementary Fig. 7 | `26_random_signature_null.tsv`, `26_published_signatures_models.tsv`, `26_published_signatures_correlations.tsv`, `26_headline_ranges_by_provenance.tsv` |
| Supplementary Fig. 8 | `24_module_preservation_lncRNA.tsv`, `24_module_preservation_mRNA.tsv`, `24_module_preservation_summary.tsv`, `02_lncRNA_module_sizes.tsv`, `02_mRNA_module_sizes.tsv` |
| Supplementary Fig. 9 | `10_lncRNA_guilt_by_association_GO_signed.tsv`, `10_lncRNA_module_composition.tsv`, `28_module_composition_by_class.tsv` |
| Supplementary Fig. 10 | `36_sc_module_localisation_combined.tsv` |
| Supplementary Fig. 11 | `05_cv_cindex_per_repeat.tsv`, `05_hub_module_selection_frequency.tsv` |
| Supplementary Fig. 12 | `12_calibration_bins.tsv`, `12_calibration_summary.tsv`, `12_validation_decision_curve.tsv` |
| Supplementary Fig. 13 | `30_endpoint_modules.tsv`, `30_endpoint_exposure.tsv`, `30_endpoint_cv_increment.tsv`, `30_endpoint_foldwise.tsv` |
| Supplementary Table 1 | `28_mutation_coverage.tsv`, `28_mutation_frequency.tsv`, `28_noFeature_vs_mutations.tsv`, `28_mutation_survival_models.tsv`, `28_lncRNA_positional_classes.tsv`, `28_quality_correlation_by_class.tsv`, `28_class_comparison_tests.tsv`, `28_module_composition_by_class.tsv` |
| Supplementary Table 2 | `29_normalisation_pc1.tsv`, `29_size_factors_vs_quality.tsv`, `29_normalisation_agreement.tsv`, `29_normalisation_module_summary.tsv`, `29_normalisation_module_contrast.tsv`, `29_bootstrap_seed_stability.tsv` |
| Supplementary Table 3 | `25_normal_consort.tsv`, `25_subtype_normal_availability.tsv`, `25_normal_axis_projection.tsv`, `25_tumour_vs_normal_distribution.tsv`, `25_paired_tumour_normal.tsv`, `25_paired_variance_partition.tsv`, `25_normal_metric_survival.tsv` |
| Supplementary Table 4 | `25_normal_axis_projection.tsv`, `35_axis_projection_subtypes.tsv` |
| Supplementary Table 5 | `23_unadjusted_lncRNA_network_summary.tsv`, `37_network_partition_agreement.tsv`, `37_adjusted_module_fate.tsv`, `37_dominant_module_positional_class.tsv` |
| Supplementary Table 6 | `38_pancancer_cohorts.tsv`, `39_pancancer_axis.tsv`, `39_pancancer_class_gradient.tsv`, `40_pancancer_network.tsv`, `41_pancancer_metric_survival.tsv`, `41_pancancer_metric_meta.tsv`, `43_pancancer_dose_response.tsv`, `43_pancancer_dose_response_projects.tsv` |
| Supplementary Table 7 | `23_axis_nested_extended.tsv`, `23_noFeature_exposure_extended.tsv`, `07_noFeature_as_exposure.tsv`, `23_axis_vs_biospecimen.tsv`, `23_batch_variance.tsv`, `23_noFeature_determinants.tsv` |
| Supplementary Table 8 | `30_cdr_coverage_and_os_concordance.tsv`, `30_endpoint_modules.tsv`, `30_endpoint_exposure.tsv`, `30_endpoint_axis_nested.tsv`, `30_endpoint_cv_increment.tsv`, `30_endpoint_locked_model.tsv`, `30_endpoint_foldwise.tsv` |
| Supplementary Table 9 | `27_metric_biology_correlations.tsv`, `27_module_gene_level_metric_correlation.tsv`, `06_GO_enrichment_top10_per_module.tsv`, `27_noFeature_variance_partition.tsv`, `27_metric_attenuation.tsv`, `27_collinearity_diagnostics.tsv`, `27_noFeature_reverse_check.tsv`, `07_noFeature_vs_clinical.tsv` |
| Supplementary Table 10 | `42_sc_truth_concordance.tsv`, `42_sc_truth_score_coherence.tsv`, `42_sc_truth_technical_candidates.tsv`, `42_sc_truth_gene_level.tsv` |
| Supplementary Table 11 | `44_prps_design.tsv`, `44_prps_axis.tsv`, `44_prps_network.tsv`, `44_prps_sc_truth.tsv`, `44_prps_module_hr.tsv` |
| Supplementary Table 12 | `12_technical_adjustment_sensitivity.tsv`, `12_technical_adjustment_arms.tsv`, `23_technical_adjustment_with_rin.tsv` |
| Supplementary Table 13 | `26_signature_provenance_key.tsv`, `26_published_signatures_genes.tsv`, `26_published_signatures_correlations.tsv`, `26_published_signatures_models.tsv`, `26_headline_ranges_by_provenance.tsv`, `26_random_signature_null.tsv` |
| Supplementary Table 14 | `03_nodal_metastasis_sensitivity.tsv` |
| Supplementary Table 15 | `33_sex_correlations.tsv`, `33_noFeature_by_sex.tsv`, `33_sex_interaction_models.tsv`, `33_sex_stratified_hr.tsv`, `33_sex_contrast_reconciliation.tsv`, `33_xy_lncRNA_counts_by_module.tsv`, `33_sensitivity_analyses.tsv`, `33_sex_label_check.tsv` |
| Supplementary Table 16 | `11_subtype_module_effects.tsv`, `24_module_preservation_lncRNA.tsv`, `24_module_preservation_mRNA.tsv` |
| Supplementary Table 17 | `24_module_preservation_summary.tsv`, `24_module_preservation_lncRNA.tsv`, `24_module_preservation_mRNA.tsv` |
| Supplementary Table 18 | `34_pt_marker_set.tsv`, `34_green_hubs_pt_flag.tsv`, `34_pt_score_correlations.tsv`, `34_pt_adjusted_models.tsv`, `34_ccAB_class_summary.tsv`, `34_green_estimate_overlap.tsv`, `34_normal_admixture.tsv` |
| Supplementary Table 19 | `36_sc_module_gene_coverage.tsv`, `36_sc_module_localisation_per_dataset.tsv`, `36_sc_module_localisation_combined.tsv`, `36_sc_module_lncRNA_pooled_meta_summary.tsv` |
| Supplementary Table 20 | `22_locked_model_coefficients.tsv`, `05_final_model_coefficients.tsv`, `22_module_score_constants.tsv`, `22_locked_model_standardisation.tsv`, `22_locked_model_baseline.tsv` |
| Supplementary Table 21 | `09_validation_cindex.tsv`, `12_calibration_summary.tsv`, `12_brier_scores.tsv`, `12_validation_decision_curve_summary.tsv`, `20_bootstrap_delta_cindex.tsv`, `20_score_spread_by_standardisation.tsv`, `20_replication_unitmatched.tsv` |
| Supplementary Data 1 | `SupplementaryData1_case_list.tsv`, `SupplementaryData1_column_definitions.tsv`, `31_supplementary_data1_summary.tsv` |
| Supplementary Data 2 | `02_lncRNA_module_genes.tsv`, `02_mRNA_module_genes.tsv`, `22_module_loadings.tsv` |
| Supplementary Data 3 | `03_all_modules_principal.tsv` |
| Supplementary Data 4 | `06_GO_enrichment_all_modules.tsv` |
| Supplementary Data 5 | `38_pancancer_file_manifest.tsv` |

## Locked decisions

- **Pan-cancer decision rules.** Set in `src/R/00_config.R`. `results/38_decision_rules_lock.json` records their values, when they were fixed and the development-repository commit at that time. Stages 39 to 41 stop if a rule differs from the lock.
- **Single-cell candidates.** `results/singlecell/candidate_lock.json` and the bulk locks record checksums of the locked tables. Stages 10 and 11 stop if a locked table changes.

## Licence

The code in `src/` is released under the MIT licence (see `LICENSE`). The article, figures and result tables remain the copyright of the authors.
