# Scripts for running diff. accessible chromatin regions (DACRs) analysis.

* ```Calculate_subtype_DACRs.R``` -- subtype specific DACRs for HTAN BRCA project, using cohort merged object.

* ```Run_DACRs_cancer_tissue_specific.FindMaarkers.R``` -- tissue and cancer-cell specific DACRs (pan-cancer snATAC paper). Compare cancer cells from each cohort vs cancer cells from all other cohorts using Tumor-Normal merged pan-cancer object.

* ```Tumor_normal_DACRs.CohortObj.R``` -- cancer-cell specific DACRs (pan-cancer snATAC paper). Compare cancer cells vs cell-of-origin cell type for each cohort, using cohort level object, or 2-cohort merged objs for HNSCC and OV.