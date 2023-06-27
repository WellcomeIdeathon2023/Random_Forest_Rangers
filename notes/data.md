.
├── common_files (these are lookup tables explaining the values)
│   ├── lk_adverse_event_severity.csv
│   ├── lk_age_event.csv
│   ├── lk_analyte.csv
│   ├── lk_ancestral_population.csv
│   ├── lk_arm_type.csv
│   ├── lk_cell_population.csv
│   ├── lk_cell_population_marker.csv
│   ├── lk_compound_role.csv
│   ├── lk_criterion_category.csv
│   ├── lk_data_completeness.csv
│   ├── lk_disease.csv
│   ├── lk_disease_stage.csv
│   ├── lk_ethnicity.csv
│   ├── lk_exp_measurement_tech.csv
│   ├── lk_exposure_material.csv
│   ├── lk_exposure_process.csv
│   ├── lk_expsample_result_schema.csv
│   ├── lk_file_detail.csv
│   ├── lk_gender.csv
│   ├── lk_hmdb.csv
│   ├── lk_lab_test_name.csv
│   ├── lk_lab_test_panel_name.csv
│   ├── lk_locus_name.csv
│   ├── lk_mass_spectrometry_type.csv
│   ├── lk_organization.csv
│   ├── lk_personnel_role.csv
│   ├── lk_plate_type.csv
│   ├── lk_protein_name.csv
│   ├── lk_protocol_type.csv
│   ├── lk_public_repository.csv
│   ├── lk_race.csv
│   ├── lk_reagent_type.csv
│   ├── lk_research_focus.csv
│   ├── lk_sample_type.csv
│   ├── lk_source_type.csv
│   ├── lk_species.csv
│   ├── lk_study_file_type.csv
│   ├── lk_study_panel.csv
│   ├── lk_subject_location.csv
│   ├── lk_t0_event.csv
│   ├── lk_time_unit.csv
│   ├── lk_transcript_type.csv
│   ├── lk_unit_of_measure.csv
│   ├── lk_user_role_type.csv
│   ├── lk_virus_strain.csv
│   └── lk_visibility_category.csv
├── readme.md
├── sdy180
│   ├── protocols
│   │   └── Study Design.PTL5803.docx
│   ├── resultfiles
│   │   ├── gene_expression_result
│   │   │   └── Nanostring_norm_data_DS10_ESIDs_SDY180.587719.txt
│   │   ├── hai_result
│   │   │   └── Serology_by_Cohort.423738.xlsx
│   │   ├── mbaa_result.csv
│   │   ├── neutralizing_antibody_titer_result
│   │   │   └── Serology_by_Cohort.423737.xlsx
│   │   ├── sdy180-dr47_subject_2_gene_expression_result.txt
│   │   ├── sdy180-dr47_subject_2_hai_result.txt
│   │   ├── sdy180-dr47_subject_2_illumina_beadarray.txt
│   │   ├── sdy180-dr47_subject_2_mbaa_result.txt
│   │   ├── sdy180-dr47_subject_2_neutralizing_antibody_titer_result.txt
│   │   └── sdy180-dr47_subject_2_virus_neutralization_result.txt
│   ├── sdy180-dr47_tab
│   │   ├── arm_2_subject.csv
│   │   ├── arm_or_cohort.csv
│   │   ├── biosample.csv
│   │   ├── control_sample_2_file_info.csv
│   │   ├── control_sample.csv
│   │   ├── experiment_2_protocol.csv
│   │   ├── experiment.csv
│   │   ├── expsample_2_biosample.csv
│   │   ├── expsample_2_file_info.csv
│   │   ├── expsample_2_reagent.csv
│   │   ├── expsample_2_treatment.csv
│   │   ├── expsample.csv
│   │   ├── expsample_mbaa_detail.csv
│   │   ├── expsample_public_repository.csv
│   │   ├── hai_result.csv
│   │   ├── immune_exposure.csv
│   │   ├── inclusion_exclusion.csv
│   │   ├── lab_test.csv
│   │   ├── lab_test_panel_2_protocol.csv
│   │   ├── lab_test_panel.csv
│   │   ├── neut_ab_titer_result.csv
│   │   ├── planned_visit.csv
│   │   ├── protocol.csv
│   │   ├── reagent.csv
│   │   ├── reagent_set_2_reagent.csv
│   │   ├── standard_curve_2_file_info.csv
│   │   ├── standard_curve.csv
│   │   ├── study_2_condition_or_disease.csv
│   │   ├── study_2_panel.csv
│   │   ├── study_2_protocol.csv
│   │   ├── subject.csv
│   │   └── treatment.csv
│   └── studyfiles
│       ├── SDY180_BISC_Data_Curation_Notes_030614.doc
│       ├── SDY180_CBC_Results_and_Dictionary.xlsx
│       └── StudyFile.Obsolete.Period.StudyGlossary.ID4633-34.2799.SDY180.txt
├── sdy296
│   ├── protocols
│   │   └── Study_Design.PTL6569.pdf
│   ├── resultfiles (there are 37 patients, 13 with data in all 4 datasets)
│   │   ├── gene_expression_result (a 60 genes expression on n = 36 from nanostring, intersect of 13 with rna_seq)
│   │   │   └── Nanostring_norm_data_DS10_ESIDs_SDY296.587721.txt
│   │   ├── hai_result.csv (haemaglutin influenza antibody results n = 37)
│   │   ├── neut_ab_titer_result.csv (0, 28 days, n = 37)
│   │   ├── rna_sequencing_result (15573 genes, n = 42, intersect of 13 with nanostring)
│   │   │   └── SDY296_EXP13760_RNA_seq.703270.tsv
│   │   ├── sdy296-dr47_subject_2_gene_expression_result.txt (n = 37, refs files)
│   │   ├── sdy296-dr47_subject_2_illumina_beadarray.txt (n = 220, refs files)
│   │   └── sdy296-dr47_subject_2_rna_sequencing_result.txt (n = 43, refs files)
│   ├── sdy296-dr47_tab
│   │   ├── arm_2_subject.csv (age)
│   │   ├── arm_or_cohort.csv
│   │   ├── assessment_component.csv
│   │   ├── assessment_panel.csv
│   │   ├── biosample.csv
│   │   ├── experiment_2_protocol.csv
│   │   ├── experiment.csv
│   │   ├── expsample_2_biosample.csv
│   │   ├── expsample_2_file_info.csv
│   │   ├── expsample_2_reagent.csv
│   │   ├── expsample_2_treatment.csv
│   │   ├── expsample.csv
│   │   ├── expsample_public_repository.csv
│   │   ├── immune_exposure.csv
│   │   ├── inclusion_exclusion.csv
│   │   ├── lab_test.csv
│   │   ├── lab_test_panel_2_protocol.csv
│   │   ├── lab_test_panel.csv
│   │   ├── planned_visit.csv
│   │   ├── reagent.csv
│   │   ├── reagent_set_2_reagent.csv
│   │   ├── study_2_condition_or_disease.csv
│   │   ├── study_2_panel.csv
│   │   ├── subject.csv
│   │   └── treatment.csv
│   └── studyfiles
│       ├── SDY296_BISC_Data_Curation_Notes.docx
│       └── StudyFile.Obsolete.Period.StudyGlossary.ID4633-34.2799.SDY296.txt
└── sdy301
    ├── protocols
    │   └── Study_Design.PTL6595.pdf
    ├── resultfiles
    │   ├── gene_expression_result
    │   │   └── Nanostring_norm_data_DS10_ESIDs_SDY301.587720.txt
    │   ├── hai_result.csv
    │   ├── neut_ab_titer_result.csv
    │   ├── rna_sequencing_result
    │   │   └── SDY301_EXP13728_RNA_seq.703279.tsv
    │   ├── sdy301-dr47_subject_2_gene_expression_result.txt
    │   ├── sdy301-dr47_subject_2_illumina_beadarray.txt
    │   └── sdy301-dr47_subject_2_rna_sequencing_result.txt
    ├── sdy301-dr47_tab
    │   ├── arm_2_subject.csv
    │   ├── arm_or_cohort.csv
    │   ├── assessment_component.csv
    │   ├── assessment_panel.csv
    │   ├── biosample.csv
    │   ├── experiment_2_protocol.csv
    │   ├── experiment.csv
    │   ├── expsample_2_biosample.csv
    │   ├── expsample_2_file_info.csv
    │   ├── expsample_2_reagent.csv
    │   ├── expsample_2_treatment.csv
    │   ├── expsample.csv
    │   ├── expsample_public_repository.csv
    │   ├── immune_exposure.csv
    │   ├── inclusion_exclusion.csv
    │   ├── lab_test.csv
    │   ├── lab_test_panel_2_protocol.csv
    │   ├── lab_test_panel.csv
    │   ├── planned_visit.csv
    │   ├── reagent.csv
    │   ├── reagent_set_2_reagent.csv
    │   ├── study_2_condition_or_disease.csv
    │   ├── study_2_panel.csv
    │   ├── subject.csv
    │   └── treatment.csv
    └── studyfiles
        ├── SDY301_BISC_Data_Curation_Notes.docx
        └── StudyFile.Obsolete.Period.StudyGlossary.ID4633-34.2799.SDY301.txt

23 directories, 163 files

