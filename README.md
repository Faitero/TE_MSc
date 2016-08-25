# General structure

├── run.sh
└── scr
    ├── R
    │   ├── Alus_variations.R
    │   ├── Bootstraping_final.R
    │   ├── Data_Prepocesig.R
    │   ├── Extract_GRanges.R
    │   ├── Fig_4.R
    │   ├── Fig_5.R
    │   ├── Fig_6.R
    │   ├── Poly_U_by_distanceR.R
    │   └── bootstraping_alusexons_clasification.R
    ├── bash
    │   ├── Alu_motives.sh
    │   ├── Fig_7.motifs.sh
    │   ├── Fig_8.Xlinks.sh
    │   ├── Join_final_tables_toPublish.sh
    │   ├── get_tables_for CLIP.sh
    │   ├── lift_and_procces.sh
    │   └── xlinks_to_coverage.sh
    └── python
        ├── 3SS_distance_from_alu.py
        ├── BED2BED_no_counts_G&B.py
        ├── BEDgraph2BED.py
        ├── add_ID_to_bed.py
        ├── check_test_sequence.py
        ├── findLongestStrech.py
        ├── flankBEDpositionsStrandSpecific.py
        ├── get3SS_from_random_alu.py
        ├── get_aluexon_corrected-2nt_positive strand.py
        ├── get_aluexon_from_distance_from_alu2.py
        ├── get_best_3SS.py
        ├── get_fasta_species.py
        ├── get_start_position_from_bed_corrected.py
        ├── lift_over_specie.py
        ├── liftover_alus.py
        └── split_bed_record.py
