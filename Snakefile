# Process input data
rule all_by_all:
    input:
        "downloads/challenge_data/drug_synergy_data/ch1_train_combination_and_monotherapy.csv/ch1_train_combination_and_monoTherapy.csv",
        "downloads/challenge_data/drug_synergy_data/ch1_leaderboard_monotherapy.csv/ch1_leaderBoard_monoTherapy.csv",
        "downloads/challenge_data/drug_synergy_data/ch1_test_monotherapy.csv/ch1_test_monoTherapy.csv",
        "downloads/challenge_data/drug_synergy_data/ch2_leaderboard_monotherapy.csv/ch2_leaderBoard_monoTherapy.csv",
        "downloads/challenge_data/drug_synergy_data/ch2_test_monotherapy.csv/ch2_test_monoTherapy.csv"
    output:
        "notebooks/all_by_all/DRUG_PAIRS.tsv.bz2",
        "notebooks/all_by_all/DRUG_PAIRS_CL.tsv.bz2",
        "notebooks/all_by_all/ALL_TRAINING_DATA.tsv.bz2",
        "notebooks/all_by_all/ALL_TRAINING_DATA_GBD.tsv.bz2",
        "notebooks/all_by_all/ALL_TRAINING_DATA_GBDCL.tsv.bz2",
        "notebooks/all_by_all/ALL_TRAINING_DATA_PAIR.tsv.bz2",
        "notebooks/all_by_all/ALL_TRAINING_DATA_PAIR_GBDPP.tsv.bz2",
        "notebooks/all_by_all/ALL_TRAINING_DATA_PAIR_GBCLP.tsv.bz2",
        "notebooks/all_by_all/done"
    shell:
        "jupyter nbconvert -y --execute --to html --stdout './notebooks/all_by_all.ipynb' > logs/all_by_all.html"


# Map drugs to HGNC and Ensembl targets
rule drug_to_target:
    input:
        "downloads/challenge_data/drug_synergy_data/drug_info_release.csv/Drug_info_release.csv"
    output:
        "notebooks/drug_to_target/done"
    shell:
        "jupyter nbconvert -y --execute --to html --stdout './notebooks/drug_to_target.ipynb' > logs/drug_to_target.html"


# Map pairs of drugs to pairs of targets
rule drug_pair_to_target_pair:
    input:
        "notebooks/drug_to_target/done"
    output:
        "notebooks/drug_pair_to_target_pair/done"
    shell:
        "jupyter nbconvert -y --execute --to html --stdout './notebooks/drug_to_target.ipynb' > logs/drug_to_target.html"


# Map pairs of drugs to pairs of targets (STITCH)
rule drug_pair_to_target_pair_stitch:
    input:
        "notebooks/drug_to_target/done"
    output:
        "notebooks/drug_pair_to_target_pair/done_stitch"
    shell:
        "jupyter nbconvert -y --execute --to html --stdout './notebooks/drug_to_target.ipynb' > logs/drug_to_target.html"


# Clare's target level features
rule original_target_features_avg:
    input:
        "notebooks/drug_pair_to_target_pair/done"
    shell:
        "mysql -u strokach -h 192.168.6.19 az_dream_2015 < sql/original_target_features_avg.sql"

rule original_target_features_max:
    input:
        "notebooks/drug_pair_to_target_pair/done"
    shell:
        "mysql -u strokach -h 192.168.6.19 az_dream_2015 < sql/original_target_features_max.sql"

# Clare's target level features (STITCH)
rule original_target_features_avg:
    input:
        "notebooks/drug_pair_to_target_pair/done_stitch"
    shell:
        "mysql -u strokach -h 192.168.6.19 az_dream_2015 < sql/original_target_features_stitch_avg.sql"

rule original_target_features_max:
    input:
        "notebooks/drug_pair_to_target_pair/done_stitch"
    shell:
        "mysql -u strokach -h 192.168.6.19 az_dream_2015 < sql/original_target_features_stitch_max.sql"



### Sanger molecular data ###

# downloads/challenge_data/sanger_molecular_data/cell_info.csv

# Gene expression
rule cnv:
    input:
        "downloads/challenge_data/sanger_molecular_data/cnv/cnv_gene.csv/cnv_gene.csv",
        "downloads/challenge_data/sanger_molecular_data/cnv/cnv_segment.csv/cnv_segment.csv"

# Mutation
rule gex:
    input:
        "downloads/challenge_data/sanger_molecular_data/gex.csv/gex.csv"


# Methylation
rule methyl:
    input:
        "downloads/challenge_data/sanger_molecular_data/methyl/cpg_isle_level/methyl_ilse_beta.csv/methyl_ilse_beta.csv",
        "downloads/challenge_data/sanger_molecular_data/methyl/cpg_isle_level/methyl_ilse_m.csv/methyl_ilse_m.csv",
        "downloads/challenge_data/sanger_molecular_data/methyl/cpg_probe_level/methyl_probe_beta.csv/methyl_probe_beta.csv",
        "downloads/challenge_data/sanger_molecular_data/methyl/cpg_probe_level/methyl_probe_m.csv/methyl_probe_m.csv",
        "downloads/challenge_data/sanger_molecular_data/methyl/cpg_probe_level/probe_info.csv/probe_info.csv"

# Cosmic
rule mutations
    input:
        "downloads/challenge_data/sanger_molecular_data/mutations.csv/mutations.csv"



### External sources ###

# Harmonizome


# Lincs


####
# Machine learning
