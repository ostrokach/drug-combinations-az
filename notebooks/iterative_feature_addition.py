# -*- coding: utf-8 -*-
import sys
import pickle
import subprocess
import time
import shlex

import numpy as np
import pandas as pd
import sqlalchemy as sa

from sklearn.ensemble import GradientBoostingRegressor

from common import submitjob

import tmp_xval


if __name__ == '__main__':
    current_features = [
        'synergy_score_gbdd_mean',
        'synergy_score_gbd_mean_mean',
        'synergy_score_gbdc_mean_mean',
        'synergy_score_gbc_mean',
        'h_mean',
        'gex_gbg_std_stitch_mean_diff',
        'bg_degree_max',
        'histology_subtype_mixed_adenosquamous_carcinoma',
        'f_mutation_maybe_bad_density_gbg_stitch_mean_diff',
        'gex_gbc_mean',
        'FingerprintMol_Kulczynski',
        'gene_essentiality_max',
        's_shortest_path_length_min',
        'methyl_probe_m_gbg_min_min_diff',
        'min_cn_gbc',
        'histology_subtype_NS',
        'feature_hydrophobe_count3_d_mean',
        'ic50_gbc_max',
        'gi_clustering_coef_max',
        'achilles_rnai_gs_mean_diff',
        'f_mutation_maybe_bad_density_gbgc_mean_diff',
        'methyl_probe_beta_drug_to_hgnc_target_min',
    ]

    clf_options = {
        'loss': 'huber',
        'learning_rate': 0.1,
        'n_estimators': 1000,
        'subsample': 0.3,
        'max_features': 1,
        'max_depth': 2,
        'random_state': np.random.RandomState(42)
    }
    clf = GradientBoostingRegressor(**clf_options)
    clf.weight_nonsyn = 0.3

    engine = sa.create_engine('mysql://strokach:@192.168.6.19:3306/az_dream_2015_ml')

    with open('machine_learning/chunks.pickle', 'rb') as ifh:
        chunks = pickle.load(ifh)

    all_features = tmp_xval.get_feature_columns(chunks[0][0])

    n_iterations = 0
    while True:
        n_iterations += 1
        remaining_features = list(set(all_features) - set(current_features))
        if not remaining_features:
            print("Done!")
            break
        for feature in remaining_features:
            system_command = """\
submitjob 1 -m 2 python tmp_xval.py --loss huber --learning_rate 0.1 --n_estimators 1000 \
--subsample 0.3 --max_features 1 --max_depth 2 --weight_nonsyn 0.3 --features '{}'
""".format(','.join(current_features + [feature]))
            print(system_command.strip())
            counter = 0
            while counter <= 5:
                try:
                    sp = subprocess.run(
                        shlex.split(system_command), stdout=sys.stdout, stderr=sys.stderr,
                        universal_newlines=True
                    )
                    sp.check_returncode()
                    break
                except Exception as e:
                    print(e)
                    time.sleep(60)
                    counter += 1
            time.sleep(0.5)
        running_jobs = submitjob.get_running_jobs()
        while len(running_jobs) > 3:
            print("running_jobs: {}".format(running_jobs))
            time.sleep(60)
            running_jobs = submitjob.get_running_jobs()
        current_features = pd.read_sql_query("""\
select features \
from iterative_feature_addition_2 \
order by nfeatures desc, score1 desc, score2 desc limit 1
""", engine)['features'].tolist()[0].split(',')
