import sys
import time
import pickle

import scipy as sp
import numpy as np
import pandas as pd
import sqlalchemy as sa

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import Imputer

import functions as fn


features_to_exclude = [
    'synergy_score', 'synergy_score_y', 'synergy_score_diff',
    'd_1', 'd_2', 'c', 'qa', 'source', 'unique_id',
    'index'
]


def get_feature_columns(df):
    feature_columns = [
        c for c in df.columns
        if (df[c].dtype in (int, float) and c not in features_to_exclude) or print(c)
    ]
    return feature_columns


def worker(clf, chunks, feature_columns=None):
    imputer = Imputer(copy=False)

    xval_results = []
    for chunk_i, (train, test_wsyn) in enumerate(chunks):
        print("\nWorking on chunk {}...".format(chunk_i))
        print("train: {}".format(train.shape))
        print("test: {}".format(test_wsyn.shape))
        sys.stdout.flush()

        print("train_sg")
        train_sg = fn.get_synergy_groups(train, pairwise=False)
        if feature_columns is None:
            feature_columns = get_feature_columns(train_sg)
        #X = train_sg[feature_columns].values
        #X = imputer.fit_transform(X)
        X = train_sg[feature_columns].fillna(0).values
        Y = train_sg['synergy_score'].values

        if hasattr(clf, 'weight_nonsyn'):
            print("Using weight {}".format(clf.weight_nonsyn))
            sample_weight = np.ones(Y.shape)
            sample_weight[Y < 20] = clf.weight_nonsyn
            clf.fit(X, Y, sample_weight=sample_weight)
        else:
            clf.fit(X, Y)

        # VALIDATION
        test = test_wsyn.drop('synergy_score', axis=1)

        # B
        test_sg = fn.get_synergy_groups(train, test, pairwise=False)
        #X_test = test_sg[feature_columns].values
        #X_test = imputer.transform(X_test)
        X_test = test_sg[feature_columns].fillna(0).values

        train_sg['synergy_score_pred'] = clf.predict(X)
        test_sg['synergy_score_pred'] = clf.predict(X_test)
        test_sg = (
            test_sg
            .drop('synergy_score', axis=1)
            .merge(test_wsyn[['d_1', 'd_2', 'c', 'synergy_score']], on=['d_1', 'd_2', 'c'])
        )
        xval_results.append(test_sg[['d_1', 'd_2', 'c', 'synergy_score', 'synergy_score_pred']])

        # Some remakrs
        print("train_sg: {}".format(train_sg.shape))
        print("test_sg: {}".format(test_sg.shape))
        print(set(feature_columns) - set(test_sg.columns))
        print(sp.stats.pearsonr(train_sg['synergy_score'], train_sg['synergy_score_pred']))
        print(sp.stats.pearsonr(test_sg['synergy_score'], test_sg['synergy_score_pred']))
        sys.stdout.flush()

    xval_results_final = pd.concat(xval_results, ignore_index=False)
    return xval_results_final


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--loss', type=str)
    parser.add_argument('--learning_rate', type=float)
    parser.add_argument('--n_estimators', type=int)
    parser.add_argument('--subsample', type=float)
    parser.add_argument('--max_features', type=float)
    parser.add_argument('--max_depth', type=int)
    parser.add_argument('--weight_nonsyn', type=float)
    parser.add_argument('--features', type=str)
    args = parser.parse_args()

    clf_options = {
        'loss': args.loss,
        'learning_rate': args.learning_rate,
        'n_estimators': args.n_estimators,
        'subsample': args.subsample,
        'max_features': args.max_features,
        'max_depth': args.max_depth,
        'random_state': np.random.RandomState(42)
    }
    clf = GradientBoostingRegressor(**clf_options)
    clf.weight_nonsyn = args.weight_nonsyn

    with open('machine_learning/chunks.pickle', 'rb') as ifh:
        chunks = pickle.load(ifh)

    features = args.features.split(',')
    xval_results_final = worker(clf, chunks, args.features.split(','))

    score1, score2 = fn.get_score_ch1(xval_results_final)
    print(score1, score2)

    df = pd.DataFrame(dict(
        **args.__dict__,
        **dict(score1=score1, score2=score2, nfeatures=len(features), last_feature=features[-1])
    ), index=[0])

    n_tries = 0
    while n_tries < 5:
        try:
            df.to_sql(
                'iterative_feature_addition_2',
                sa.create_engine('mysql://strokach:@192.168.6.19:3306/az_dream_2015_ml'),
                if_exists='append',
                index=False
            )
            break
        except Exception as e:
            print("Try number {}...".format(n_tries))
            print(e)
            time.sleep(1)
            n_tries += 1
    print("Done!")
