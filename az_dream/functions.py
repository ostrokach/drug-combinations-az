import sys
import os
import os.path as op
import re
import warnings
from functools import partial

import numpy as np
import scipy as sp
import pandas as pd
import requests
from requests.utils import quote

from sklearn.ensemble import GradientBoostingRegressor
from IPython.display import display

import logging

logger = logging.getLogger()


def interquantile_range(x):
    return np.subtract(*np.percentile(x, [75, 25]))

STATS = ['max', 'min', 'mean', 'std']

# %%
metrics = [
    'euclidean',
    'minkowski',
    'cityblock',
    'seuclidean',
    'sqeuclidean',
    'cosine',
    'correlation',
    'hamming',
    'jaccard',
    'chebyshev',
    'canberra',
    'braycurtis',
    # 'mahalanobis',
    'yule',
    'matching',
    'dice',
    'kulsinski',
    'rogerstanimoto',
    'russellrao',
    'sokalmichener',
    'sokalsneath',
    # 'wminkowski',
]


def get_best_metric(data, columns, scores, feature_subset_name, row_subset_name=None):
    """
    Parameters
    ----------
    data : Array
        m: labels
        n: features
    columns : list
        Name of each row in `data`.
    scores : Series
        index: labels
        values: features
    feature_subset_name : str
        Name of the feature that we are focusing on.
        (Useful if you have multiple DataFrames).
    row_subset_name : str
        Name of the rows that we are focusing on.
    """
    index = pd.Index(columns, name='x')
    columns = pd.Index(columns, name='y')

    # Get a correlation score using all available metrics
    results = []
    for metric in metrics:
        feature_column_name = (
            feature_subset_name +
            (('_' + 'ss' + row_subset_name) if row_subset_name else '') +
            '_' + metric
        )
        try:
            PW = sp.spatial.distance.squareform(
                sp.spatial.distance.pdist(data, metric=metric)
            )
        except ValueError:
            print("Skipping metric {}".format(metric))
        PW = (
            pd.DataFrame(PW, index=index, columns=columns)
            .unstack()
            .reset_index()
            .rename(columns={0: feature_column_name})
        )
        PW = PW[PW['x'] != PW['y']]
        if '|' in scores.index[0]:
            PW['xy'] = PW['x'] + '|' + PW['y']
            PW['score_diff'] = PW['xy'].map(scores)
        else:
            PW['score_x'] = PW['x'].map(scores)
            PW['score_y'] = PW['y'].map(scores)
            PW['score_diff'] = PW['score_y'] - PW['score_x']
        PW = PW.reindex_axis(
            ['x', 'y', feature_column_name, 'score_diff'],
            axis='columns')
        # Pearson
        tmp = PW[[feature_column_name, 'score_diff']].dropna()
        corr, pvalue = sp.stats.pearsonr(
            tmp[feature_column_name],
            tmp['score_diff'].abs()
        )
        count = len(tmp)
        results.append((metric, corr, pvalue, count, PW, ))
    results_df = pd.DataFrame(results, columns=['metric', 'corr', 'pvalue', 'count', 'PW'])
    results_df['feature_subset_name'] = feature_subset_name
    results_df['row_subset_name'] = row_subset_name
    results_df['corr_abs'] = results_df['corr'].abs()
    results_df.sort_values('corr_abs', ascending=False, inplace=True, na_position='last')
    return results_df


def get_dd_distance(df, d2g, dss, name):
    """
    Parameters
    ----------
    df : DataFrame
        m: genes
        n: features
    d2g : dictionary
        keys: drugs
        values: genes that they target
    dss : Series
        index: drugs
        values: average synergy score

    Example
    -------
    fn.get_dd_distance(
        df[columns].set_index('g_1', 'g_2'),
        d2g,
        dss.set_index('xy')['score'],
        'mydf')

    """
    drug_data = dict()
    for drug, targets in d2g.items():
        df_drug = df.loc[df.index.isin(targets), :]
        if not df_drug.empty:
            drug_data[drug] = df_drug.std(axis=0)
    df_gbd = pd.DataFrame(drug_data).T
    #
    results = get_best_metric(
        data=df_gbd.values,
        columns=df_gbd.index.tolist(),
        scores=dss,
        feature_subset_name=name
    )
    return results, df_gbd


def get_cc_distance(df, css, row_subsets=None, column_subsets=None):
    """
    pa : Panel
        Typically ``Genes x Cell Lines x Features``.
    css : Series
        Column synergy scores.
        Typically ``Cell Line`` synergy scores
    """
    # Rows
    if row_subsets is None:
        row_subsets = {'': df.index}
    else:
        for key in row_subsets:
            if row_subsets[key] is None:
                row_subsets[key] = df.index
    # Columns
    if column_subsets is None:
        if isinstance(df.columns, pd.MultiIndex):
            column_subsets = {c: (c, ) for c in df.columns.levels[0]}
            columns = list(df.columns.levels[1])
        else:
            column_subsets = {'': None}
            columns = list(df.columns)
    print("column_subsets: {}".format(column_subsets))
    # Run!
    results_all = []
    for row_subset_name, row_subset in row_subsets.items():
        for feature_subset_name, column_subset in column_subsets.items():
            # Create data array
            if column_subset is None:
                data = df.loc[df.index.isin(row_subset), :].values.T
            else:
                data = np.vstack([
                    df.loc[df.index.isin(row_subset), c].values
                    for c in column_subset
                ]).T
            # Get a correlation score using all available metrics
            results_df = get_best_metric(data, columns, css, feature_subset_name, row_subset_name)
            # Save best result
            results_all.append(results_df[:1])
    # Final
    results_all_df = pd.concat(results_all, ignore_index=False)
    return results_all_df


# %%
def get_correlation(df, one_column, all_columns=None):
    """Calculate an all-by-one correlation matrix.
    """
    if all_columns is None:
        all_columns = [c for c in df.columns if df[c].dtype in (int, float) and c != one_column]
    results = []
    for column in all_columns:
        df_tmp = df[[column, one_column]].dropna(how='any')
        count = df_tmp.shape[0]
        corr, pvalue = sp.stats.pearsonr(df_tmp[column], df_tmp[one_column])
        results.append([column, corr, pvalue, count])
    # Don't sort the list because NaNs are not handled properly
    results_df = pd.DataFrame(results, columns=['feature', 'correlation', 'pvalue', 'count'])
    results_df['correlation_abs'] = results_df['correlation'].abs()
    results_df.sort_values(by='correlation_abs', ascending=False, inplace=True, na_position='last')
    return results_df


# %%
def my_groupby(df, by, features=None, stats=STATS, name=''):
    """
    all_by_all
    """
    grouped = df.groupby(by)
    if features is None:
        df = grouped.agg(stats)
    else:
        df = grouped[features].agg(stats)
    df[('count', '')] = grouped[df.columns[0][0]].count().values
    df.columns = ['_'.join([c[0], name, c[1]]).strip('_').lower() for c in df.columns]
    df.reset_index(inplace=True)
    return df


def loo_groupby(df, by, features, stats, name=''):
    """
    Perfor the equivalent of `df.groupby(by)[features].agg(stats)`,
    with with the `.agg` not including the row in question.

    Parameters
    ----------
    df : DataFrame
        Data to analyse.
    on : list
        Columns to group on.
    features : list
        Columns to analyse.
    Stats : list
        Statistics to calculate.
    """
    def get_fn(name):
        """Returns numpy function with the name `name`.
        """
        return getattr(np, name)

    results = []
    if 'index' in df:
        ic = ['index']
    elif 'index_1' in df and 'index_2' in df:
        ic = ['index_1', 'index_2']
    elif 'index_x' in df and 'index_y' in df:
        ic = ['index_x', 'index_y']
    else:
        df = df.reset_index()
        assert 'index' in df
        ic = ['index']

    df = df.sort_values(ic + by)

    grouped = df.groupby(by)
    for key, group in grouped:
        for idx, row in group.iterrows():
            # Index
            if len(ic) == 1:
                group_other = group[group[ic[0]] != row[ic[0]]]
                idx = [row[ic[0]]]
            else:
                group_other = group[
                    (group[ic[0]] != row[ic[0]]) &
                    (group[ic[0]] != row[ic[1]]) &
                    (group[ic[1]] != row[ic[0]]) &
                    (group[ic[1]] != row[ic[1]])
                ]
                idx = [row[ic[0]], row[ic[1]]]
            # Statistics
            if group_other.empty:
                feature_stats = [np.nan for i in range(len(features) * len(stats))]
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    feature_stats = [
                        get_fn(stat)(group_other[feature].values)
                        for feature in features
                        for stat in stats
                    ]
            count = group_other[features[0]].count()
            # Combine
            results.append(
                idx +
                (list(key) if isinstance(key, (list, tuple)) else [key]) +
                [count] +
                feature_stats
            )

    # All together
    columns = (
        ic +
        by +
        ['count_' + name] +
        ['_'.join([feature, name, stat.replace('nan', '')])
         for feature in features for stat in stats]
    )
    results_df = pd.DataFrame(results, columns=columns)
    results_df = results_df.sort_values(ic + by)
    results_df.index = df.index
    assert (df[ic + by] == results_df[ic + by]).all().all()
    return results_df


def get_differences(df, keep=None):
    """
    Calculate differences for a pair of drugs (ends with _1 and _2)
    or for a pair of cell lines (ends with _x and _y).
    """
    df = df.copy()
    columns = df.columns.copy()
    columns_to_drop = []
    for column in columns:
        if not isinstance(column, str):
            continue
        elif column.endswith('_1'):
            suffixes = ('_1', '_2')
        elif column.endswith('_x'):
            suffixes = ('_x', '_y')
        elif column.endswith('_a'):
            suffixes = ('_a', '_b')
        else:
            continue
        if df[column].dtype not in [int, float]:
            # print("Skipping column '{}' because it appears to be a string...".format(column))
            continue
        column_2 = column.replace(suffixes[0], suffixes[1])
        if (column_2 == column) or (column_2 not in df.columns):
            continue
        if column.startswith('synergy_score'):
            # continue
            if column == 'synergy_score_x':
                df[column.replace(suffixes[0], '_diff')] = df[column_2] - df[column]
                # keep synergy_score_x
            elif column.endswith('_1'):
                df[column.replace(suffixes[0], '_mean')] = (df[column_2] + df[column]) / 2
                df[column.replace(suffixes[0], '_diff')] = (df[column_2] - df[column]).abs()
                columns_to_drop.extend([column, column_2])
            elif column.endswith('_x'):
                df[column.replace(suffixes[0], '_diff')] = df[column_2] - df[column]
                columns_to_drop.extend([column])
            else:
                raise Exception
        else:
            if column.endswith('_1'):
                df[column.replace(suffixes[0], '_mean')] = (df[column_2] + df[column]) / 2
                df[column.replace(suffixes[0], '_diff')] = (df[column_2] - df[column]).abs()
                if keep is None:
                    # Default
                    columns_to_drop.extend([column, column_2])
                elif keep == -1:
                    pass
                elif keep == 0:
                    columns_to_drop.extend([column, column_2])
                elif keep == 1:
                    columns_to_drop.extend([column_2])
                elif keep == 2:
                    columns_to_drop.extend([column])
            elif column.endswith('_x'):
                df[column.replace(suffixes[0], '_diff')] = (df[column_2] - df[column])
                if keep is None:
                    # Default
                    columns_to_drop.extend([column])
                elif keep == -1:
                    pass
                elif keep == 0:
                    columns_to_drop.extend([column, column_2])
                elif keep == 1:
                    columns_to_drop.extend([column_2])
                elif keep == 2:
                    columns_to_drop.extend([column])
            else:
                raise Exception
    df.drop(pd.Index(columns_to_drop), axis=1, inplace=True)
    return df


def get_synergy_groups(train_df, test_df, pairwise_data=None):
    """
    """
    STATS = ['median']
    this_groupby = partial(my_groupby, features=['synergy_score'], stats=STATS)

    if 'synergy_score' in test_df:
        print("`df_test` contains a 'synergy_score' column. "
              "This should not be done in production!")
    else:
        test_df['synergy_score'] = np.nan

    # Group by drug pair (only for ch1)
    df_dd = this_groupby(train_df, by=['d_1', 'd_2'], name='gbdd')
    # Group by cell line
    df_c = this_groupby(train_df, by=['c'], name='gbc')
    # Group by drug / cell line
    df_stack = pd.concat([
        train_df[['d_1', 'c', 'synergy_score']].rename(columns={'d_1': 'd'}),
        train_df[['d_2', 'c', 'synergy_score']].rename(columns={'d_2': 'd'}),
    ])
    df_dc = this_groupby(df_stack, by=['d', 'c'], name='gbdc')
    df_d = this_groupby(df_stack, by=['d'], name='gbd')

    # Combine
    test_df_wsyn = (
        test_df
        .merge(
            df_dd,
            on=['d_1', 'd_2'], suffixes=('', '_dup'), how='left')
        .merge(
            df_c,
            on=['c'], suffixes=('', '_dup'), how='left')
        .merge(
            df_dc.rename(columns={'d': 'd_1'}),
            on=['d_1', 'c'], suffixes=('', '_dup'), how='left')
        .merge(
            df_dc.rename(columns={'d': 'd_2'}),
            on=['d_2', 'c'], suffixes=('_1', '_2'), how='left')
        .merge(
            df_d.rename(columns={'d': 'd_1'}),
            on=['d_1'], suffixes=('', '_dup'), how='left')
        .merge(
            df_d.rename(columns={'d': 'd_2'}),
            on=['d_2'], suffixes=('_1', '_2'), how='left')
    )
    test_df_wsyn = get_differences(test_df_wsyn)

    # No test_df assertions
    assert (test_df[['d_1', 'd_2', 'c']] == test_df_wsyn[['d_1', 'd_2', 'c']]).all().all()
    assert test_df.shape[0] == test_df_wsyn.shape[0]

    if pairwise_data is None:
        return test_df_wsyn

    train_df['unique_id'] = train_df['d_1'] + '.' + train_df['d_2'] + '.' + train_df['c']
    test_df_wsyn['unique_id'] = test_df['d_1'] + '.' + test_df['d_2'] + '.' + test_df['c']

    tmp = train_df[['d_1', 'd_2', 'c', 'unique_id', 'synergy_score']].copy()
    df_pair = tmp.merge(tmp, on=['d_1', 'd_2'], suffixes=('_x', '_y'))
    df_pair['synergy_score_diff'] = (df_pair['synergy_score_y'] - df_pair['synergy_score_x'])
    df_pair = df_pair[df_pair['unique_id_x'] != df_pair['unique_id_y']]
    df_pair.drop(
        pd.Index(['synergy_score_y', 'synergy_score_x', 'unique_id_x', 'unique_id_y']),
        axis=1, inplace=True)

    df_pair_gbcc = my_groupby(
        df_pair, by=['c_x', 'c_y'], name='gbcc', features=['synergy_score_diff'], stats=STATS)

    test_df_pair = (
        train_df
        .merge(test_df_wsyn, on=['d_1', 'd_2'], suffixes=('_x', '_y'), how='right')
        .merge(df_pair_gbcc, on=['c_x', 'c_y'], how='left')
        .merge(pairwise_data, on=['c_x', 'c_y'], how='left')
        .rename(columns={'c_y': 'c', 'synergy_score_y': 'synergy_score'})
    )
    test_df_pair = test_df_pair[test_df_pair['unique_id_x'] != test_df_pair['unique_id_y']]
    test_df_pair['synergy_score_gbcc_median'] = (
        test_df_pair['synergy_score_x'] + df_pair_gbcc['synergy_score_diff_gbcc_median'])
    test_df_pair.drop(pd.Index(['unique_id_x', 'unique_id_y']), axis=1, inplace=True)

    test_df_pair_diff = get_differences(test_df_pair)
    return test_df_pair_diff


def get_synergy_groups_obsolete(df, df_test=None, pairwise=False):
    """
    .. note::

        Obsolete. Do not use.
    """
    iqr = lambda x: np.subtract(*np.percentile(x, [75, 25]))
    STATS = ['median']
    if df_test is None:
        # stats = [('nan' + s if isinstance(s, str) else s) for s in STATS]
        stats = ['median']
        # Returns DF of the same shape as the original
        this_groupby = partial(loo_groupby, features=['synergy_score'], stats=stats)
        ic = ['index']
        df = df.reset_index()
    else:
        this_groupby = partial(my_groupby, features=['synergy_score'], stats=STATS)
        ic = []
        if 'synergy_score' in df_test:
            print("`df_test` contains a 'synergy_score' column. "
                  "This should not be done in production!")
        else:
            df_test['synergy_score'] = np.nan

    # Group by drug pair (only for ch1)
    df_dd = this_groupby(df, by=['d_1', 'd_2'], name='gbdd')
    # Group by cell line
    df_c = this_groupby(df, by=['c'], name='gbc')
    # Group by drug / cell line
    df_stack = pd.concat([
        df[ic + ['d_1', 'c', 'synergy_score']].rename(columns={'d_1': 'd'}),
        df[ic + ['d_2', 'c', 'synergy_score']].rename(columns={'d_2': 'd'}),
    ])
    df_dc = this_groupby(df_stack, by=['d', 'c'], name='gbdc')
    df_d = this_groupby(df_stack, by=['d'], name='gbd')

    # Combine
    def _get_df_wsyn(df):
        df_wsyn = (
            df
            .merge(
                df_dd,
                on=ic + ['d_1', 'd_2'], suffixes=('', '_dup'), how='left')
            .merge(
                df_c,
                on=ic + ['c'], suffixes=('', '_dup'), how='left')
            .merge(
                df_dc.rename(columns={'d': 'd_1'}),
                on=ic + ['d_1', 'c'], suffixes=('', '_dup'), how='left')
            .merge(
                df_dc.rename(columns={'d': 'd_2'}),
                on=ic + ['d_2', 'c'], suffixes=('_1', '_2'), how='left')
            .merge(
                df_d.rename(columns={'d': 'd_1'}),
                on=ic + ['d_1'], suffixes=('', '_dup'), how='left')
            .merge(
                df_d.rename(columns={'d': 'd_2'}),
                on=ic + ['d_2'], suffixes=('_1', '_2'), how='left')
        )
        return df_wsyn

    df_wsyn = _get_df_wsyn(df)
    df_wsyn = get_differences(df_wsyn)
    # No test_df assertions
    assert (
        df[ic + ['d_1', 'd_2', 'c']] ==
        df_wsyn[ic + ['d_1', 'd_2', 'c']]
    ).all().all()
    assert df.shape[0] == df_wsyn.shape[0]

    if df_test is None:
        df_test_wsyn = df_wsyn.copy()
    else:
        df_test_wsyn = _get_df_wsyn(df_test)
        df_test_wsyn = get_differences(df_test_wsyn)
        # test_df assertions
        assert (
            df_test[ic + ['d_1', 'd_2', 'c']] ==
            df_test_wsyn[ic + ['d_1', 'd_2', 'c']]
        ).all().all()
        assert df_test.shape[0] == df_test_wsyn.shape[0]

    if not pairwise:
        return df_test_wsyn

    # Pair ===================================================================
    STATS = ['median']
    if df_test is None:
        stats = [('nan' + s) for s in STATS]
        this_groupby = partial(loo_groupby, features=['synergy_score_diff'], stats=stats)
        ic = ['index_x', 'index_y']
    else:
        this_groupby = partial(my_groupby, features=['synergy_score_diff'], stats=STATS)

    print('1...')
    sys.stdout.flush()
    df_tmp = df.merge(df, on=['d_1', 'd_2'], suffixes=('_x', '_y'))
    df_tmp['synergy_score_diff'] = (df_tmp['synergy_score_y'] - df_tmp['synergy_score_x'])
    df_gbcc = this_groupby(df_tmp, by=['c_x', 'c_y'], name='gbcc')
    df_gbcc.drop('count_gbcc', axis=1)

    print('2...')
    sys.stdout.flush()
    df_tmp = df.merge(df, on=['c'], suffixes=('_x', '_y'))
    df_tmp['synergy_score_diff'] = (df_tmp['synergy_score_y'] - df_tmp['synergy_score_x'])
    df_gbdddd = this_groupby(df_tmp, by=['d_1_x', 'd_1_y', 'd_2_x', 'd_2_y'], name='gbdddd')
    df_gbdddd.drop('count_gbdddd', axis=1)

    print('3...')
    sys.stdout.flush()
    df_tmp = df_stack.merge(df_stack, on=['c'], suffixes=('_x', '_y'))
    df_tmp['synergy_score_diff'] = (df_tmp['synergy_score_y'] - df_tmp['synergy_score_x'])
    df_gbdcdc = this_groupby(df_tmp, by=['d_x', 'd_y'], name='gbdcdc')
    df_gbdcdc.drop('count_gbdcdc', axis=1)

    print('4...')
    sys.stdout.flush()
    df_wsyn['tmp'] = 1
    df_test_wsyn['tmp'] = 1
    df_wsyn_pair = df_wsyn.merge(df_test_wsyn, on=['tmp'], suffixes=('_x', '_y'))

    print('5...')
    sys.stdout.flush()
    df_wsyn_pair = (
        df_wsyn_pair
        .merge(
            df_gbcc, on=ic + ['c_x', 'c_y'], how='left')
        .merge(
            df_gbdddd, on=ic + ['d_1_x', 'd_1_y', 'd_2_x', 'd_2_y'], how='left')
        .merge(
            df_gbdcdc.rename(columns={'d_x': 'd_1_x', 'd_y': 'd_1_y'}),
            on=ic + ['d_1_x', 'd_1_y'], how='left')
        .merge(
            df_gbdcdc.rename(columns={'d_x': 'd_1_x', 'd_y': 'd_2_y'}),
            on=ic + ['d_1_x', 'd_2_y'], suffixes=('_11', '_12'), how='left')
        .merge(
            df_gbdcdc.rename(columns={'d_x': 'd_2_x', 'd_y': 'd_1_y'}),
            on=ic + ['d_2_x', 'd_1_y'], how='left')
        .merge(
            df_gbdcdc.rename(columns={'d_x': 'd_2_x', 'd_y': 'd_2_y'}),
            on=ic + ['d_2_x', 'd_2_y'], suffixes=('_21', '_22'), how='left')
    )

    print('6...')
    sys.stdout.flush()
    df_wsyn_pair['synergy_score_gbcc_mean'] = (
        df_wsyn_pair['synergy_score_x'] + df_wsyn_pair['synergy_score_diff_gbcc_mean']
    )
    # df_wsyn_pair.drop('synergy_score_diff_gbcc_mean', axis=1)

    df_wsyn_pair['synergy_score_gbdddd_mean'] = (
        df_wsyn_pair['synergy_score_x'] + df_wsyn_pair['synergy_score_diff_gbdddd_mean']
    )
    # df_wsyn_pair.drop('synergy_score_diff_gbdddd_mean', axis=1)

    df_wsyn_pair['synergy_score_diff_gbdcdc_mean'] = (
        df_wsyn_pair['synergy_score_diff_gbdcdc_mean_11'] +
        df_wsyn_pair['synergy_score_diff_gbdcdc_mean_12'] +
        df_wsyn_pair['synergy_score_diff_gbdcdc_mean_21'] +
        df_wsyn_pair['synergy_score_diff_gbdcdc_mean_22']
    ) / 4
    # TODO: do something about '_max', '_min', '_mean', '_std' as well
    df_wsyn_pair.drop(
        pd.Index([c for c in df_wsyn_pair.columns if c[-3:] in ['_11', '_12', '_21', '_22']]),
        axis=1, inplace=True)
    df_wsyn_pair['synergy_score_gbdcdc_mean'] = (
        df_wsyn_pair['synergy_score_x'] + df_wsyn_pair['synergy_score_diff_gbdcdc_mean']
    )
    # df_wsyn_pair.drop('synergy_score_diff_gbdcdc_mean', axis=1)

    print('7...')
    sys.stdout.flush()
    df_wsyn_pair_diff = get_differences(df_wsyn_pair)
    return df_wsyn_pair_diff


# %% OLD STUFF
def get_weighted_stats(
        df, group_by_column, feature_column, weight_column,
        mean_range=(-100, 100), std_range=(0, 100)):
    """
    """
    def get_average(x):
        try:
            return np.average(x[feature_column].values, weights=x[weight_column])
        except ZeroDivisionError as e:
            print('The following error occured:\n{}'.format(e))
            return np.average(x[x[weight_column] > 0][feature_column].values)

    def get_std(x):
        try:
            return np.sqrt(np.average(
                (x[feature_column].values - get_average(x))**2, weights=x[weight_column]
            ))
        except ZeroDivisionError as e:
            print('The following error occured:\n{}'.format(e))
            return np.sqrt(np.average(
                (x[x[weight_column] > 0][feature_column].values - get_average(x))**2
            ))

    df_gp = df.groupby(group_by_column)

    df_2 = pd.DataFrame(
        data={
            'MEAN': df_gp.apply(get_average),
            'STD': df_gp.apply(get_std),
            'LEN': df_gp[weight_column].apply(np.sum),
            'MEAN_UNWEIGHTED': df_gp[feature_column].apply(np.mean),
            'STD_UNWEIGHTED': df_gp[feature_column].apply(np.std),
            'LEN_UNWEIGHTED': df_gp[feature_column].apply(len),
        },
    )
    df_2 = df_2.reset_index()

    # Normalize the mean
    df_2.loc[(df_2['MEAN'] <= mean_range[0]), 'MEAN'] = mean_range[0]
    df_2.loc[(df_2['MEAN'] >= mean_range[1]), 'MEAN'] = mean_range[1]

    # Normalize the standard deviation
    df_2['STE'] = df_2['STD']
    df_2.loc[(df_2['STD'] <= std_range[0]), 'STD'] = std_range[0]
    df_2.loc[(df_2['STD'] >= std_range[1]), 'STD'] = std_range[1]
    # Standard error
    df_2.loc[(df_2['STE'] <= std_range[0]), 'STE'] = np.nan
    df_2.loc[(df_2['STE'] >= std_range[1]), 'STE'] = np.nan
    temp_len = df_2['LEN'].apply(lambda x: max(x, 1))
    df_2['STE'] = df_2['STE'] / np.sqrt(temp_len)
    df_2.loc[pd.isnull(df_2['STE']), 'STE'] = np.nanmax(df_2['STE'])
    df_2['STE_UNWEIGHTED'] = df_2['STD_UNWEIGHTED'] / np.sqrt(df_2['LEN_UNWEIGHTED'])
    # Confidence score
    df_2['log_ste'] = np.log(df_2['STE'].values)
    df_2['norm_log_ste'] = (
        (df_2['log_ste'] - df_2['log_ste'].min()) /
        (df_2['log_ste'].max() - df_2['log_ste'].min())
    )
    # The higher the normalized standard error, the lower the confidence
    df_2['CONFIDENCE'] = 1 - df_2['norm_log_ste']

    df_2.drop(['log_ste', 'norm_log_ste'], axis=1, inplace=True)

    return df_2


# %%
def cross_validate(clf, df, x_columns, y_column):
    combination_ids = set(df['COMBINATION_ID'])
    for combination_id in combination_ids:
        print('Cross validating {}...'.format(combination_id))
        x_train = df[df['COMBINATION_ID'] != combination_id][x_columns].values
        y_train = df[df['COMBINATION_ID'] != combination_id][y_column].values
        clf.fit(x_train, y_train)
        x_test = df[df['COMBINATION_ID'] == combination_id][x_columns].values
        df.loc[df['COMBINATION_ID'] == combination_id, 'prediction'] = clf.predict(x_test)


# %%
class CrossValidate:

    def __init__(self, clf_options, df, x_columns, y_column):
        self.clf_options = clf_options
        self.df = df
        self.x_columns = x_columns
        self.y_column = y_column

    def run(self):
        self.clf_options['random_state'] = np.random.RandomState(42)
        clf = GradientBoostingRegressor(**self.clf_options)
        cross_validate(clf, self.df, self.x_columns, self.y_column)

    def save(self, filename):
        self.df['prediction'].to_pickle(filename)


# %%
pug = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/'
properties = [
    'MolecularFormula',
    'MolecularWeight',
    'CanonicalSMILES',
    'IsomericSMILES',
    'InChI',
    'InChIKey',
    'IUPACName',
    'XLogP',
    'ExactMass',
    'MonoisotopicMass',
    'TPSA',
    'Complexity',
    'Charge',
    'HBondDonorCount',
    'HBondAcceptorCount',
    'RotatableBondCount',
    'HeavyAtomCount',
    'IsotopeAtomCount',
    'AtomStereoCount',
    'DefinedAtomStereoCount',
    'UndefinedAtomStereoCount',
    'BondStereoCount',
    'DefinedBondStereoCount',
    'UndefinedBondStereoCount',
    'CovalentUnitCount',
    'Volume3D',
    'XStericQuadrupole3D',
    'YStericQuadrupole3D',
    'ZStericQuadrupole3D',
    'FeatureCount3D',
    'FeatureAcceptorCount3D',
    'FeatureDonorCount3D',
    'FeatureAnionCount3D',
    'FeatureCationCount3D',
    'FeatureRingCount3D',
    'FeatureHydrophobeCount3D',
    'ConformerModelRMSD3D',
    'EffectiveRotorCount3D',
    'ConformerCount3D',
    'Fingerprint2D'
]
properties_string = ','.join(properties)


class NotFoundError(Exception):
    pass


def _get_properties(qtype, key):
    """
    Parameters
    ----------
    qtype : str
        {'cid', 'smiles'}
    """
    print(key)
    url = (
        pug + 'compound/{}/{}/property/{}/JSON'
        .format(qtype, quote(key, safe=''), properties_string)
    )
    r = requests.get(url)
    # print(url)
    try:
        properties = r.json()
    except ValueError as e:
        print(e)
        print(r.reason)
        if r.reason not in ['Not Found']:
            print(url)
        raise NotFoundError
    try:
        properties = properties['PropertyTable']['Properties']
    except KeyError as e:
        print(e)
        print(properties)
        print(url)
        raise NotFoundError
    assert len(properties) == 1
    return properties[0]


def get_drug_features(df):
    """Query the PubChem REST API for drug properties.

    https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html
    """
    result_list = []
    for idx, row in df.iterrows():
        key = row['SMILES or PubChem ID']
        print()
        print(key)
        print('-' * 80)
        if pd.isnull(key):
            result_list.append(row)
            continue
        elif key.isnumeric():
            keys = [key]
            qtype = 'cid'
        else:
            qtype = 'smiles'
            keys = key.split(';')

        for key in keys:
            try:
                properties = _get_properties(qtype, key)
            except NotFoundError:
                continue
            row = row.append(pd.Series(properties))
            break
        result_list.append(row)
    return pd.concat(result_list, axis=1).transpose()


# %%
def get_drug_features_2(cids):
    cids = {int(c) for c in cids if pd.notnull(c)}
    url_template = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON'
    gene_names = {cid: [] for cid in cids}
    atc_codes = {cid: [] for cid in cids}
    for cid in cids:
        url = url_template.format(cid)
        r = requests.get(url)
        result = r.json()
        try:
            records = result['Record']['Section']
        except KeyError:
            print(result.keys())
            print(r.reason)
            print(url)
            raise
        return records
        for record in records:
            print(record)
            if 'TOCHeading' not in record:
                print(record)
                break

            if record['TOCHeading'] == 'Gene Name':
                print(record)
                for info in record['Information']:
                    if info['Name'] == 'Gene Name':
                        gene_names[cid].append(info['StringValue'])

            if record['TOCHeading'] == 'ATC Code':
                print(record)
                for info in record['Information']:
                    if info['Name'] == 'ATC Code':
                        atc_code = info['StringValue'].split()[0]
                        atc_code_2 = (
                            info['URL'][
                                info['URL'].index('?code=') + 6:
                                info['URL'].index('&')
                            ]
                        )
                        assert atc_code == atc_code_2
                        atc_codes[cid].append(atc_code)
    return gene_names, atc_codes


# %%
def p3print(recs, n_spaces):
    for rec in recs:
        if type(rec) == dict:
            if 'TOCHeading' in rec:
                print(' ' * n_spaces + 'TOCHeading: {}'.format(rec['TOCHeading']))
            if 'Description' in rec:
                print(' ' * n_spaces + 'Description: {}'.format(rec['Description']))
            if 'Information' in rec:
                print(' ' * (n_spaces + 4) + 'Information')
                for info in rec['Information']:
                    p3print(info, n_spaces+4)
            if 'Section' in rec:
                print(' ' * (n_spaces + 4) + 'Section')
                for section in rec['Section']:
                    p3print(section, n_spaces+4)
        else:
            print(rec)


# %%
def get_cids(cid):
    """
    https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html
    """
    similar_cids_url = (
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/cids/JSON'
        '?cids_type=same_connectivity'
    )
    url = similar_cids_url.format(cid)
    r = requests.get(url)
    all_cids = r.json()['IdentifierList']['CID']
    return all_cids


def get_cid_website(cid):
    """
    https://pubchem.ncbi.nlm.nih.gov/
    """
    cid_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON/'
    url = cid_url.format(cid)
    r = requests.get(url)
    return r.json()['Record']['Section']


def get_sids(cids):
    similar_sids_url = (
        'http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/sids/JSON'
    )
    chunk_size = 20
    sids = []
    for i in list(range(0, len(cids), chunk_size))[:-1]:
        cids_chunk = cids[i:i + chunk_size]
        url = similar_sids_url.format(','.join(str(cid) for cid in cids_chunk))
        r = requests.get(url)
        for info in r.json()['InformationList']['Information']:
            try:
                cid = info['CID']
                for sid in info['SID']:
                    sids.append((cid, sid))
            except KeyError:
                if 'CID' not in info or len(info) > 1:
                    raise
    return sids


def get_sid_data(sids):
    dbs_to_keep = [
        'ChEMBL', 'DrugBank', 'BindingDB', 'Broad Institute', 'KEGG', 'IUPHAR-DB',
    ]
    chunk_size = 20
    result = []
    for i in list(range(0, len(sids), chunk_size))[:-1]:
        sids_chunk = sids[i:i + chunk_size]
        sid_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{}/JSON'
        url = sid_url.format(','.join(str(sid) for sid in sids_chunk))
        r = requests.get(url)
        for info in r.json()['PC_Substances']:
            sid = info['sid']['id']
            db_name = info['source']['db']['name']
            db_id = info['source']['db']['source_id']['str']
            if db_name in dbs_to_keep:
                result.append((sid, db_name, db_id))
    return result


def get_atcs(cids):
    names = []
    mesh_synonyms = []
    depositor_supplied_synonyms = []
    atcs = []
    targets = []
    general_functions = []
    #
    remove_links = re.compile('[<>\(\)]')
    for cid in cids:
        recs = get_cid_website(cid)
        for rec in recs:
            if rec.get('TOCHeading') == 'Names and Identifiers':
                for subrec in rec['Section']:
                    if subrec.get('TOCHeading') == 'Record Title':
                        for info in subrec['Information']:
                            if info.get('Name') == 'Record Title':
                                names.append((cid, info['StringValue']))
                    elif subrec.get('TOCHeading') == 'Synonyms':
                        for subsubrec in subrec['Section']:
                            for info in subsubrec['Information']:
                                if info.get('Name') == 'MeSH Synonyms':
                                    for mesh_synonym in info['StringValueList']:
                                        mesh_synonyms.append((cid, mesh_synonym))
                                if info.get('Name') == 'Depositor-Supplied Synonyms':
                                    for syn in info['StringValueList']:
                                        depositor_supplied_synonyms.append((cid, syn))
            if rec.get('TOCHeading') == 'Pharmacology and Biochemistry':
                for subrec in rec['Section']:
                    if subrec.get('TOCHeading') == 'ATC Code':
                        for info in subrec['Information']:
                            if info.get('Name') == 'ATC Code':
                                atc = info['StringValue'].split()[0]
                                atcs.append((cid, atc))
            if rec.get('TOCHeading') == 'Biomolecular Interactions and Pathways':
                for subrec in rec['Section']:
                    if subrec.get('TOCHeading') == 'DrugBank Interactions':
                        for subsubrec in subrec['Section']:
                            if subsubrec.get('TOCHeading') in [
                                    'Target', 'Enzyme', 'Transporter', 'General Function']:
                                for info in subsubrec['Information']:
                                    if info.get('Name') in ['Target', 'Enzyme', 'Transporter']:
                                        name = ''.join(
                                            remove_links.split(info['StringValue'])[::2])
                                        bioentity = info['URL'].split('/')[-1]
                                        targets.append((cid, info['Name'], name, bioentity))
                                    elif info.get('Name') == 'General Function':
                                        g_function = ''.join(
                                            remove_links.split(info['StringValue'])[::2])
                                        general_functions.append((cid, g_function))

    return names, mesh_synonyms, depositor_supplied_synonyms, atcs, targets, general_functions


# %%
#if __name__ == '__main__':
#    data = [
#        ('CNC(=O)c1ccc(cc1F)N2C(=S)N(C(=O)C2(C)C)c3ccc(C#N)c(c3)C(F)(F)F;'
#         'CC1(C(=O)N(C(=S)N1c2ccc(c(c2)F)C(=O)NC)c3ccc(c(c3)C(F)(F)F)C#N)C'),
#        ('CN(C)C/C=C/C(=O)Nc1cc2c(cc1O[C@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F'),
#
#    ]
#    columns = ['SMILES or PubChem ID']
#    df = pd.DataFrame(data, columns=columns)
#    get_drug_features(df)
