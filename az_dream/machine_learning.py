import numpy as np
import scipy as sp
import pandas as pd

from sklearn.preprocessing import Imputer, StandardScaler, RobustScaler, MaxAbsScaler
from sklearn.decomposition import PCA

import tensorflow as tf

from az_dream import functions as fn


class ML:

    def __init__(
            self,
            clf_type='regressor',
            feature_format='pca',
            pairwise=False,
            do_syn=False,
            do_pca=False,
            n_features=None,
            n_components=300,
            syn_cutoff=20.0,
            p_cutoff=0.5,
            scale_factor=0.01):
        """
        """
        self.clf_type = clf_type
        self.feature_format = feature_format
        self.pairwise = pairwise
        self.do_syn = do_syn
        self.do_pca = do_pca
        self.n_features = n_features
        self.n_components = n_components
        self.syn_cutoff = syn_cutoff
        self.p_cutoff = p_cutoff
        self.scale_factor = scale_factor

    features_to_exclude = [
        'synergy_score', 'synergy_score_y', 'synergy_score_diff',
        'd_1', 'd_2', 'c', 'qa', 'source', 'unique_id',
        'index'
    ]

    def get_feature_columns(self, df):
        feature_columns = [
            c for c in df.columns
            if (df[c].dtype in (int, float) and c not in self.features_to_exclude) or print(c)
        ]
        return feature_columns

    def get_this_feature_columns(self, df):
        print("Input type: {}".format(self.feature_format))
        if self.feature_format == 'pca':
            this_feature_columns = list(range(self.n_components))
        else:
            this_feature_columns = self.get_feature_columns(df)
        return this_feature_columns

    def pca(self, train_df, test_df, fc):
        pca = PCA()
        standard_scaler = StandardScaler()
        maxabs_scaler = MaxAbsScaler()

        pca.fit(np.vstack([train_df[fc].values, test_df[fc].values]))
        train_df.loc[:, fc] = pca.transform(train_df[fc].values)
        test_df.loc[:, fc] = pca.transform(test_df[fc].values)

        standard_scaler.fit(np.vstack([train_df[fc].values, test_df[fc].values]))
        train_df.loc[:, fc] = standard_scaler.transform(train_df[fc].values)
        test_df.loc[:, fc] = standard_scaler.transform(test_df[fc].values)

        maxabs_scaler.fit(np.vstack([train_df[fc].values, test_df[fc].values]))
        train_df.loc[:, fc] = maxabs_scaler.transform(train_df[fc].values)
        test_df.loc[:, fc] = maxabs_scaler.transform(test_df[fc].values)

        train_dataset = train_df[fc].values[:, :self.n_components]
        valid_dataset = test_df[fc].values[:, :self.n_components]

        return train_dataset, valid_dataset

    def only_meaningful_features(self, train_dataset, valid_dataset):
        from sklearn.feature_selection import VarianceThreshold
        selector = VarianceThreshold()
        selector.fit(train_dataset)
        train_dataset = selector.transform(train_dataset)
        valid_dataset = selector.transform(valid_dataset)
        return train_dataset, valid_dataset

    def select_features(self, train_dataset, valid_dataset, train_labels):
        from sklearn.feature_selection import SelectKBest, f_regression
        selector = SelectKBest(f_regression, k=self.n_features)
        selector.fit(train_dataset, train_labels)
        train_dataset = selector.transform(train_dataset)
        valid_dataset = selector.transform(valid_dataset)
        return train_dataset, valid_dataset

    def get_data_for_step(self, train_df, test_df, pairwise_data=None):
        imputer = Imputer(strategy='median')
        robust_scaler = RobustScaler()
        maxabs_scaler = MaxAbsScaler()

        fc = self.get_this_feature_columns(train_df)

        if self.do_syn:
#            train_sg = fn.get_synergy_groups(
#                train_df, train_df.drop('synergy_score', axis=1), pairwise_data)
#            test_sg = fn.get_synergy_groups(
#                train_df, test_df.drop('synergy_score', axis=1), pairwise_data)
            train_sg = fn.get_synergy_groups(train_df, train_df, pairwise_data)
            test_sg = fn.get_synergy_groups(train_df, test_df, pairwise_data)

            if pairwise_data is not None:
                nfc = [
                    c for c in train_sg.columns
                    if c not in fc and
                    train_sg[c].dtype in (int, float) and
                    c not in ['synergy_score', 'synergy_score_diff']
                ]
            else:
                nfc = [
                    c for c in train_sg.columns
                    if any([isinstance(c, str) and c.startswith(prefix)
                           for prefix in ['count_', 'synergy_score_']])
                    # for prefix in ['synergy_score_']])
                ]
            print(nfc)

            assert not any([c in (fc + nfc) for c in ['synergy_score', 'synergy_score_diff']])

            imputer.fit(
                np.vstack([train_sg[nfc].values, test_sg[nfc].values]))
            train_sg.loc[:, nfc] = imputer.transform(train_sg[nfc].values)
            test_sg.loc[:, nfc] = imputer.transform(test_sg[nfc].values)

            robust_scaler.fit(
                np.vstack([train_sg[nfc].values, test_sg[nfc].values]))
            train_sg.loc[:, nfc] = robust_scaler.transform(train_sg[nfc].values)
            test_sg.loc[:, nfc] = robust_scaler.transform(test_sg[nfc].values)

            maxabs_scaler.fit(
                np.vstack([train_sg[nfc].values, test_sg[nfc].values]))
            train_sg.loc[:, nfc] = maxabs_scaler.transform(train_sg[nfc].values)
            test_sg.loc[:, nfc] = maxabs_scaler.transform(test_sg[nfc].values)

            fc = [c for c in fc if c in train_sg.columns]
            if self.do_pca:
                train_dataset, valid_dataset = self.pca(train_sg, test_sg, fc + nfc)
            else:
                train_dataset = train_sg[fc + nfc].values
                valid_dataset = test_sg[fc + nfc].values
                train_dataset, valid_dataset = self.only_meaningful_features(train_dataset, valid_dataset)

            train_info = train_sg[['d_1', 'd_2', 'c']]
            valid_info = test_sg[['d_1', 'd_2', 'c']]

            train_labels = train_sg['synergy_score'].values * self.scale_factor
            valid_labels = test_sg['synergy_score'].values * self.scale_factor

            train_sg.drop('synergy_score', axis=1, inplace=True)
            test_sg.drop('synergy_score', axis=1, inplace=True)
            # assert (train_df[['d_1', 'd_2', 'c']] == train_sg[['d_1', 'd_2', 'c']]).all().all()
            # assert (test_df[['d_1', 'd_2', 'c']] == test_sg[['d_1', 'd_2', 'c']]).all().all()

        else:
            if self.do_pca:
                train_dataset, valid_dataset = self.pca(train_df, test_df, fc)
            else:
                train_dataset = train_df[fc].values
                valid_dataset = test_df[fc].values

            train_info = train_df[['d_1', 'd_2', 'c']]
            valid_info = test_df[['d_1', 'd_2', 'c']]

            train_labels = train_df['synergy_score'].values * self.scale_factor
            valid_labels = test_df['synergy_score'].values * self.scale_factor

        if self.n_features:
            train_dataset, valid_dataset = (
                self.select_features(train_dataset, valid_dataset, train_labels)
            )

        return train_info, train_dataset, train_labels, valid_info, valid_dataset, valid_labels

    def accuracy(self, y, y_pred):
        # acc = (y == (y_pred > 20)).sum() / y.shape[0]
        acc = ((y > 0.2) == (y_pred > self.p_cutoff)).sum() / y.shape[0]
        return acc

    def format_labels(self, labels):
        if self.clf_type == 'classifier':
            classes = np.stack([
                (labels > (self.syn_cutoff * self.scale_factor)).astype(int),
                (labels <= (self.syn_cutoff * self.scale_factor)).astype(int)
            ]).T
        elif self.clf_type == 'regressor':
            classes = np.reshape(labels, (-1, 1))
        # classes = (labels > (20.0 * scale_factor)).astype(int)
        return classes

    def get_score_ch1(self, df):
        """
        There is a slight deviation from the R script, but not sure where it comes from...
        It is in the nomenator of the

            DrugCombiScore_ch1_Wsq <- sum(R$V1*sqrt(R$N-1))/sum(sqrt(R$N-1))  # primary metric

        fraction, and at least some of the `R$V1` values are the same...
        """
        by = ['d_1', 'd_2']
        columns = ['d_1', 'd_2', 'c', 'synergy_score', 'synergy_score_pred']

        results = []
        grouped = df[columns].groupby(by)
        n_skipped = 0
        for key, group in grouped:
            if group.shape[0] <= 1:
                # results.append(list(key) + [0, 0, 0])
                # print("Skipping key '{}' because it has only 1 value...".format(key))
                n_skipped += 1
                continue
            else:
                count = group['synergy_score'].count()
                ss_max = group['synergy_score'].max()
                P_r, P_p = sp.stats.pearsonr(group['synergy_score'], group['synergy_score_pred'])
                results.append(list(key) + [count, ss_max, P_r])

        # print("Skipped {} out of {} groups...".format(n_skipped, len(grouped)))

        results_df = pd.DataFrame(
            results, columns=by + ['count', 'synergy_score_max', 'pearson_r'])

        primary_metric = (
            np.sum(results_df['pearson_r'] * np.sqrt(results_df['count'].values - 1)) /
            np.sum(np.sqrt(results_df['count'].values - 1))
        )

        results_df_high = results_df[results_df['synergy_score_max'] > self.syn_cutoff]
        tiebreak_metric = (
            np.sum(results_df_high['pearson_r'] * np.sqrt(results_df_high['count'].values - 1)) /
            np.sum(np.sqrt(results_df_high['count'].values - 1))
        )

        return primary_metric, tiebreak_metric

    def get_score_ch2(self, df, probas=True):
        from statsmodels.formula.api import ols
        from statsmodels.stats.anova import anova_lm

        df = df.copy()
        df['dd'] = df['d_1'] + '.' + df['d_2']
        if probas:
            df['synergy_score_pred'] = (df['synergy_score_pred'] > self.p_cutoff).astype(int)
        nosum_lm = ols('synergy_score ~ dd + c + synergy_score_pred', data=df).fit()
        pval = -np.log10(anova_lm(nosum_lm).loc['synergy_score_pred', 'PR(>F)'])
        # Positive or negative?
        if (df[df['synergy_score_pred'] == 1]['synergy_score'].mean() >=
                df[df['synergy_score_pred'] == 0]['synergy_score'].mean()):
            sign = 1
        else:
            sign = -1
        pval = sign * pval
        # BAC
        df['synergy_score'] = (
            # (df['synergy_score'] > (self.syn_cutoff * self.scale_factor))
            (df['synergy_score'] > (self.syn_cutoff))
        ).astype(int)
        pos = df[df['synergy_score'] == 1]
        neg = df[df['synergy_score'] == 0]
        # Sensitivity (true positive rate)
        tpr = (pos['synergy_score_pred'] == 1).sum() / pos.shape[0]
        # Specificity (true negative rate)
        tnr = (neg['synergy_score_pred'] == 0).sum() / neg.shape[0]
        # print(tpr, tnr)
        bac = (tpr + tnr) / 2
        return pval, bac
