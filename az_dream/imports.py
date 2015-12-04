# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:27:49 2015

@author: strokach
"""

import os
import os.path as op

import numpy as np
import scipy as sp
import pandas as pd
import sqlalchemy as sa

import matplotlib.pyplot as plt
import seaborn as sns

from IPython.display import display, HTML

from common import dat


#%% Load basic data
if False:
    drug_info_release = pd.read_csv('challenge_data/drug_synergy_data/Drug_info_release.csv', sep=',')

    ch1_train_combination_and_monotherapy = (
        pd.read_csv('challenge_data/drug_synergy_data/ch1_train_combination_and_monoTherapy.csv', sep=',')
    )

    ch1_leaderboard_monotherapy = (
        pd.read_csv('challenge_data/drug_synergy_data/ch1_leaderBoard_monoTherapy.csv', sep=',')
    )

    ch1_test_monotherapy = (
        pd.read_csv('challenge_data/drug_synergy_data/ch1_test_monoTherapy.csv', sep=',')
    )

    ch2_leaderboard_monotherapy = (
        pd.read_csv('challenge_data/drug_synergy_data/ch2_leaderBoard_monoTherapy.csv', sep=',')
    )

    ch2_test_monotherapy = (
        pd.read_csv('challenge_data/drug_synergy_data/ch2_test_monoTherapy.csv', sep=',')
    )


#%%
if False:
    store = pd.HDFStore('challenge_data/challenge_data.h5')
    store['drug_info_release'] = drug_info_release
    store['ch1_train_combination_and_monotherapy'] = ch1_train_combination_and_monotherapy
    store['ch1_leaderboard_monotherapy'] = ch1_leaderboard_monotherapy
    store['ch1_test_monotherapy'] = ch1_test_monotherapy
    store['ch2_leaderboard_monotherapy'] = ch2_leaderboard_monotherapy
    store['ch2_test_monotherapy'] = ch2_test_monotherapy
    store.close()
else:
    store = pd.HDFStore('challenge_data/challenge_data.h5')
    



