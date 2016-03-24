import pandas as pd
# import pytest

from az_dream import functions as fn


# %% 1
# df = pd.read_csv('test_data/1.csv', index_col=0)
df = pd.read_csv('test_data/2.csv', index_col=0)

df_sg = fn.get_synergy_groups(df)

results_df = fn.loo_groupby(
    df, ['d_1', 'd_2'], features=['synergy_score'], stats=['mean'], name='xxx')

assert (df.reset_index()['index'] == results_df['index']).all().all()


#assert (
#    df_gp_1[-1:][['0_xxx_mean', '1_xxx_mean', '2_xxx_mean']].values ==
#    df[:-1].groupby(['d_1', 'd_2'])[['synergy_score']].agg('mean').values
#).all()

# assert (df[ic] == results_df[ic]).all().all()
