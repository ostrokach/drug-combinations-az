import sys
import os.path as op
import time
from importlib import reload

import numpy as np
import scipy as sp
import pandas as pd
import sqlalchemy as sa
from sklearn.utils import shuffle, check_random_state
import tensorflow as tf

from az_dream import machine_learning


def main(
    train_df,
    valid_df,
    pairwise_data=None,
    #
    BATCH_SIZE=64,
    CLF_TYPE='classifier',
    FEATURE_FORMAT='pca',
    DO_SYN=True,
    DO_PCA=False,
    N_FEATURES=None,
    N_COMPONENTS=400,
    #
    LAYER_1_SIZE=32,
    WEIGHT_STDDEV=1,
    LAMBDA=1e-4,
    LEARNING_RATE=0.5,
    LEARNING_RATE_DECAY=0.9999999,
    #
    N_STEPS=203,
    STEP_OPTIONS=[1],
    SUBSAMPLE=False,
    DROPOUT=0.5,
    P_CUTOFF=0.5,
    #
    WITH_TENSORBOARD=False,
    RANDOM_STATE=42,
):
    """
    """
    random_state = check_random_state(RANDOM_STATE)

    ml = machine_learning.ML(
        clf_type=CLF_TYPE,
        feature_format=FEATURE_FORMAT,
        do_syn=DO_SYN,
        do_pca=DO_PCA,
        n_features=N_FEATURES,
        n_components=N_COMPONENTS,
        syn_cutoff=20.0,
        p_cutoff=P_CUTOFF,
        scale_factor=1,
    )

    # N_STEPS_PER_CHECK = N_STEPS // 5
    N_STEPS_PER_CHECK = 20

    train_info, train_dataset, train_labels, valid_info, valid_dataset, valid_labels = (
        ml.get_data_for_step(train_df, valid_df, pairwise_data)
    )

    train_info.index = range(len(train_info))
    valid_info.index = range(len(valid_info))

    grouped_dd_skipped = []
    grouped_dd = [
        g for g in train_info.groupby(['d_1', 'd_2'])
        if (len(g[1]) > 1 or grouped_dd_skipped.append(g))]
    print("grouped_dd_skipped: {}".format(len(grouped_dd_skipped)))

    grouped_c_skipped = []
    grouped_c = [
        g for g in train_info.groupby(['c'])
        if (len(g[1]) > 1 or grouped_c_skipped.append(g))]
    print("grouped_c_skipped: {}".format(len(grouped_c_skipped)))

    # TF graph parameters
    import az_dream.tf.settings
    reload(az_dream.tf.settings)
    import az_dream.tf.settings as s
    s.CLF_TYPE = CLF_TYPE
    s.N_ROWS = None
    s.N_COLUMNS = train_dataset.shape[1]

    s.LAYER_1_SIZE = LAYER_1_SIZE
    s.WEIGHT_STDDEV = WEIGHT_STDDEV
    s.LAMBDA = LAMBDA

    s.LEARNING_RATE = LEARNING_RATE
    s.LEARNING_RATE_DECAY = LEARNING_RATE_DECAY
    s.DROPOUT = DROPOUT
    s.WITH_TENSORBOARD = WITH_TENSORBOARD
    s.RANDOM_STATE = RANDOM_STATE

    print("Number of features: {}".format(s.N_COLUMNS))

    import az_dream.tf.graph
    reload(az_dream.tf.graph)
    import az_dream.tf.graph as g

    session_config = tf.ConfigProto(
        intra_op_parallelism_threads=1,
        inter_op_parallelism_threads=1,
        use_per_session_threads=True,
    )
    with tf.Session(graph=g.graph, config=session_config) as session:
        if WITH_TENSORBOARD:
            writer = tf.train.SummaryWriter("/tmp/az_dream", session.graph_def)
        tf.initialize_all_variables().run()
        print("Initialized\n")

        step_idxs = []
        ch1_scores_max = (-np.inf, -np.inf)
        ch2_scores_max = (-np.inf, -np.inf)
        for step in range(N_STEPS):
            batch_predictions_r_all = []
            batch_predictions_c_all = []
            batch_labels_all = []
            batch_losses = []

            def run_step(feed_dict):
                _, l, predictions_r, predictions_c = session.run(
                    [g.optimizer, g.loss, g.train_prediction_r, g.train_prediction_c],
                    feed_dict=feed_dict)
                batch_predictions_r_all.append(predictions_r)
                batch_predictions_c_all.append(predictions_c)
                batch_labels_all.append(batch_labels)
                batch_losses.append(l)

            if (step % N_STEPS_PER_CHECK) in STEP_OPTIONS:
                which = (step % N_STEPS_PER_CHECK)
            else:
                which = random_state.choice(STEP_OPTIONS)

            if which == 0:
                train_info_this_step, train_dataset_this_step, train_labels_this_step = (
                    shuffle(train_info, train_dataset, train_labels, random_state=random_state)
                )
                step_idxs.append(train_info_this_step.index[0])

                # Collect predictions from every batch
                for offset in range(0, train_dataset_this_step.shape[0] - BATCH_SIZE, BATCH_SIZE):
                    batch_data = train_dataset_this_step[offset:(offset + BATCH_SIZE), :]
                    batch_labels = train_labels_this_step[offset:(offset + BATCH_SIZE)]
                    feed_dict = {
                        g.tf_train_dataset: batch_data,
                        g.tf_train_labels: ml.format_labels(batch_labels),
                        g.tf_step: step,
                    }
                    run_step(feed_dict)
            # Group by drug pair
            elif which == 1:
                grouped = shuffle(grouped_dd, random_state=random_state)
                for key, group in grouped:
                    index = shuffle(group.index.tolist(), random_state=random_state)
                    step_idxs.append(index[0])

                    if SUBSAMPLE:
                        index_subset_size = np.random.randint(2, len(index) + 1)
                        index = index[:index_subset_size]
                    batch_data = train_dataset[index, :]
                    batch_labels = train_labels[index]
                    feed_dict = {
                        g.tf_train_dataset: batch_data,
                        g.tf_train_labels: ml.format_labels(batch_labels),
                        g.tf_step: step,
                    }
                    run_step(feed_dict)
            # Group by cell line
            elif which == 2:
                grouped = shuffle(grouped_c, random_state=random_state)
                for key, group in grouped:
                    index = shuffle(group.index.tolist(), random_state=random_state)
                    step_idxs.append(index[0])

                    if SUBSAMPLE:
                        index_subset_size = np.random.randint(2, len(index) + 1)
                        index = index[:index_subset_size]
                    batch_data = train_dataset[index, :]
                    batch_labels = train_labels[index]
                    feed_dict = {
                        g.tf_train_dataset: batch_data,
                        g.tf_train_labels: ml.format_labels(batch_labels),
                        g.tf_step: step,
                    }
                    run_step(feed_dict)

            step_loss = np.mean(np.hstack(batch_losses))
            step_labels = np.hstack(batch_labels_all)
            step_predictions_r = np.hstack(batch_predictions_r_all)
            step_predictions_c = np.hstack(batch_predictions_c_all)

            if step % N_STEPS_PER_CHECK in STEP_OPTIONS or (step == (N_STEPS - 1)):
                print("Step: {} ({})".format(step, which))
                print("Minibatch loss: {:.3}".format(step_loss))
                if pd.isnull(step_loss):
                    break
                # Correlation
                valid_predictions_r = g.valid_prediction_r.eval(
                    feed_dict={g.tf_valid_dataset: valid_dataset})
                print("Minibatch correlation: {:.3f} / {:.3f}".format(
                    sp.stats.pearsonr(step_labels, step_predictions_r)[0],
                    sp.stats.pearsonr(valid_labels, valid_predictions_r)[0],
                ))
                tmp_r = valid_info.copy()
                tmp_r['synergy_score'] = valid_labels
                tmp_r['synergy_score_pred'] = valid_predictions_r
                tmp_r = (
                    tmp_r
                    .groupby(['d_1', 'd_2', 'c'])
                    [['synergy_score', 'synergy_score_pred']]
                    .agg('mean')
                    .reset_index()
                )
                print('tmp_r shape: {}'.format(tmp_r.shape))
                ch1_scores = ml.get_score_ch1(tmp_r)
                if ch1_scores > ch1_scores_max:
                    ch1_scores_max = ch1_scores
                print("Ch1 scores: ({:.4f}, {:.4f})".format(*ch1_scores))
                # Accuracy
                valid_predictions_c = g.valid_prediction_c.eval(
                    feed_dict={g.tf_valid_dataset: valid_dataset})
                print("Minibatch accuracy: {:.3f} / {:.3f}".format(
                    ml.accuracy(step_labels, step_predictions_c),
                    ml.accuracy(valid_labels, valid_predictions_c),
                ))
                tmp_c = valid_info.copy()
                tmp_c['synergy_score'] = valid_labels
                tmp_c['synergy_score_pred'] = valid_predictions_c
                tmp_c = (
                    tmp_c
                    .groupby(['d_1', 'd_2', 'c'])
                    [['synergy_score', 'synergy_score_pred']]
                    .agg('mean').reset_index()
                )
                try:
                    ch2_scores = ml.get_score_ch2(tmp_c)
                except (np.linalg.LinAlgError, ZeroDivisionError, ValueError) as e:
                    print(e)
                    ch2_scores = (np.nan, np.nan)
                if ch2_scores > ch2_scores_max:
                    ch2_scores_max = ch2_scores

                print("Ch2 scores: ({:.2f}, {:.3f})".format(*ch2_scores))
                print('')
                if WITH_TENSORBOARD:
                    result_merged = g.merged.eval()
                    writer.add_summary(result_merged, step)
                sys.stdout.flush()

        # Checking that the random number generator is working...
        # print(step_idxs)

    return list(ch1_scores_max) + list(ch2_scores_max), tmp_r, tmp_c


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--subchallenge', type=str)

    parser.add_argument('--BATCH_SIZE', type=int, default=64)
    parser.add_argument('--FEATURE_FORMAT', type=str, default='scaled')
    parser.add_argument('--DO_SYN', type=int, default=1)
    parser.add_argument('--DO_PCA', type=int, default=0)
    parser.add_argument('--N_FEATURES', type=int, default=None)
    parser.add_argument('--N_COMPONENTS', type=int, default=400)

    parser.add_argument('--LAYER_1_SIZE', type=int, default=32)
    parser.add_argument('--WEIGHT_STDDEV', type=float, default=1)
    parser.add_argument('--LAMBDA', type=float, default=1e-4)
    parser.add_argument('--LEARNING_RATE', type=float, default=0.5)
    parser.add_argument('--LEARNING_RATE_DECAY', type=float, default=1)

    parser.add_argument('--N_STEPS', type=int, default=203)
    parser.add_argument('--STEP_OPTIONS', type=str, default='1')
    parser.add_argument('--SUBSAMPLE', type=int, default=0)
    parser.add_argument('--DROPOUT', type=float, default=0.5)
    parser.add_argument('--P_CUTOFF', type=float, default=0.3)
    parser.add_argument('--RANDOM_STATE', type=int, default=None)

    args = parser.parse_args()

    DATA_PATH = '/home/kimlab1/strokach/biodata/recipes/az_dream_2015/notebooks/machine_learning'
    if args.FEATURE_FORMAT == 'imputed':
        input_df = pd.read_hdf(
            op.join(DATA_PATH, 'ddc_data_{}_imputed_sqrt_allqa.h5'.format(args.subchallenge)))
    elif args.FEATURE_FORMAT == 'scaled':
        input_df = pd.read_hdf(
            op.join(DATA_PATH, 'ddc_data_{}_scaled_sqrt_allqa.h5'.format(args.subchallenge)))
    elif args.FEATURE_FORMAT == 'pca':
        input_df = pd.read_hdf(
            op.join(DATA_PATH, 'ddc_data_{}_pca_scaled_sqrt_allqa.h5'.format(args.subchallenge)))

    if args.subchallenge.startswith('ch1'):
        CLF_TYPE = 'regressor'
        df_train = (
            input_df[
                (input_df['source'] == 'train') |
                (input_df['source'] == 'ch2_validate')
            ]
        ).copy()

        df_validate = (
            input_df[
                (input_df['source'] == 'ch1_validate')
            ]
        ).copy()

        df_test = (
            input_df[
                (input_df['source'] == 'ch1_test')
            ]
        ).copy()
        pairwise_data = pd.read_hdf(op.join(DATA_PATH, 'pairwise_data.h5'))

    elif args.subchallenge.startswith('ch2'):
        CLF_TYPE = 'classifier'
        df_train = (
            input_df[
                (input_df['source'] == 'train') |
                (input_df['source'] == 'ch1_validate')
            ]
        ).copy()

        df_validate = (
            input_df[
                (input_df['source'] == 'ch2_validate')
            ]
        ).copy()

        df_test = (
            input_df[
                (input_df['source'] == 'ch2_test')
            ]
        ).copy()
        pairwise_data = None

    df_train.index = range(len(df_train))
    df_validate.index = range(len(df_validate))

    opts = args.__dict__.copy()
    opts['STEP_OPTIONS'] = [int(i) for i in opts['STEP_OPTIONS'].split(',')]
    opts['CLF_TYPE'] = CLF_TYPE
    opts.pop('subchallenge')

    scores_list = []
    if opts['RANDOM_STATE'] is None:
        print("Using default RANDOM_STATE...")
        for i in range(3):
            opts['RANDOM_STATE'] = i
            scores, tmp_r, tmp_c = main(df_train, df_validate, pairwise_data, **opts)
            scores_list.append(scores)
    else:
        scores, tmp_r, tmp_c = main(df_train, df_validate, pairwise_data, **opts)
        scores_list.append(scores)
    scores = [np.mean(x) for x in zip(*scores_list)]

    df = pd.DataFrame(dict(
        **args.__dict__,
        **dict(ch1_score1=scores[0], ch1_score2=scores[1],
               ch2_score1=scores[2], ch2_score2=scores[3])
    ), index=[0])

    n_tries = 0
    while n_tries < 5:
        try:
            df.to_sql(
                'tf_optimization_3',
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
