
import tensorflow as tf

import az_dream.tf.settings as s

NUM_LABELS = 2 if s.CLF_TYPE == 'classifier' else 1

tf.set_random_seed(s.RANDOM_SEED)

graph = tf.Graph()
with graph.as_default():
    # Input data. For the training data, we use a placeholder that will be fed
    # at run time with a training minibatch.
    tf_train_dataset = tf.placeholder(
        tf.float32, shape=(s.N_ROWS, s.N_COLUMNS), name='tf_train_dataset')
    tf_train_labels = tf.placeholder(
        tf.float32, shape=(s.N_ROWS, NUM_LABELS), name='tf_train_labels')
    tf_step = tf.placeholder(tf.int64, name='tf_step')

    if s.tf_valid_dataset is not None:
        tf_valid_dataset = tf.constant(s.tf_valid_dataset, dtype=tf.float32)
    else:
        tf_valid_dataset = tf.placeholder(
            tf.float32, shape=(s.N_ROWS, s.N_COLUMNS), name='tf_valid_dataset')

    if s.tf_test_dataset is not None:
        tf_test_dataset = tf.constant(s.tf_test_dataset, dtype=tf.float32)
    else:
        tf_test_dataset = tf.placeholder(
            tf.float32, shape=(s.N_ROWS, s.N_COLUMNS), name='tf_test_dataset')

    # Variables.
    layer1_weights = tf.Variable(
        tf.truncated_normal(
            [s.N_COLUMNS, s.LAYER_1_SIZE],
            stddev=s.WEIGHT_STDDEV,
            seed=s.RANDOM_SEED),
        name='layer1_weights')
    layer1_biases = tf.Variable(tf.zeros([s.LAYER_1_SIZE]), name='layer1_biases')
    layer2_weights = tf.Variable(
        tf.truncated_normal(
            [s.LAYER_1_SIZE, NUM_LABELS],
            stddev=s.WEIGHT_STDDEV,
            seed=s.RANDOM_SEED),
        name='layer2_weights')
    layer2_biases = tf.Variable(tf.zeros([NUM_LABELS]), name='layer2_biases')

    def model(data, train=False):
        if train:
            data = tf.nn.dropout(data, s.DROPOUT, seed=s.RANDOM_SEED)
        hidden = tf.matmul(data, layer1_weights) + layer1_biases
        hidden = tf.nn.relu6(hidden)
        hidden = tf.matmul(hidden, layer2_weights) + layer2_biases
        return hidden

    def correlation(X, Y):
        X = X[:, 0:1]
        Y = Y[:, 0:1]
        X_mean = tf.reduce_mean(X)
        Y_mean = tf.reduce_mean(Y)
        X_std = tf.sqrt(tf.reduce_sum(tf.square(X - X_mean)))
        Y_std = tf.sqrt(tf.reduce_sum(tf.square(Y - Y_mean)))
        cov = tf.matmul(tf.transpose(X - X_mean), (Y - Y_mean))
        corr = cov / (X_std * Y_std)
        return tf.squeeze(corr)

    logits = model(tf_train_dataset, train=True)
    loss_c = tf.nn.softmax_cross_entropy_with_logits(logits, tf_train_labels)
    loss_r = -correlation(logits, tf_train_labels)

    if s.CLF_TYPE == 'classifier':
        loss = loss_c
    elif s.CLF_TYPE == 'regressor':
        loss = loss_r

    # Add regularization
    regularizers = (
        tf.nn.l2_loss(layer1_weights) + tf.nn.l2_loss(layer1_biases) +
        tf.nn.l2_loss(layer2_weights) + tf.nn.l2_loss(layer2_biases)
        # tf.nn.l2_loss(layer3_weights) + tf.nn.l2_loss(layer3_biases)
    )
    loss += s.LAMBDA * regularizers

    # Optimizer.
    # tf_step = tf.Variable(0, trainable=False)

    learning_rate = tf.train.exponential_decay(
        s.LEARNING_RATE, tf_step, 1, s.LEARNING_RATE_DECAY)
    optimizer = tf.train.AdagradOptimizer(learning_rate).minimize(loss)  # global_step=tf_step)

    # Predictions for the training, validation, and test data.
    train_prediction_c = tf.nn.softmax(logits)[:, 0]
    valid_prediction_c = tf.nn.softmax(model(tf_valid_dataset))[:, 0]
    test_prediction_c = tf.nn.softmax(model(tf_test_dataset))[:, 0]

    # Predictions for the training, validation, and test data.
    train_prediction_r = logits[:, 0]
    valid_prediction_r = model(tf_valid_dataset)[:, 0]
    test_prediction_r = model(tf_test_dataset)[:, 0]

    if s.WITH_TENSORBOARD:
        loss_r_sum = tf.scalar_summary('loss_r_summary', loss_r)
        loss_c_sum = tf.scalar_summary('loss_c_summary', tf.reduce_mean(loss_c))
        layer1_weights_hist = tf.histogram_summary("layer1_weights", layer1_weights)
        layer2_weights_hist = tf.histogram_summary("layer2_weights", layer2_weights)
        layer1_biases_hist = tf.histogram_summary("layer1_biases", layer1_biases)
        layer2_biases_hist = tf.histogram_summary("layer2_biases", layer2_biases)
        merged = tf.merge_all_summaries()
