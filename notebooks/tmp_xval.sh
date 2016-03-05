#!/bin/bash

for loss in 'ls' 'huber' 'lad' ; do
  for learning_rate in 0.01 0.05 0.1 0.2 ; do  #  ...
    for n_estimators in 20 50 100 200 500 1000 2000 3000 4000 ; do  #  ...
      for subsample in 0.1 0.2 0.3 0.5 0.7 0.9 1 ; do
        echo $learning_rate $n_estimators $subsample
        for max_features in 0.1 0.3 0.5 0.7 0.8 0.9 1 ; do
          for max_depth in 1 2 3 4 ; do
            for weight_nonsyn in 0.1 0.3 0.5 0.7 0.9 1 ; do
              status=1
              while [[ $status != 0 ]]; do
                submitjob 1 -c 1 -m 2 python tmp_xval.py --loss $loss --learning_rate $learning_rate --n_estimators $n_estimators --subsample $subsample --max_features $max_features --max_depth $max_depth --weight_nonsyn $weight_nonsyn
                status=$?
                echo 'status:' $status
                if [[ $status != 0 ]] ; then
                  sleep 300
                fi
              done
              sleep 0.5
            done
          done
        done
      done
    done
  done
done
