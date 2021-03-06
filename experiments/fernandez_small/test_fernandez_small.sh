#!/usr/bin/env bash

Rscript experiments/fernandez_small/fernandez_small.R \
--method_conf=experiments/fernandez_small/method_conf.json \
--dim_conf=experiments/fernandez_small/dim_conf.json \
--result_path=experiments/fernandez_small/test_results \
--n_reps=4 \
--seed=2022 \
--n_cores=1 \
--n_blas_threads=4 & disown