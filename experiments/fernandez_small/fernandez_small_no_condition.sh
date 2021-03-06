#!/usr/bin/env bash

Rscript experiments/fernandez_small/fernandez_small.R \
--method_conf=experiments/fernandez_small/no_condition.json \
--dim_conf=experiments/fernandez_small/dim_conf.json \
--result_path=experiments/fernandez_small/results \
--n_reps=18 \
--seed=2022 \
--n_cores=18 \
--n_blas_threads=1 & disown