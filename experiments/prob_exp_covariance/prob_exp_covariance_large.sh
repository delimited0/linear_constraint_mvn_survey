#!/usr/bin/env bash

Rscript experiments/prob_exp_covariance/prob_exp_covariance.R \
--method_conf=experiments/prob_exp_covariance/method_conf.json \
--dim_conf=experiments/prob_exp_covariance/dim_large.json \
--dep=.1 \
--result_path=experiments/prob_exp_covariance/results \
--n_reps=18 \
--seed=2022 \
--n_cores=18 \
--blas_threads=1 & disown