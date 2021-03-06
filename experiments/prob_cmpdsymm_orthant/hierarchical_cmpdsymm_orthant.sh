#!/usr/bin/env bash

Rscript experiments/prob_cmpdsymm_orthant/prob_cmpdsymm_orthant.R \
--method_conf=experiments/prob_cmpdsymm_orthant/hierarchical_conf.json \
--dim_conf=experiments/prob_cmpdsymm_orthant/hierarchical_dim.json \
--corr=.5 \
--result_path=experiments/prob_cmpdsymm_orthant/results \
--n_reps=18 \
--seed=2022 \
--n_cores=18