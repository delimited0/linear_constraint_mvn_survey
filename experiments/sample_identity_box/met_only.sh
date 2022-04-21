#!/usr/bin/env bash

Rscript experiments/sample_identity_box/identity_box.R \
--method_conf=experiments/sample_identity_box/met_only.json \
--dim_conf=experiments/sample_identity_box/dim_conf.json \
--half_width=3 \
--sample_path=experiments/sample_identity_box/samples \
--seed=2022 \
--n_threads=18