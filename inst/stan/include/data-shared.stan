  // model type
  int<lower=1,upper=4> family;
  // family == 1: bernoulli
  // family == 2: cumulative
  // family == 3: gaussian
  // family == 4: beta
  int<lower=1> N;  // total number of observations
  // response variable
  array[family == 1 || family == 2 ? N : 0] int Yint;
  array[family == 3 || family == 4 ? N : 0] real Yreal;
  int<lower=1> N_item;
  int<lower=1> N_itemC;  // item pairs
  int<lower=1> N_person;
  int<lower=1> N_trait;
  int<lower=0> N_item_fix;
  int<lower=0> N_item_est;
  // indices over N
  array[N] int<lower=1> J_item1;
  array[N] int<lower=1> J_item2;
  array[N] int<lower=1> J_itemC;
  array[N] int<lower=1> J_person;
  array[N] int<lower=1> J_trait1;
  array[N] int<lower=1> J_trait2;
  array[N_item_fix] int<lower=1> J_item_fix;
  array[N_item_est] int<lower=1> J_item_est;
  // indicate inverted items
  int<lower=0> N_item_pos;
  int<lower=0> N_item_neg;
  array[N_item_pos] int<lower=1> J_item_pos;
  array[N_item_neg] int<lower=1> J_item_neg;
  // indicate items used in multiple blocks
  int<lower=0> N_item_equal;
  array[N_item_equal] int<lower=1> J_item_equal;
  array[N_item_equal] int<lower=1> J_item_orig;
  // number of response categories in ordinal models
  // should be set to 2 for other families
  int<lower=2> ncat;
