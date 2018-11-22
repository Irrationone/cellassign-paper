import os
import sys
import numpy as np
import pandas as pd
import pymc3 as pm
import scipy as sc
import theano
import theano.tensor as tt
import copy
import feather
import re
import argparse

import matplotlib.pyplot as plt
import statsmodels.stats.stattools
import statsmodels.graphics.gofplots

# Helper functions
def zvalue(value, sd=1, mu=0):
    """
    Calculate the z-value for a normal distribution. By default standard normal.
    """
    return (value - mu) / tt.sqrt(2. * sd ** 2.)

def cdf(mu, sd, value):
    z = zvalue(value, mu=mu, sd=sd)
    return tt.erfc(-z / tt.sqrt(2.))/2.

# Run null model
def run_null_model(nulls, iter_count = 2000, tune_iters = 2000):
    with pm.Model() as model:
        sd_null = pm.Gamma('sd_null', alpha = .1, beta = 1.)
        b_null = pm.Gamma('b_null', alpha = 1., beta = .1)
        
        dispersed_prob = pm.Beta('dispersed_prob', alpha = 1., beta = 1.)

        pm.Mixture('null',
                  comp_dists = [pm.Normal.dist(mu = 0., sd = sd_null),
                               pm.Laplace.dist(mu = 0., b = b_null)],
                  w = tt.as_tensor([1.-dispersed_prob, dispersed_prob]),
                  observed = nulls)

        pm.Deterministic('log_prob', model.logpt)

        trace = pm.sample(iter_count, tune = tune_iters, chains = 4)
        ppc = pm.sample_ppc(trace, samples=iter_count, model=model)
    
    return({'trace': trace, 'ppc': ppc})

# Splatter model with null model
def splatter_model(observations, nulls, null_sd, null_b, null_dispersed_prob, iter_count = 2000, tune_iters = 2000):
    with pm.Model() as model:
        # Probability of being a DE gene
        de_prob = pm.Uniform('de_prob', lower = 0., upper = 1.)

        # Probability of being downregulated
        down_prob = pm.Beta('down_prob', alpha = 1., beta = 1.)

        # Mean and sd for Gaussian for DE genes
        mu_pos = pm.Lognormal('mu_pos', mu = 0., sd = 1.)
        mu_neg = pm.Lognormal('mu_neg', mu = 0., sd = 1.)

        sd_pos = pm.Gamma('sd_pos', alpha = 1., beta = 1.)
        sd_neg = pm.Gamma('sd_neg', alpha = 1., beta = 1.)
        
        dispersed_prob = null_dispersed_prob
        spike_component = pm.Normal.dist(mu = 0., sd = null_sd) 
        slab_component = pm.Laplace.dist(mu = 0., b = null_b)

        # Sample from Gaussian-Laplace mixture for null (spike-and-slab mixture)
        pm.Mixture('null',
                  comp_dists = [spike_component,
                               slab_component],
                  w = tt.as_tensor([1.-dispersed_prob, dispersed_prob]),
                  observed = nulls)
        

        pos_component = pm.Bound(pm.Normal, lower=0.).dist(mu=mu_pos, sd=sd_pos)
        neg_component = pm.Bound(pm.Normal, upper=0.).dist(mu=-1*mu_neg, sd=sd_neg)
        pos_component_abs = pm.Bound(pm.Normal, lower=0.).dist(mu=-1*mu_pos, sd=sd_pos)
        neg_component_abs = pm.Bound(pm.Normal, upper=0.).dist(mu=mu_neg, sd=sd_neg)
        
        cdf_pos = cdf(mu = mu_pos, sd = sd_pos, value = 0.)
        cdf_neg = cdf(mu = -1*mu_neg, sd = sd_neg, value = 0.)

        pm.Mixture('obs',
                  w = tt.as_tensor([(1.-de_prob) * (1.-dispersed_prob),
                                    (1.-de_prob) * dispersed_prob, 
                                    de_prob * (1.-down_prob) * (1.-cdf_pos), 
                                    de_prob * down_prob * cdf_neg,
                                    de_prob * (1.-down_prob) * cdf_pos, 
                                    de_prob * down_prob * (1.-cdf_neg)]),
                  comp_dists = [spike_component,
                                slab_component, 
                                pos_component, 
                                neg_component, 
                                pos_component_abs,
                                neg_component_abs],
                  observed = observations)

        pm.Deterministic('log_prob', model.logpt)

        trace = pm.sample(iter_count, tune = tune_iters, chains = 4)
        ppc = pm.sample_ppc(trace, samples=iter_count, model=model)
    
    return({'trace': trace, 'ppc': ppc})

# V2 model
def v2_model(observations, nulls, null_sd, null_b, null_dispersed_prob, iter_count = 2000, tune_iters = 2000):
    with pm.Model() as model:
        # Probability of being a DE gene
        de_prob = pm.Beta('de_prob', alpha = 1., beta = 5.)
        
        # Probability of being downregulated
        down_prob = pm.Beta('down_prob', alpha = 1., beta = 1.)
        
        dispersed_prob = null_dispersed_prob

        mu_pos = pm.Lognormal('mu_pos', mu = -3, sd = 1.)
        mu_neg = pm.Lognormal('mu_neg', mu = -3, sd = 1.) 
        sd_pos = pm.Gamma('sd_pos', alpha = 0.01, beta = 1.) 
        sd_neg = pm.Gamma('sd_neg', alpha = 0.01, beta = 1.) 
        nu_pos = pm.Gamma('nu_pos', alpha = 5., beta = 1.)
        nu_neg = pm.Gamma('nu_neg', alpha = 5., beta = 1.)
        
        spike_component = pm.Normal.dist(mu = 0., sd = null_sd)
        slab_component = pm.Laplace.dist(mu = 0., b = null_b)

        # Sample from Gaussian-Laplace mixture for null (spike-and-slab mixture)
        pm.Mixture('null',
                  comp_dists = [spike_component,
                               slab_component],
                  w = tt.as_tensor([1.-dispersed_prob, dispersed_prob]),
                  observed = nulls)
    
        pos_component = pm.Bound(pm.StudentT, lower = 0.).dist(mu = mu_pos, sd = sd_pos, nu = nu_pos)
        neg_component = pm.Bound(pm.StudentT, upper = 0.).dist(mu = -mu_neg, sd = sd_neg, nu = nu_neg)
    
        pm.Mixture('obs',
                  w = tt.as_tensor([(1.-de_prob) * (1.-dispersed_prob),
                                    (1.-de_prob) * dispersed_prob,
                                    de_prob * (1.-down_prob),
                                   de_prob * down_prob]),
                  comp_dists = [spike_component, slab_component, pos_component, neg_component],
                  observed = observations)


        pm.Deterministic('log_prob', model.logpt)
        
        for RV in model.basic_RVs:
            print(RV.name, RV.logp(model.test_point))

        trace = pm.sample(iter_count, tune = tune_iters, chains = 4)
        ppc = pm.sample_ppc(trace, samples=iter_count, model=model)
    
    return({'trace': trace, 'ppc': ppc})

# V3 model
def v3_model(observations, nulls, null_sd, null_b, null_dispersed_prob, iter_count = 2000, tune_iters = 2000, max_fc = 1000):
    with pm.Model() as model:
        # Probability of being a DE gene
        de_prob = pm.Beta('de_prob', alpha = 1., beta = 5.)
        
        # Probability of being downregulated
        down_prob = pm.Beta('down_prob', alpha = 1., beta = 1.)
        
        dispersed_prob = null_dispersed_prob

        mu_pos = pm.Lognormal('mu_pos', mu = -3, sd = 1.)
        mu_neg = pm.Lognormal('mu_neg', mu = -3, sd = 1.) 
        sd_pos = pm.Gamma('sd_pos', alpha = 0.01, beta = 1.) 
        sd_neg = pm.Gamma('sd_neg', alpha = 0.01, beta = 1.) 
        nu_pos = pm.Gamma('nu_pos', alpha = 5., beta = 1.)
        nu_neg = pm.Gamma('nu_neg', alpha = 5., beta = 1.)
        
        spike_component = pm.Normal.dist(mu = 0., sd = null_sd)
        slab_component = pm.Laplace.dist(mu = 0., b = null_b)

        # Sample from Gaussian-Laplace mixture for null (spike-and-slab mixture)
        pm.Mixture('null',
                  comp_dists = [spike_component,
                               slab_component],
                  w = tt.as_tensor([1.-dispersed_prob, dispersed_prob]),
                  observed = nulls)
    
        pos_component = pm.Bound(pm.StudentT, lower = 0., upper = np.log(max_fc)).dist(mu = mu_pos, sd = sd_pos, nu = nu_pos)
        neg_component = pm.Bound(pm.StudentT, upper = 0., lower = -np.log(max_fc)).dist(mu = -mu_neg, sd = sd_neg, nu = nu_neg)
    
        pm.Mixture('obs',
                  w = tt.as_tensor([(1.-de_prob) * (1.-dispersed_prob),
                                    (1.-de_prob) * dispersed_prob,
                                    de_prob * (1.-down_prob),
                                   de_prob * down_prob]),
                  comp_dists = [spike_component, slab_component, pos_component, neg_component],
                  observed = observations)


        pm.Deterministic('log_prob', model.logpt)
        
        for RV in model.basic_RVs:
            print(RV.name, RV.logp(model.test_point))

        trace = pm.sample(iter_count, tune = tune_iters, chains = 4)
        ppc = pm.sample_ppc(trace, samples=iter_count, model=model)
    
    return({'trace': trace, 'ppc': ppc})

# Extract logFCs from input
def get_logfcs_comparison(df, type1, type2):
    type2 = re.sub("[ \+\\/]", ".", type2)
    return(np.array(df[df.cluster == type1]["logFC." + type2]))

# Plot QQplot
def plot_qq(res, logfcs):
    pred_samples = np.random.choice(res['ppc']['obs'].flatten(), size = logfcs.shape[0])
    statsmodels.graphics.gofplots.qqplot_2samples(pred_samples, logfcs, line = '45')

# Read in arguments
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--logfc_file",
        dest="logfc_file",
        help="File containing logFCs for each gene and class, i.e. merged output of findMarkers",
        metavar="FILE")
    parser.add_argument(
        "--null_logfc_file",
        dest="null_logfc_file",
        help="File containing permutations of null logFCs",
        metavar="FILE")
    parser.add_argument(
        "--class1",
        dest="class1",
        help="First class to compare")
    parser.add_argument(
        "--class2",
        dest="class2",
        help="Second class to compare")
    parser.add_argument(
        "-o", "--outdir", dest="outdir", help="Output directory", metavar="DIR")
    args = parser.parse_args()
    
    return(args)


def main(args):
    # Read inputs
    logfc_table = feather.read_dataframe(args.logfc_file)
    null_logfc_table = feather.read_dataframe(args.null_logfc_file)
    
    logfcs_diff = get_logfcs_comparison(logfc_table, args.class1, args.class2)
    logfcs_same = np.array(null_logfc_table[(null_logfc_table['celltype'] == args.class1) & (null_logfc_table['permutation'] <= 10)].logfc)
    
    # Run models
    null_model_results = run_null_model(logfcs_same, iter_count=4000, tune_iters = 2000)
    null_sd_mean = pm.summary(null_model_results['trace'])['mean']['sd_null']
    null_b_mean = pm.summary(null_model_results['trace'])['mean']['b_null']
    null_dispersed_prob_mean = pm.summary(null_model_results['trace'])['mean']['dispersed_prob']
    
    splatter_model_results = splatter_model(logfcs_diff, logfcs_same, null_sd = null_sd_mean, null_b = null_b_mean, null_dispersed_prob = null_dispersed_prob_mean, iter_count = 1000, tune_iters = 1000)
    
    v3_model_results = v3_model(logfcs_diff, logfcs_same, null_sd = null_sd_mean, null_b = null_b_mean, null_dispersed_prob = null_dispersed_prob_mean, iter_count = 2500, tune_iters = 1000)
    
    # Write outputs
    
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    
    splatter_model_trace_dir = os.path.join(args.outdir, "splatter")
    v3_model_trace_dir = os.path.join(args.outdir, "v3_model")
    
    ## Save traces
    pm.save_trace(splatter_model_results['trace'], directory = splatter_model_trace_dir)
    pm.save_trace(v3_model_results['trace'], directory = v3_model_trace_dir)
    
    ## Save ppc's
    splatter_model_results['ppc']
    
    ## Save logfcs used
    
    


if __name__ == '__main__':
    args = read_args()
    
    main(args)