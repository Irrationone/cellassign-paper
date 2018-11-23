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
import pickle

import matplotlib.pyplot as plt
import statsmodels.stats.stattools
import statsmodels.graphics.gofplots


# Plot QQplot
def plot_qq(res, logfcs, outfname, seed = 1481, variable = 'obs'):
    np.random.seed(seed)
    pred_samples = np.random.choice(res['ppc'][variable].flatten(), size = logfcs.shape[0])
    
    if variable == 'obs':
        xlab = "Quantiles of observed logFC"
    elif variable == 'null':
        xlab = "Quantiles of null logFC"
    else:
        raise Exception("No such variable option.")
    
    fig, ax = plt.subplots()
    fig = statsmodels.graphics.gofplots.qqplot_2samples(pred_samples, logfcs, line = '45',
                                                       xlabel = xlab,
                                                       ylabel = "Quantiles of posterior predictive samples")
    plt.savefig(outfname)
    
# Plot posterior predictive summaries
def plot_pp_summary(res, logfcs, stat = 'mean', nbins = 30):
    _, ax = plt.subplots(figsize=(12, 6))
    
    if stat == 'mean':
        # From pymc3 docs
        ax.hist([n.mean() for n in res['ppc']['obs']], bins=nbins, alpha=0.5)
        ax.axvline(logfcs.mean())
        ax.set(title='Posterior predictive of the mean', xlabel='mean(x)', ylabel='Frequency')
    else:
        raise Exception("Not implemented.")

# Plot logFC values
def plot_logfcs(logfcs, outfname):
    plt.hist(logfcs, bins = 'auto') 
    plt.title("logFC")
    plt.xlim((-0.02, 0.02))
    plt.savefig(outfname)
    
# Extract logFCs from input
def get_logfcs_comparison(df, type1, type2):
    type2 = re.sub("[ \+\\/]", ".", type2)
    return(np.array(df[df.cluster == type1]["logFC." + type2]))

# Read in arguments
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model_dir",
        dest="model_dir",
        help="Model result directory",
        metavar="DIR")
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
    
    model_file = os.path.join(args.model_dir, "model_results.pkl")
    with open(model_file, "rb") as input_file:
        model_result = pickle.load(input_file)
    
    # Write outputs
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    
    null_qq_plot = os.path.join(args.outdir, "qq_null.pdf")
    obs_qq_plot = os.path.join(args.outdir, "qq_obs.pdf")
    logfc_diff_plot = os.path.join(args.outdir, "logfc_diff_dist.pdf")
    logfc_same_plot = os.path.join(args.outdir, "logfc_same_dist.pdf")
    
    ## Create plots
    plot_qq(model_result, logfcs_same, outfname = null_qq_plot, variable = 'null')
    plot_qq(model_result, logfcs_diff, outfname = obs_qq_plot, variable = 'obs')
    plot_logfcs(logfcs_diff, logfc_diff_plot)
    plot_logfcs(logfcs_same, logfc_same_plot)

if __name__ == '__main__':
    args = read_args()
    
    main(args)
    
    