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
import string

import matplotlib.pyplot as plt
import statsmodels.stats.stattools
import statsmodels.graphics.gofplots


def combined_plot(splatter_res, v3_res, logfcs_diff, logfcs_same, outfname, seed = 1481, xmin, xmax):
    np.random.seed(seed)
    splatter_pred_samples = np.random.choice(splatter_res['ppc']['obs'].flatten(), size = logfcs_diff.shape[0])
    
    np.random.seed(seed)
    v3_pred_samples = np.random.choice(v3_res['ppc']['obs'].flatten(), size = logfcs_diff.shape[0])
    
    fig, axs = plt.subplots(2, 2, figsize=(12,10))
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.3, hspace=0.3)
    axs = axs.flat
    
    axs[0].hist(logfcs_same, bins = 'auto')
    axs[0].set_xlim((xmin, xmax))
    axs[0].set_xlabel("logFC")
    axs[0].set_ylabel("Frequency")
    
    axs[1].hist(logfcs_same, bins = 'auto')
    axs[1].set_xlim((xmin, xmax))
    axs[1].set_xlabel("logFC")
    axs[1].set_ylabel("Frequency")
    
    
    statsmodels.graphics.gofplots.qqplot_2samples(splatter_pred_samples, logfcs_diff, line = '45',
                                                       xlabel = "Quantiles of observed logFC",
                                                       ylabel = "Quantiles of posterior predictive samples",
                                                 ax = axs[2])
    
    statsmodels.graphics.gofplots.qqplot_2samples(v3_pred_samples, logfcs_diff, line = '45',
                                                       xlabel = "Quantiles of observed logFC",
                                                       ylabel = "Quantiles of posterior predictive samples",
                                                 ax = axs[3])
    
    for n, ax in enumerate(axs):
        ax.text(-0.1, 1.05, string.ascii_lowercase[n], transform=ax.transAxes, 
                size=20, weight='bold')
        
    plt.savefig(outfname)


# Extract logFCs from input
def get_logfcs_comparison(df, type1, type2):
    type2 = re.sub("[ \+\\/]", ".", type2)
    return(np.array(df[df.cluster == type1]["logFC." + type2]))

# Read in arguments
def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--splatter_model_dir",
        dest="splatter_model_dir",
        help="Splatter model result directory",
        metavar="DIR")
    parser.add_argument(
        "--v3_model_dir",
        dest="v3_model_dir",
        help="v3 model result directory",
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
        "--xmin",
        type=float,
        dest="xmin",
        help="Minimum logFC value to plot", 
        default = 2)
    parser.add_argument(
        "--xmax",
        type=float,
        dest="xmax",
        help="Maximum logFC value to plot", 
        default = 2)
    parser.add_argument(
        "-o", "--outdir", dest="outdir", help="Output directory", metavar="DIR")
    args = parser.parse_args()
    
    return(args)

def main(args):
    # Read inputs
    print("Reading logFCs ...")
    sys.stdout.flush()
    logfc_table = feather.read_dataframe(args.logfc_file)
    null_logfc_table = feather.read_dataframe(args.null_logfc_file)
    
    logfcs_diff = get_logfcs_comparison(logfc_table, args.class1, args.class2)
    logfcs_same = np.array(null_logfc_table[(null_logfc_table['celltype'] == args.class1) & (null_logfc_table['permutation'] <= 10)].logfc)
    
    print("Reading model files ...")
    sys.stdout.flush()
    splatter_model_file = os.path.join(args.splatter_model_dir, "model_results.pkl")
    with open(splatter_model_file, "rb") as input_file:
        splatter_model_result = pickle.load(input_file)
        
    v3_model_file = os.path.join(args.v3_model_dir, "model_results.pkl")
    with open(v3_model_file, "rb") as input_file:
        v3_model_result = pickle.load(input_file)
    
    # Write outputs
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
        
    combined_plot_file = os.path.join(args.outdir, "combined_plot.pdf")
    
    ## Create plots
    print("Creating combined plot ...")
    sys.stdout.flush()
    combined_plot(splatter_model_result, v3_model_result, 
                  logfcs_diff, logfcs_same, combined_plot_file, seed = 1481,
                  xmin = args.xmin, xmax=args.xmax)

if __name__ == '__main__':
    args = read_args()
    
    main(args)
    
    