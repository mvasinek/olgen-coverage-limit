from scipy.stats import binom, nbinom

import sys

INITIAL_COVERAGE = 30

def relaxTPFP(coverage, vaf, err, err_lim, det_lim):
    """Returns a tupple with minimal recomended coverage and variant allele frequency

    Parameters
    ----------
    coverage : int
        initial coverage to start with
    vaf : float
        variant allele frequency, assumes number from interval (0;1)
    err : float
        sequencing error
    err_lim : float
        cumulative probability limit for sequencing error
    det_lim : float
        cumulative probability limit for variant allele

    Returns
    -------
    tupple
        tupple containing coverage depth and variant allele frequency (absolute)
    """
    
    res_var = 0
    while True:
    
        #Relaxation condition no 1.
        while True:        
            errs = 0
            while binom.cdf(errs, coverage, err) < err_lim:
                errs += 1
            
            break

        #Relaxation condition no 2.
        variants = errs
        aux = 1 - binom.cdf(variants-1, coverage, vaf)
        
        if aux >= det_lim:
            res_var = variants
            break
        coverage += 1

    return (coverage,res_var)

def relaxTPFPMin(coverage, vaf, err, err_lim, det_lim, min_var):
    """Returns a tupple with minimal recomended coverage and variant allele frequency 
    with additional parameter minimal variant frequency(absolute)

    Parameters
    ----------
    coverage : int
        initial coverage to start with
    vaf : float
        variant allele frequency, assumes number from interval (0;1)
    err : float
        sequencing error
    err_lim : float
        cumulative probability limit for sequencing error
    det_lim : float
        cumulative probability limit for variant allele
    min_var : int
        minimal variant allele frequency

    Returns
    -------
    tupple
        tupple containing coverage depth and variant allele frequency (absolute)
    """
    res_var = 0
    while True:
        #Relaxation condition no 1.
        while True:        
            errs = 0
            while binom.cdf(errs, coverage, err) < err_lim:
                errs += 1
            
            break

        #Relaxation condition no 2.
        variants = errs
        aux = 1 - binom.cdf(variants-1, coverage, vaf)
        
        if aux >= det_lim and min_var <= variants:
            res_var = variants
            break
        coverage += 1

    return (coverage,res_var)
