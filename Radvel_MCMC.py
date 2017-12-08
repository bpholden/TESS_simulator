from astropy.io import ascii
from astropy.table import Table
import numpy as np
import scipy as sp
import pandas as pd
import scipy.constants as sc
import scipy.optimize as op
import matplotlib.pyplot as plt
import re
import sys
import optparse
import copy
import os.path
from glob import glob
import re
sys.path.append("../simulator/")
import getpriority
import Generate_Velocities
from consts import JD0, ECC, OMEGA
import radvel
import ConfigParser



def readin_velsfile(fn):
    indata = ascii.read(fn)
    names=["time","mnvel","errvel","i2val","pphase"]
    for i in range(0,len(names)):
        oldn = "col%d" % (i+1)
        newn = names[i]
        indata.rename_column(oldn,newn)
    good = (indata['errvel'] > 0) & (indata['errvel'] < 100)
    return indata[good]

def initialize_model(planets,TESSAPFdata,addextra=False):
    time_base = 2458362.5

    nplan = len(planets)
    if addextra:
        nplan = nplan + 1
    params = radvel.Parameters(nplan,basis='per tp e w k') # number of planets = nplan
    n = 1
    for idx in planets:
        c_period = TESSAPFdata['period'][idx]
        k = radvel.utils.semi_amplitude(TESSAPFdata['est_mass'][idx],c_period,TESSAPFdata['mstar'][idx],0.)
        tpstart = time_base + c_period * (TESSAPFdata['phase'][idx] + 0.25)
        params['per%d' % (n)] = radvel.Parameter(value=c_period)
        params['tp%d' % (n)] = radvel.Parameter(value=tpstart)
        params['e%d' % (n)] = radvel.Parameter(value=0.0001)
        params['w%d' % (n)] = radvel.Parameter(value=np.pi/2.)
        params['k%d' % (n)] = radvel.Parameter(value=k)
        n = n + 1

    if addextra:
        params['per%d' % (n)] = radvel.Parameter(value=30.)
        params['tp%d' % (n)] = radvel.Parameter(value=time_base + 15.)
        params['e%d' % (n)] = radvel.Parameter(value=0.0001)
        params['w%d' % (n)] = radvel.Parameter(value=np.pi/2.)
        params['k%d' % (n)] = radvel.Parameter(value=2.)
        
    mod = radvel.RVModel(params, time_base=time_base)
    mod.params['dvdt'] = radvel.Parameter(value=0.0)
    mod.params['curv'] = radvel.Parameter(value=0.0)
    return mod


def make_like(mod,invels,planets,addextra=False):
    like = radvel.likelihood.RVLikelihood(mod, invels['time'], invels['mnvel'], invels['errvel'])
    like.params['gamma'] = radvel.Parameter(value=0.01)
    like.params['jit']= radvel.Parameter(value=3.)

    for n in range(1,like.params.num_planets+1):
        like.params['w%d' % (n)].vary = True
        like.params['e%d' % (n)].vary = False
        like.params['per%d' % (n)].vary = False
        like.params['tp%d' % (n)].vary = False

    if addextra:
        n = like.params.num_planets
        like.params['e%d' % (n)].vary = True
        like.params['per%d' % (n)].vary = True
        like.params['tp%d' % (n)].vary = True
        
    like.params['curv'].vary = False
    like.params['jit'].vary = True
    like.params['gamma'].vary = True        
    like.params['dvdt'].vary = False    
    
    return like

def init_posterior(like,planets,addextra=False):

    post = radvel.posterior.Posterior(like)
    post.priors += [radvel.prior.Gaussian( 'jit', np.log(3), 0.5)]
    post.priors += [radvel.prior.Gaussian( 'gamma', 0, 10)]
    n = 1
    ntot = len(planets)
    if addextra:
        ntot = ntot = 1
        
    for idx in range(0,ntot):
        n = idx + 1
        post.priors += [radvel.prior.Gaussian( 'k%d' % (n), np.log(5), 10)]
        
    return post

def plot_results(post,chains,outdir,sname):

    radvel.plotting.corner_plot(post, chains,saveplot=os.path.join(outdir,sname+"_corner.pdf"))
    radvel.plotting.rv_multipanel_plot(post,saveplot=os.path.join(outdir,sname+"_mp.pdf"))

    return

def writefit(sname,outdir,periods,Ks,err_Ks,Ms,err_Ms,rplanets,nvels):

    outname = "%s.tdradvelfit" % (sname)
    outname = os.path.join(outdir,outname)
    outfp = open(outname,"a+")
    for i in range(0,len(Ks)):
        ostr = "%f %f %f %f %f %f %d\n" % (periods[i],Ks[i],err_Ks[i],Ms[i],err_Ms[i],rplanets[i],nvels)
        # period (days) K (m/s) err_K (m/s) planet mass (M_earth) error (M_earth) R (R_earth) #vs 
        outfp.write(ostr)

    outfp.close()

def mcmc_planets(post,outdir,sname,mstars,addextra=False):
    
    conf_base = sname
    nwalkers = 20
    nsteps = 1000
    ensembles = 8
    msg = "Running MCMC for {}, N_ensembles = {}, N_walkers = {}, N_steps = {} ...".format(
        conf_base, ensembles, nwalkers, nsteps)
    print msg

    chains = radvel.mcmc(
        post, nwalkers=nwalkers, nrun=nsteps,
        ensembles=ensembles
    )

    Ks = []
    err_Ks = []
    Ms = []
    err_Ms = []
    # Get quantiles and update posterior object
    post_summary=chains.quantile([0.159, 0.5, 0.841])
    for n in range(1,post.params.num_planets+1):

        post_summary['Mpsini%d' % (n)] = radvel.utils.Msini(post_summary['k%d' % (n)],post.params['per%d' % (n)].value,mstars[0],0.)
        if writefit:
            Ks.append(post_summary['k%d' % (n)][0.5])
            err_Ks.append(( post_summary['k%d' % (n)][0.841] - post_summary['k%d' % (n)][0.159]) / 2.)
            Ms.append(post_summary['Mpsini%d' % (n)][0.5])
            err_Ms.append(( post_summary['Mpsini%d' % (n)][0.841] - post_summary['Mpsini%d' % (n)][0.159]) / 2.)

    if addextra:
        n = post.params.num_planets 
        post_summary['Mpsini%d' % (n)] = radvel.utils.Msini(post_summary['k%d' % (n)],post.params['per%d' % (n)].value,mstars[0],0.)
        if writefit:
            Ks.append(post_summary['k%d' % (n)][0.5])
            err_Ks.append(( post_summary['k%d' % (n)][0.841] - post_summary['k%d' % (n)][0.159]) / 2.)
            Ms.append(post_summary['Mpsini%d' % (n)][0.5])
            err_Ms.append(( post_summary['Mpsini%d' % (n)][0.841] - post_summary['Mpsini%d' % (n)][0.159]) / 2.)


    print "Saving output files..."
    saveto = os.path.join(outdir, sname+'_post_summary.csv')
    post_summary.to_csv(saveto, sep=',')
    
    return chains,Ks,err_Ks,Ms,err_Ms
