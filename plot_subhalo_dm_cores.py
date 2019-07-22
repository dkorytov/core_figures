#!/usr/bin/env python2.7


from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk
import sys

from util import *

def plot_diff(data,bins,label):
    h,xbins = np.histogram(data,bins=bins)
    plt.plot(dtk.bins_avg(xbins), h, label = label)

def plot_subhalo_dm_cores(param_file_name, plot_show = True):
    param = dtk.Param(sys.argv[1])
    core_loc = param.get_string('core_loc')
    fof_loc = param.get_string('fof_loc')
    subhalo_loc = param.get_string('subhalo_loc')
    bighalo_loc = param.get_string('bighalo_loc')
    sod_loc = param.get_string('sod_loc')
    step = param.get_int('step')
    core_data_raw = load_cores(core_loc.replace('${step}',str(step)))
    fof_data = load_fof(fof_loc.replace('${step}',str(step)),1e14)
    sod_data = load_sod(sod_loc.replace('${step}',str(step)))
    fof_data = combine_dics(fof_data, sod_data, 'fof_tag')
    core_data = get_diff(core_data_raw, fof_data)
    fitted_core_data = get_diff(
        select_dic(core_data_raw, 
                   (core_data_raw['radius']>0.04) & (core_data_raw['infall_mass']>10**12.08)),
        fof_data)
    subhalo_data = get_diff(load_subhalo(subhalo_loc.replace('${step}',str(step))), fof_data)
    bighalo_data = get_diff(load_bighalo(bighalo_loc.replace('${step}',str(step))), fof_data)
    plt.figure()
    r_bins = np.linspace(0,1.0,16)
    r_bins_avg = dtk.bins_avg(r_bins)
    r_bins_vol = 4.0/3.0*np.pi* ( r_bins[1:]**3 - r_bins[:-1]**3)

    h, _ = np.histogram(bighalo_data['del_r_r200'], bins = r_bins)
    h=h/np.sum(h)
    plt.plot(r_bins_avg, h/r_bins_vol, '-', label='Dark Matter')
    h, _ = np.histogram(subhalo_data['del_r_r200'], bins = r_bins)
    h=h/np.sum(h)
    plt.plot(r_bins_avg, h/r_bins_vol, '-', label='Subhalo')
    h, _ = np.histogram(core_data['del_r_r200'], bins = r_bins)
    h=h/np.sum(h)
    plt.plot(r_bins_avg, h/r_bins_vol, '-', label='All Cores')
    # h, _ = np.histogram(fitted_core_data['del_r_r200'], bins = r_bins)
    # h=h/np.sum(h)
    # plt.plot(r_bins_avg, h/r_bins_vol, '-o', label='cores fit')
    
    plt.yscale('log')
    plt.grid()
    plt.ylabel('Density [h^3 Mpc^-3]')
    plt.xlabel('R/R$_{200}$')
    plt.legend(loc='best')

    plt.figure()
    plt.hist(core_data['del_r_r200'])
    plt.hist(subhalo_data['del_r_r200'])

    dtk.save_figs('figs/'+param_file_name+'/'+__file__+'/')
    plt.show()

def plot_cores_in_halos(param_file_name,plot_show = True):
    param = dtk.Param(sys.argv[1])
    core_loc = param.get_string('core_loc')
    fof_loc = param.get_string('fof_loc')
    subhalo_loc = param.get_string('subhalo_loc')
    bighalo_loc = param.get_string('bighalo_loc')
    sod_loc = param.get_string('sod_loc')
    step = param.get_int('step')
    core_data_raw = load_cores(core_loc.replace('${step}',str(step)))
    fof_data = load_fof(fof_loc.replace('${step}',str(step)),1e14)
    sod_data = load_sod(sod_loc.replace('${step}',str(step)))
    fof_data = combine_dics(fof_data, sod_data, 'fof_tag')
    core_data = get_diff(core_data_raw, fof_data)
    fitted_core_data = get_diff(
        select_dic(core_data_raw, 
                   (core_data_raw['radius']>0.04) & (core_data_raw['infall_mass']>10**12.08)),
        fof_data)
    subhalo_data = get_diff(load_subhalo(subhalo_loc.replace('${step}',str(step))), fof_data)
    bighalo_data = get_diff(load_bighalo(bighalo_loc.replace('${step}',str(step))), fof_data)
    
    if plot_show:
        plt.show()
if __name__ == "__main__":
    plot_subhalo_dm_cores(sys.argv[1])
