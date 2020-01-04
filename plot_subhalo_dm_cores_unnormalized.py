#!/usr/bin/env python3


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
    fof_mass_cut = param.get_float('fof_mass_cut')
    core_data_raw = load_cores(core_loc.replace('${step}',str(step)))
    fof_data = load_fof(fof_loc.replace('${step}',str(step)),fof_mass_cut)
    sod_data = load_sod(sod_loc.replace('${step}',str(step)))
    fof_data = combine_dics(fof_data, sod_data, 'fof_tag')
    core_data = get_diff(core_data_raw, fof_data, test=True)

    subhalo_data = get_diff(
        combine_dics( fof_data, 
                      load_subhalo( subhalo_loc.replace('${step}', str(step) )), 
                      'fof_tag'), 
        fof_data)
    subhalo_tmp = {'fof_tag': np.sort(np.unique(subhalo_data['fof_tag']))}
    fof_data = combine_dics(  fof_data, subhalo_tmp,'fof_tag')
    core_data = combine_dics(fof_data, core_data, 'fof_tag')
    fitted_core_data = get_diff(
        select_dic(core_data_raw, 
                   (core_data_raw['radius']<0.04) & (core_data_raw['central']==False)),
        fof_data)
    
    # bighalo_data = get_diff(load_bighalo(bighalo_loc.replace('${step}',str(step))), fof_data)
    plt.figure()
    plt.title("FoF Mass > {:.1e} and Has Subhalo Info".format(fof_mass_cut))
    r_bins = np.linspace(0,1.0,32)
    r_bins_avg = dtk.bins_avg(r_bins)
    r_bins_vol = 4.0/3.0*np.pi* ( r_bins[1:]**3 - r_bins[:-1]**3)
    subhalo_cnt = len(np.unique(subhalo_data['fof_tag']))
    core_cnt =  len(np.unique(core_data['fof_tag']))
    print("subhalo  num: ", subhalo_cnt)
    print("core  num: ", core_cnt)
    print("Unique halos in subhalos: ", np.shape(np.unique(subhalo_data['fof_tag'])))
    print("Unique halos in cores: ", np.shape(np.unique(core_data['fof_tag'])))
    print("min max subhalo halo mass: ", dtk.min_max(np.log10(subhalo_data['fof_mass'])))
    print("min max core halo mass: ", dtk.min_max(np.log10(core_data['fof_mass'])))
    # h, _ = np.histogram(bighalo_data['del_r_r200'], bins = r_bins)
    # # h=h/np.sum(h)
    # plt.plot(r_bins_avg, h, '-', label='Dark Matter')
    h, _ = np.histogram(subhalo_data['del_r_r200'], bins = r_bins)
    # h=h/np.sum(h)
    plt.plot(r_bins_avg, h, '-', label='Subhalo')

    h, _ = np.histogram(core_data['del_r_r200'], bins = r_bins)
    # h=h/np.sum(h)
    plt.plot(r_bins_avg, h, '-', label='All Cores')

    slct = core_data['central'] == False
    h, _ = np.histogram(core_data['del_r_r200'][slct], bins = r_bins)
    # h=h/np.sum(h)
    plt.plot(r_bins_avg, h, '-', label='Satellite Cores')

    h, _ = np.histogram(fitted_core_data['del_r_r200'], bins = r_bins)
    # h=h/np.sum(h)
    plt.plot(r_bins_avg, h, '-', label='Compact Satellites Cores')
    
    # h, _ = np.histogram(fitted_core_data['del_r_r200'], bins = r_bins)
    # h=h/np.sum(h)
    # plt.plot(r_bins_avg, h/r_bins_vol, '-o', label='cores fit')

    # plt.yscale('log')
    plt.grid()
    plt.ylabel('Count')
    plt.xlabel('R/R$_{200}$')
    plt.legend(loc='best')

    # plt.figure()
    # plt.hist(core_data['del_r_r200'])
    # plt.hist(subhalo_data['del_r_r200'])

    plt.figure()
    plt.plot(subhalo_data['fof_mass'], subhalo_data['subhalo_mass'], 'x', alpha=0.3)
    
    dtk.save_figs('figs/'+param_file_name+'/'+__file__+'/')
    plt.show()

def plot_cores_in_halos(param_file_name, plot_show = True):
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
    # fitted_core_data = get_diff(
    #     select_dic(core_data_raw, 
    #                (core_data_raw['radius']>0.04) & (core_data_raw['central']==False)),
    #     fof_data)
    subhalo_data = load_subhalo(subhalo_loc.replace('${step}',str(step)))
    print("subhalo_unique num: ", np.shape(np.unique(subhalo_data['fof_tag'])))
    subhalo_data = get_diff(subhalo_data, fof_data)
    bighalo_data = get_diff(load_bighalo(bighalo_loc.replace('${step}',str(step))), fof_data)
    fof_tags = np.unique(subhalo_data['fof_tag'])
    if plot_show:
        for fof_tag in fof_tags:
            print(fof_tag)
            slct_cores  = core_data_raw['fof_tag'] == fof_tag
            slct_subhalos = subhalo_data['fof_tag'] == fof_tag
            slct_fof   = fof_data['fof_tag'] == fof_tag
            slct_bh    = bighalo_data['fof_tag'] == fof_tag
            plt.figure()
            h,xbins,ybins = np.histogram2d(bighalo_data['x'][slct_bh], bighalo_data['y'][slct_bh],
                                            bins=256)
            plt.pcolor(xbins,ybins, h.T, cmap='Greys', norm=clr.LogNorm())
            plt.plot(fof_data['x'][slct_fof], fof_data['y'][slct_fof], 'x')
            plt.scatter(core_data_raw['x'][slct_cores], core_data_raw['y'][slct_cores], 
                        facecolor='None',
                        edgecolor='r',
                        label='Core')
            plt.scatter(subhalo_data['x'][slct_subhalos], subhalo_data['y'][slct_subhalos],  
                        facecolor='None',
                        edgecolor='b',
                        label='Subhalos')
            plt.axis('equal')
            plt.legend(loc='best')
            plt.show()

if __name__ == "__main__":
    plot_subhalo_dm_cores(sys.argv[1])
    plot_cores_in_halos(sys.argv[1])
