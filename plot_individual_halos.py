#!/usr/bin/env python2.7

from __future__ import print_function, division

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk

from util import *


def print_min_max(data,label,slct):
    print(label,":", np.min(data[label][slct]), np.max(data[label][slct]))

def plot_halo(fof_tag, fof_data, sod_data, core_data, subhalo_data, bighalo_data, 
              show_plot = True, plot_bighalo=False, plot_subhalos=False):
    plt.figure()
    slct_cores = core_data['fof_tag'] == fof_tag
    slct_fof = fof_data['fof_tag'] == fof_tag
    slct_sod = sod_data['fof_tag'] == fof_tag

    if(np.sum(slct_fof) == 0):
        print("skip")
        return
    mass_fof = fof_data['fof_mass'][slct_fof][0]
    print("{:.2e}".format(mass_fof))
    x = fof_data['x'][slct_fof][0]
    y = fof_data['y'][slct_fof][0]
    z = fof_data['z'][slct_fof][0]
    mass_sod = sod_data['sod_halo_mass'][slct_sod][0]
    rad_sod  = sod_data['sod_halo_radius'][slct_sod][0]
    if plot_bighalo:
        slct_bighalo = bighalo_data['fof_tag']==fof_tag
        if(np.sum(slct_bighalo)==0):
            print('skip')
            return
        h,xbins,ybins = np.histogram2d(bighalo_data['x'][slct_bighalo], 
                                       bighalo_data['y'][slct_bighalo], 
                                       bins=100)
        plt.pcolor(xbins, ybins, h.T, cmap='Greys',norm=clr.LogNorm())
    sc = plt.scatter(core_data['x'][slct_cores], core_data['y'][slct_cores], 
                facecolors = 'none', c = core_data['infall_step'][slct_cores], label='Cores',
                s=80)
    
    plt.colorbar(sc)
    if plot_subhalos:
        slct_subhalo = subhalo_data['fof_tag'] == fof_tag    
        plt.scatter(subhalo_data['x'][slct_subhalo], subhalo_data['y'][slct_subhalo], 
                    facecolors = 'none', edgecolors = 'm', label='subhalos',
                    marker='^',s=90)

    sod_circle = plt.Circle((x,y),rad_sod, color='k', fill=False)
    ax = plt.gca()
    ax.add_artist(sod_circle)
    plt.plot(x,y,'xb',label='halo_center')
    plt.title("htag: {}\nM_fof: {:.2e}  M_sod: {:.2e}".format(fof_tag, mass_fof, mass_sod))

    plt.axes().set_aspect('equal','datalim')
    plt.legend(loc='best',framealpha=0.3)
    print_min_max(core_data,'radius',slct_cores)
    print_min_max(core_data,'infall_step',slct_cores)
    print_min_max(core_data,'infall_mass',slct_cores)
    print_min_max(core_data,'central',slct_cores)
    assert np.max(core_data['infall_step'][slct_cores]) == 499
    plt.grid()
    dtk.save_figs("figs/"+sys.argv[1]+"/"+__name__+"/")
    if show_plot:
        plt.show()
    plt.close()


def plot_individual_halos(param_file_name):
    param = dtk.Param(param_file_name)
    core_loc = param.get_string('core_loc')
    fof_loc = param.get_string('fof_loc')
    subhalo_loc = param.get_string('subhalo_loc')
    bighalo_loc = param.get_string('bighalo_loc')
    sod_loc = param.get_string('sod_loc')
    step = param.get_int('step')
    core_data_raw = load_cores(core_loc.replace('${step}',str(step)))
    fof_data = load_fof(fof_loc.replace('${step}',str(step)),1e14)
    sod_data = load_sod(sod_loc.replace('${step}',str(step)))
    #fof_data = combine_dics(fof_data, sod_data, 'fof_tag')
    core_data = core_data_raw


    plot_bighalo = True
    plot_subhalos = True

    if plot_subhalos:
        subhalo_data = load_subhalo(subhalo_loc.replace('${step}',str(step)))
    else:
        subhalo_data = None
    if plot_bighalo:
        bighalo_data = load_bighalo(bighalo_loc.replace('${step}',str(step)))
    else:
        bighalo_data = None
    #bighalo_data = None
    tags = np.unique(sod_data['fof_tag'])
    tags = np.unique(fof_data['fof_tag'])
    #tags = np.unique(bighalo_data['fof_tag'])
    for tag in tags:
        print(tag)
        plot_halo(tag, fof_data, sod_data, core_data, subhalo_data, bighalo_data, show_plot=True, plot_bighalo=plot_bighalo, plot_subhalos=plot_subhalos)

if __name__ == "__main__":
    plot_individual_halos(sys.argv[1])
