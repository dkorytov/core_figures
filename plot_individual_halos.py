#!/usr/bin/env python3

from __future__ import print_function, division

import sys
import os
import matplotlib
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk

from util import *


def print_min_max(data,label,slct):
    print(label,":", np.min(data[label][slct]), np.max(data[label][slct]))

def plot_halo(fof_tag, fof_data, sod_data, core_data, subhalo_data, bighalo_data, 
              show_plot = True, plot_bighalo=False, plot_subhalos=False, title=False):
    plt.figure(figsize=(7,5))
    slct_cores = (core_data['fof_tag'] == fof_tag) & (core_data['infall_mass'] > 2e11)

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
        h,xbins,ybins = np.histogram2d(bighalo_data['x'][slct_bighalo]-x, 
                                       bighalo_data['y'][slct_bighalo]-y, 
                                       bins=100)
        plt.pcolor(xbins, ybins, h.T, cmap='Greys',norm=clr.LogNorm())

    if plot_subhalos:
        slct_subhalo = (subhalo_data['fof_tag'] == fof_tag) & (subhalo_data['subhalo_mass'] > 1.6e11)
        plt.scatter(subhalo_data['x'][slct_subhalo]-x, subhalo_data['y'][slct_subhalo]-y, 
                    facecolors='none', edgecolors ='tab:blue', label='subhalos', lw=2,
                    marker='s', s=90)

    sc = plt.scatter(core_data['x'][slct_cores]-x, core_data['y'][slct_cores]-y, 
                     facecolors = 'tab:orange', edgecolor='none', label='cores',
                     s=30, alpha=0.66)
    
    # plt.colorbar(sc)

    sod_circle = plt.Circle((0,0),rad_sod, color='k', fill=False)
    ax = plt.gca()
    ax.add_artist(sod_circle)
    plt.plot([], [], 'k-', label='R$_{200c}$')
    plt.plot(x-x, y-y, 'x', color='k',label='halo center', lw=1, markersize=10)

    if title:
        plt.title("htag: {}\nM\_fof: {:.2e}  M\_sod: {:.2e}".format(fof_tag, mass_fof, mass_sod))

    plt.axes().set_aspect('equal','datalim')
    plt.legend(loc='upper right',framealpha=0.0)
    plt.ylabel('y [$h^{-1}$ Mpc]')
    plt.xlabel('x [$h^{-1}$ Mpc]')

    plt.text(0.975, 0.025, 'M$_{{200c}}$ = {:2.2e}[$h_{{-1}}$M$_{{200c}}$]'.format(mass_fof),
             transform=plt.gca().transAxes, verticalalignment='bottom',
             horizontalalignment='right')
    plt.tight_layout()
    # pyplot_zoom_out(right_amount=1.5, top_amount=0.5)
    pyplot_zoom_out(right_amount=0.5)
    print_min_max(core_data,'radius',slct_cores)
    print_min_max(core_data,'infall_step',slct_cores)
    print_min_max(core_data,'infall_mass',slct_cores)
    print_min_max(core_data,'central',slct_cores)
    assert np.max(core_data['infall_step'][slct_cores]) == 499
    # plt.grid()
    dtk.save_figs("figs/"+__file__+"/"+sys.argv[1]+"/")
    if show_plot:
        plt.show()
    plt.close()
    
def pyplot_zoom_out(amount=0.0, horizontal_amount=None, vertical_amount=None, right_amount=None, left_amount=None, top_amount=None, bottom_amount=None):
    """ increases the axes by fractional `amount`"""
    # Selects the first non-None value from the 'or' list
    right_amount = right_amount or horizontal_amount or amount
    left_amount = left_amount or horizontal_amount or amount
    top_amount = top_amount or vertical_amount or amount
    bottom_amount = bottom_amount or vertical_amount or amount
    
    ylim = plt.ylim()
    ydiff_cen = (ylim[1]-ylim[0])/2
    ycen = np.average(ylim)
    plt.ylim([ycen-ydiff_cen*(1+bottom_amount), ycen+ydiff_cen*(1.0+top_amount)])
    xlim = plt.xlim()
    xdiff_cen = (xlim[1]-xlim[0])/2
    xcen = np.average(xlim)
    plt.xlim([xcen-xdiff_cen*(1+left_amount), xcen+xdiff_cen*(1.0+right_amount)])


    
def plot_individual_halos(param_file_name):
    param = dtk.Param(param_file_name)
    core_loc = param.get_string('core_loc')
    fof_loc = param.get_string('fof_loc')
    if 'subhalo_loc' in param:
        subhalo_loc = param.get_string('subhalo_loc')
    else:
        subhalo_loc = None
        plot_subhalos = False
    bighalo_loc = param.get_string('bighalo_loc')
    sod_loc = param.get_string('sod_loc')
    step = param.get_int('step')
    core_data_raw = load_cores(core_loc.replace('${step}',str(step)))
    fof_data = load_fof(fof_loc.replace('${step}',str(step)),1e14)
    sod_data = load_sod(sod_loc.replace('${step}',str(step)))
    # fof_data = combine_dics(fof_data, sod_data, 'fof_tag')
    core_data = core_data_raw


    plot_bighalo = True
    plot_subhalos = True

    if plot_subhalos:
        subhalo_data = load_subhalo(subhalo_loc.replace('${step}',str(step)))
        print('subhalo data: ', subhalo_data)
    else:
        subhalo_data = None
    if plot_bighalo:
        bighalo_data = load_bighalo(bighalo_loc.replace('${step}',str(step)))
    else:
        bighalo_data = None
    print('subhalo data: ', subhalo_data)
    tags = np.unique(fof_data['fof_tag'])
    # tags = np.unique(bighalo_data['fof_tag'])
    for tag in tags:
        print(tag)
        print(subhalo_data)
        plot_halo(tag, fof_data, sod_data, core_data, subhalo_data, bighalo_data, show_plot=True, plot_bighalo=plot_bighalo, plot_subhalos=plot_subhalos)

if __name__ == "__main__":
    plot_individual_halos(sys.argv[1])
