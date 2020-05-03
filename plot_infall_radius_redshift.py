#!/usr/bin/env python3

from __future__ import division, print_function

from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk
import sys

def sort_dic(dic, label):
    srt = np.argsort(dic[label])
    for key in dic.keys():
        dic[key] = dic[key][srt]
    return


def select_dic(dic,slct):
    new_dic = {}
    for key in dic.keys():
        new_dic[key] = dic[key][slct]
    return new_dic


def left_join_dic(dic_a,dic_b,key):
    sort_dic(dic_b,key)
    indx = dtk.search_sorted(dic_b[key],dic_a[key])
    return left_combine_dic(dic_a,select_dic(dic_b,indx))


def left_combine_dic(dic_a, dic_b):
    for key in dic_b.keys():
        dic_a[key]=dic_b[key]
    return dic_a


def combine_dic(dic_a,dic_b):
    new_dic = {}
    for key in dic_a.keys():
        new_dic[key] = dic_a[key]
    for key in dic_b.keys():
        new_dic[key] = dic_b[key]
    return dic


def load_cores(core_loc,redshift, stepz):
    result = {}
    dtk.gio_inspect(core_loc)
    print("loading in...")
    print('fof_tag')
    result['fof_tag'] = dtk.gio_read(core_loc,'fof_halo_tag')
    print('x')
    result['x']       = dtk.gio_read(core_loc,'x')
    print('y')
    result['y']       = dtk.gio_read(core_loc,'y')
    print('z')
    result['z']       = dtk.gio_read(core_loc,'z')
    print('infall_mass')
    result['infall_mass'] = dtk.gio_read(core_loc,'infall_tree_node_mass')
    print('radius')
    result['radius'] = dtk.gio_read(core_loc,'radius')
    print('infall_step')
    result['infall_step'] = dtk.gio_read(core_loc,'infall_step')
    print('central')
    result['central'] = dtk.gio_read(core_loc,'central')
    print('infall_redshift')
    result['infall_redshift'] = stepz.get_z(result['infall_step'])
    result['infall_redshift+1'] = result['infall_redshift']+0.0999
    return result


def load_fof(fof_loc):
    result = {}
    dtk.gio_inspect(fof_loc)
    print('loading fof tags')
    result['fof_tag'] = dtk.gio_read(fof_loc,'fof_halo_tag')
    print('loading fof mass')
    result['fof_mass'] = dtk.gio_read(fof_loc,'fof_halo_mass')
    return result


def process_cores(cores,cores_cen):
    #TODO finish
    indx = dtk.search_sorted(cores_cen['fof_tag'], cores['fof_tag'])
    #assert(np.sum(indx==-1) == 0)
    cores['cen_dx'] = cores['x']-cores_cen['x'][indx]
    cores['cen_dy'] = cores['y']-cores_cen['y'][indx]
    cores['cen_dz'] = cores['z']-cores_cen['z'][indx]
    cores['cen_dr'] = np.sqrt(cores['cen_dx']**2 + cores['cen_dy']**2 + cores['cen_dz']**2 )
    u,cnt = np.unique(cores['fof_tag'], return_counts = True)
    srt = np.argsort(u)
    u=u[srt]
    cnt = cnt[srt]
    
    print("unique: ",u, cnt)


def plot_a_b(core_data, a, a_label, b, b_label, a_bins=None, b_bins=None,yscale=None,xscale=None, plot_median = False):
    plt.figure()
    if(a_bins is None):
        xbins = 100
    else:
        xbins = a_bins
    if(b_bins is None):
        ybins = 100
    else:
        ybins = b_bins

    h,xbins,ybins = np.histogram2d(core_data[a], core_data[b], bins=(xbins, ybins))
    plt.pcolor(xbins, ybins, h.T, cmap='PuBu', norm=clr.LogNorm())
    cb = plt.colorbar()
    cb.set_label('Count')
    plt.xlabel(a_label);plt.ylabel(b_label)
    if xscale == 'log':
        plt.xscale('log')
    if yscale == 'log':
        plt.yscale('log')
    if plot_median:
        cnt = np.sum(h,axis=1)
        slct = cnt>10
        # def get_bounds(data):
        #     print(data)
        #     if(len(data) > 0):
        #         return np.percentile(data,[0.159,0.5,0.841])
        #     else:
        #         return np.array([0,0,0])
        # data_bounds = dtk.binned_function(core_data[a], core_data[b], xbins, get_bounds)
        # print(data_bounds)
        # plt.plot(dtk.bins_avg(xbins)[slct], data_bounds[slct][1], 'r', label='median')
        # plt.plot(dtk.bins_avg(xbins)[slct], data_bounds[slct][0], '--r', label='1 sigma')
        # plt.plot(dtk.bins_avg(xbins)[slct], data_bounds[slct][2], '--r')
        median = dtk.binned_median(core_data[a], core_data[b], xbins)
        plt.plot(dtk.bins_avg(xbins)[slct], median[slct], 'r', label='median')
        plt.legend(loc='best', framealpha=0)
    plt.tight_layout()


def plot_a(core_data, a, a_label, a_bins=None, xscale=None, yscale=None):
    plt.figure()
    if(a_bins is None):
        a_bins  =100
    h,xbins = np.histogram(core_data[a], bins=a_bins)
    bin_widths = dtk.bins_width_dex(xbins)
    h = h/bin_widths
    plt.plot(dtk.bins_avg(xbins), h, '-',label=a_label)
    plt.grid()
    plt.xlabel(a_label)
    plt.ylabel('Count')
    if not(xscale is None):
        plt.xscale(xscale)
    if not(yscale is None):
        plt.yscale(yscale)
    plt.tight_layout()

if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    core_loc = param.get_string('core_loc')
    steps = param.get_int_list('steps')
    step = param.get_int('step')

    fof_loc = param.get_string('fof_loc')
    stepz = dtk.StepZ(sim_name = 'AlphaQ')
    cores = load_cores(core_loc.replace("${step}",str(step)),0, stepz)
    fof = load_fof(fof_loc.replace("${step}",str(step)))
    print('left join dic')
    cores  = left_join_dic(cores,fof,'fof_tag')
    cores_sat = select_dic(cores,cores['central']==0)
    cores_cen = select_dic(cores,cores['central']==1)

    # redshift = stepz.get_z(steps)[::-1]
    redshift = np.sort(np.unique(cores['infall_redshift']))
    print(redshift)
    scale_factor = stepz.get_a(steps)+0.001
    redshift = redshift[redshift<8]
    print('processing cores')
    process_cores(cores, cores_cen)
    print('plot1')
    plot_a_b(cores_sat,'infall_step', 'infall step', 
             'infall_mass', 'Infall Mass', 
             np.array(steps), np.logspace(10,15,100),
             yscale='log',
             plot_median=True)
    print('plot2')
    plot_a_b(cores_sat, 'fof_mass', 'FoF Halo Mass [Msun/h]',  
             'infall_mass', 'Infall Mass [Musn/h]',
             np.logspace(10,15,100), np.logspace(10,15,100), yscale='log', xscale='log',
             plot_median=True)
    print('plot3')
    plot_a_b(cores_sat, 'fof_mass', 'FoF Halo Mass [Msun/h]',  
             'radius', 'Core Radius [Mpc/h]',
             np.logspace(10,15,100), np.logspace(-3,0,100), yscale='log', xscale='log',
             plot_median=True)
    print('plot3')
    plot_a_b(cores_sat, 'infall_mass', 'Infall Mass [Msun/h]',  
             'radius', 'Core Radius [Mpc/h]',
             np.logspace(10,15,100), np.logspace(-3,0,100), yscale='log', xscale='log',
             plot_median=True)
    print('plot3.1')
    plot_a_b(cores, 'infall_mass', 'Infall Mass [Msun/h]',  
             'radius', 'Core Radius [Mpc/h]',
             np.logspace(10,15,100), np.logspace(-4,0,100), yscale='log', xscale='log')
    print('plot4')
    plot_a_b(cores_sat, 'infall_step', 'Infall Step',  
             'radius', 'Core Radius [Mpc/h]',
             np.array(steps), np.logspace(-3,0,100), yscale='log',
             plot_median=True)
    print('plot5')


    plot_a_b(cores_sat, 'infall_redshift+1', "Infall Redshift + 0.1",
             'radius', 'Core Radius [Mpc/h]', 
             redshift+.1,
             np.logspace(-3,0,100), yscale='log', xscale='log',)

    plot_a_b(cores_sat, 'infall_redshift', "Infall redshift",
             'radius', 'Core Radius [Mpc/h]', 
             redshift,
             np.logspace(-3,0,100), yscale='log', )


    print('plot6')
    plot_a(cores_sat, 'infall_redshift', "Infall Redshift",
           a_bins=redshift)
    print('plot7')
    plot_a(cores_sat, 'radius', "Core Radius",
           a_bins = np.logspace(-3,0,100),
           xscale = 'log')
    print('plot8')
    plot_a(cores_sat, 'infall_mass', "Infall Mass",
           a_bins = np.logspace(10,15,100),
           xscale = 'log')
    print('plot9')
    plot_a(cores_sat, 'infall_redshift', "Infall Redshift",
           a_bins=redshift,
           yscale='log')
    print('plot10')
    plot_a(cores_sat, 'radius', "Core Radius",
           a_bins = np.logspace(-3,0,100),
           xscale = 'log',
           yscale='log')
    print('plot11')
    plot_a(cores_sat, 'infall_mass', "Infall Mass",
           a_bins = np.logspace(10,15,100),
           xscale = 'log',
           yscale='log'   )

    ###### Saving figs #######
    dtk.save_figs("figs/"+__file__+"/"+sys.argv[1]+"/")
    plt.show()
    
