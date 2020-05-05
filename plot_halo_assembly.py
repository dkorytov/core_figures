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


class AllCore:
    def __init__(self):
        pass
    
    def save_to_cache(self, cache='cache/all_cores.hdf5'):
        dtk.save_dict_hdf5(cache, self.master_core_dic)

    def load_from_cache(self, cache='cache/all_cores.hdf5'):
        self.master_core_dic = dtk.load_dict_hdf5(cache)
        
    def load_cores_from_file(self, core_loc, steps):
        self._load_all_cores(core_loc, steps)
        self._preprocess_cores()
        
    def _load_all_cores(self, core_loc, steps):
        all_cores_list = []
        for i, step in enumerate(steps):
            print("{}/{}".format(i, len(steps)))
            cores_data = load_cores(core_loc.replace('${step}', str(step)), sort=False, step=step)
            all_cores_list.append(cores_data)
        self.master_core_dic = concat_dics(all_cores_list)

    def _preprocess_cores(self,):
        # print(np.shape(self.master_core_dic['core_tag']))
        # print(np.shape(self.master_core_dic['step']))
        srt = np.lexsort((self.master_core_dic['step'], self.master_core_dic['core_tag'],))
        reorder_dic(self.master_core_dic, srt)
        
    def get_var(self, core_tag, variable):
        right_index = np.searchsorted(self.master_core_dic['core_tag'], core_tag, side='right')
        left_index = np.searchsorted(self.master_core_dic['core_tag'], core_tag, side='left') 
        return self.master_core_dic[variable][left_index:right_index]

def plot_individual_halo(sod_data, core_data, all_core, sod_index, zoom=False):
    fof_halo_tag = sod_data['fof_tag'][sod_index]
    x = sod_data['x'][sod_index]
    y = sod_data['y'][sod_index]
    r200 = sod_data['sod_halo_radius'][sod_index]
    m200 = sod_data['sod_halo_mass'][sod_index]
    print("mass: {:.2e}\nradius: {:.2f}".format(m200, r200))
    core_tags = core_data['core_tag'][core_data['fof_tag']==fof_halo_tag]
    print(core_tags)
    plt.figure()


    for core_tag in core_tags:
        core_x = all_core.get_var(core_tag, 'x')
        core_y = all_core.get_var(core_tag, 'y')
        xx = (core_x-x)
        xx[xx>128.0]=xx[xx>128.0]-256.0
        xx[xx<-128.0]=xx[xx<-128.0]+256.0
        yy = (core_y-y)
        yy[yy>128.0]=yy[yy>128.0]-256.0
        yy[yy<-128.0]=yy[yy<-128.0]+256.0

        plt.plot(xx, yy,'-', color='tab:orange', alpha=0.025)
        
    plt.axes().set_aspect('equal','datalim')
    
    plt.plot([], [], 'k-', label='R$_{200c}$')
    plt.plot([], [], '-', color='tab:orange', alpha=0.66, label='core trajectory' )
    plt.text(0.975, 0.025, 'M$_{{200c}}$ = {:2.2e}[$h_{{-1}}$M$_{{200c}}$]'.format(m200),
             transform=plt.gca().transAxes, verticalalignment='bottom',
             horizontalalignment='right')
    
    plt.legend(framealpha=0.0)
    c = plt.Circle((0,0), r200, color='k', fill=False)
    plt.gca().add_artist(c)
    plt.ylabel('y [$h^{-1}$ Mpc]')
    plt.xlabel('x [$h^{-1}$ Mpc]')
    plt.tight_layout()

    dtk.save_figs('figs/'+__file__+'/'+sys.argv[1]+'/', '.pdf')
    plt.close('all')

    
def plot_halo_assembly(param_file):
    param = dtk.Param(param_file)
    core_loc = param.get_string('core_loc')
    sod_loc = param.get_string('fof_loc')
    step = param.get_int('step')
    steps = param.get_int_list('steps')
    sod_loc = param.get_string('sod_loc')

    sod_data = load_sod(sod_loc.replace('${step}', str(step)))
    core_data = load_cores(core_loc.replace('${step}', str(step)))
    all_core = AllCore()
    # all_core.load_cores_from_file(core_loc,steps)
    # all_core.save_to_cache()
    all_core.load_from_cache()

    for i in range(0,len(sod_data['sod_halo_mass'])):
        if sod_data['sod_halo_mass'][i] > 1e14:
            plot_individual_halo(sod_data, core_data, all_core, i)
            plt.show()
    # for i in range(0,20):
    #     print(i, all_core.master_core_dic['step'][i], all_core.master_core_dic['core_tag'][i])

    # for core_tag in [193654783976931328, 193654783976931330, 193654783976931331, 193654783976931332]:
    #     x = all_core.get_var(core_tag, 'x')
    #     y = all_core.get_var(core_tag, 'y')
    #     plt.figure()
    #     plt.plot(x,y)
    plt.show()

if __name__ == "__main__":
    plot_halo_assembly(sys.argv[1])
