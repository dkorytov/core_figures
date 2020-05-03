#!/usr/bin/env python2.7

from __future__ import print_function, division

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk

from util import *


class AllCore:
    def __init__(self, core_loc, steps):
        self._load_all_cores(core_loc, steps)
        self._preprocess_cores()
        
    def _load_all_cores(self, core_loc, steps):
        all_cores_list = []
        for step in steps:
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


def plot_halo_assembly(param_file):
    param = dtk.Param(param_file)
    core_loc = param.get_string('core_loc')
    sod_loc = param.get_string('fof_loc')
    step = param.get_int('step')
    steps = param.get_int_list('steps')
    sod_loc = param.get_string('sod_loc')

    sod_data = load_sod(sod_loc.replace('${step}', str(step)))
    core_data = load_cores(core_loc.replace('${step}', str(step)))

    all_core = AllCore(core_loc,steps)

    for i in range(0,20):
        print(i, all_core.master_core_dic['step'][i], all_core.master_core_dic['core_tag'][i])

    for core_tag in [193654783976931328, 193654783976931330, 193654783976931331, 193654783976931332]:
        x = all_core.get_var(core_tag, 'x')
        y = all_core.get_var(core_tag, 'y')
        plt.figure()
        plt.plot(x,y)
    plt.show()

if __name__ == "__main__":
    plot_halo_assembly(sys.argv[1])
