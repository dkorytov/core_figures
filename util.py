import numpy as np
import dtk

from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)

def sort_dic(dic,sort_key):
    srt = np.argsort(dic[sort_key])
    for key in dic.keys():
        dic[key] = dic[key][srt]

def select_dic(dic,slct):
    result = {}
    for key in dic.keys():
        result[key] = dic[key][slct]
    return result
        
def combine_dics(dic1, dic2, key):
    sort_dic(dic1, key)
    srt = dtk.search_sorted(dic1[key], dic2[key])
    found = srt != -1
    dic_result = {}
    for key in dic1.keys():
        dic_result[key] = dic1[key][srt][found]
    for key in dic2.keys():
        dic_result[key] = np.copy(dic2[key][found])
    return dic_result

def plot_subhalo_dm_cores():
    pass
    return

def load_cores(core_loc):
    print("loading cores...")
    result = {}
    result['fof_tag'] = dtk.gio_read(core_loc,'fof_halo_tag')
    result['x']       = dtk.gio_read(core_loc,'x')
    result['y']       = dtk.gio_read(core_loc,'y')
    result['z']       = dtk.gio_read(core_loc,'z')
    result['infall_mass'] = dtk.gio_read(core_loc,'infall_tree_node_mass')
    result['radius']  = dtk.gio_read(core_loc,'radius')
    result['infall_step'] = dtk.gio_read(core_loc,'infall_step')
    result['central']   = dtk.gio_read(core_loc,'central')
    sort_dic(result,'fof_tag')
    return result

def load_fof(fof_loc, min_mass = None):
    print("loading fof...")
    result = {}
    result['fof_tag'] = dtk.gio_read(fof_loc,'fof_halo_tag')
    result['x'] = dtk.gio_read(fof_loc,'fof_halo_center_x')
    result['y'] = dtk.gio_read(fof_loc,'fof_halo_center_y')
    result['z'] = dtk.gio_read(fof_loc,'fof_halo_center_z')
    result['fof_mass'] = dtk.gio_read(fof_loc,'fof_halo_mass')
    if min_mass is not None:
        slct = result['fof_mass']>min_mass
        result = select_dic(result,slct)
        print("{:.2e}, {:.2e}".format(np.min(result['fof_mass']), min_mass))
    sort_dic(result,'fof_tag')
    return result

def load_sod(sod_loc):
    result = {}
    result['fof_tag'] = dtk.gio_read(sod_loc, 'fof_halo_tag')
    result['sod_halo_radius'] = dtk.gio_read(sod_loc, 'sod_halo_radius')
    result['sod_halo_mass'] = dtk.gio_read(sod_loc, 'sod_halo_mass')
    return result

def load_subhalo(sh_loc):
    print("loading subhalos")
    result = {}
    # dtk.gio_inspect(sh_loc)
    result['fof_tag'] = dtk.gio_read(sh_loc,'fof_halo_tag')
    result['x']       = dtk.gio_read(sh_loc,'subhalo_mean_x')
    result['y']       = dtk.gio_read(sh_loc,'subhalo_mean_y')
    result['z']       = dtk.gio_read(sh_loc,'subhalo_mean_z')
    result['subhalo_mass']    = dtk.gio_read(sh_loc,'subhalo_mass')
    result['subhalo_tag']     = dtk.gio_read(sh_loc, 'subhalo_tag')
    slct = result['subhalo_tag'] != 0
    result = select_dic(result, slct)
    sort_dic(result,'fof_tag')
    return result
    
def load_bighalo(bh_loc):
    result = {}
    print('loading big halo fof tag')
    result['fof_tag'] = dtk.gio_read(bh_loc,'fof_halo_tag')
    print('loading big halo x')
    result['x']       = dtk.gio_read(bh_loc,'x')
    print('loading big halo y')
    result['y']       = dtk.gio_read(bh_loc,'y')
    print('loading big halo z')
    result['z']       = dtk.gio_read(bh_loc,'z')
    print(result['fof_tag'])
    print(result['x'])
    print(result['y'])
    print(result['z'])
    # sort_dic(result,'fof_tag')
    # print(result.keys())
    # print(result['fof_tag'].size)
    # print(result['x'].size)
    return result

def fix_periodic_boundary(del_x, period):
    slct = del_x > period/2.0
    del_x[slct] = del_x[slct] -period
    slct = del_x < -period/2.0
    del_x[slct] = del_x[slct] + period

def get_diff(data, fof_data, test=False):
    print('getting difference')
    sort_dic(fof_data, 'fof_tag')
    index = dtk.search_sorted(fof_data['fof_tag'], data['fof_tag'])
    print('\tdone sorting')
    slct = index != -1
    if(test):
        print("fof size: ", np.shape(fof_data['fof_tag']))
        print("data size: ", np.shape(data['fof_tag']))
        print("slct size: ", np.shape(slct), "sum: ", np.sum(slct))

    data_new = select_dic(data,slct)
    index_new = index[slct]
    del_x = data_new['x']-fof_data['x'][index_new]
    del_y = data_new['y']-fof_data['y'][index_new]
    del_z = data_new['z']-fof_data['z'][index_new]
    fix_periodic_boundary(del_x, 256.0)
    fix_periodic_boundary(del_y, 256.0)
    fix_periodic_boundary(del_z, 256.0)
    del_r = np.sqrt(del_x**2 + del_y**2 + del_z**2)
    data_new['del_r'] = del_r
    data_new['del_r_r200'] = del_r/fof_data['sod_halo_radius'][index_new]
    print('\tdone')
    return data_new

def get_unique_count(data, key):
    unique = np.unique(data[key])
    return len(unique)
