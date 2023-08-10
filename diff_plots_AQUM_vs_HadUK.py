import iris
import iris.cube
import iris.plot as iplt
import iris.quickplot as qplt
from iris.time import PartialDateTime
from iris.analysis.cartography import rotate_pole, unrotate_pole
iris.FUTURE.datum_support = True

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import BoundaryNorm

import cartopy.crs as ccrs
import numpy as np
import cf_units
import re
import sys

import os
import time

from regrid_AQUM_to_HadUK import extract_cube_area, precision_round


output_path = 'Regridded_outputs/Diffs_plots/'

# check whether output directory already exists
if not os.path.exists(output_path):
    os.mkdir(output_path)
    print(f"\nFolder {output_path} created!")
else:
    print(f"\nFolder {output_path} already exists")


def plot_fields_diffs(cube_ctl,
                      cube_exp,
                      colormap=None,
                      cb_min = None,
                      cb_max = None,
                      cb_diff_min = None,
                      cb_diff_max = None,
                      short_str=None,
                      ctl_str=None,
                      exp_str=None,
                      filename_prefix=None,
                      output_path=None):
    """
    Contour plot the three cubes provided to the function.
    :cube_ctl:      iris cube to be plotted
    :cube_exp:      iris cube to be plotted
    :colormap:      matplotlib colormap
    :cb_min:        minimum value displayed on contour plots
    :cb_max:        maximum value displayed on contour plots
    :cb_diff_min:   minimum value displayed on diff contour plots
    :cb_diff_max:   maximum value displayed on diff contour plots
    :short_str:     used in contour plots title and plot output filename
    """

    if cube_ctl is None:
        raise ValueError(f"\n No cube provided for cube_ctl")
    if cube_exp is None:
        raise ValueError(f"\n No cube provided for cube_exp")
    
    if not isinstance(cube_ctl, iris.cube.Cube):
        raise ValueError(f"{cube_ctl} is not an Iris cube")
    
    if not isinstance(cube_exp, iris.cube.Cube):
        raise ValueError(f"{cube_exp} is not an Iris cube")
    
    cube_p1 = cube_ctl
    cube_p2 = cube_exp
    cube_p3 = 100 * (cube_p1 - cube_p2) / cube_p1

    print(f"\nCubes to be plotted:")
    print(cube_p1.summary(True))
    print(cube_p2.summary(True))

    if colormap is None:
        cmap = plt.get_cmap("Spectral")
        cmap_r = cmap.reversed()
        cmap = cmap_r
    else:
        cmap = plt.get_cmap(colormap)

    cmap_diff = plt.get_cmap("bwr")

    if cb_min is None:
        cb_min = precision_round(np.percentile(cube_p1.data, 0.1),0)
    if cb_max is None:
        cb_max = precision_round(np.percentile(cube_p1.data, 99.9),0)
    
    if (cb_diff_min is None) or (cb_diff_max is None):
        diff_min = precision_round(np.percentile(cube_p3.data, 0.1),0)
        diff_max = precision_round(np.percentile(cube_p3.data, 99.9),0)
        
        cb_diff_min = -1 * max(abs(diff_min),abs(diff_max))
        cb_diff_max =  1 * max(abs(diff_min),abs(diff_max))

    levels = np.linspace(cb_min, cb_max, 21)
    tick_lvl = levels[0::5]
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)

    levels_diff = np.linspace(cb_diff_min, cb_diff_max, 21)
    tick_lvl_diff = levels_diff[0::5]
    norm_diff = BoundaryNorm(levels_diff, ncolors=cmap.N, clip=False)

#   Start of figure
    fig = plt.figure(figsize=(12,8), dpi=100)

#   Fig 1
    plt.subplot(1, 3, 1)
    iplt.contourf(cube_p1, levels=levels, cmap=cmap, 
                  norm=norm, extend='both')
    plt.gca().coastlines(resolution='10m')

    cb_label_str = short_str + ' [' + str(cube_p1.units) + ']'
    plt.colorbar(orientation='horizontal', 
                 ticks=tick_lvl,
                 # format="%.1e",
                 pad=0.05, 
                 aspect=30, 
                 shrink=0.9).set_label(cb_label_str)
    suptitle = ctl_str
    plt.title(suptitle, fontsize=10, linespacing=1.2)

#   Fig 2
    plt.subplot(1, 3, 2)
    iplt.contourf(cube_p2, levels=levels, cmap=cmap, 
                  norm=norm, extend='both')
    plt.gca().coastlines(resolution='10m')

    cb_label_str = short_str + ' [' + str(cube_p2.units) + ']'
    plt.colorbar(orientation='horizontal', 
                 ticks=tick_lvl,
                 # format="%.1e",
                 pad=0.05, 
                 aspect=30, 
                 shrink=0.9).set_label(cb_label_str)
    suptitle = exp_str
    plt.title(suptitle, fontsize=10, linespacing=1.2)

#   Fig 3
    plt.subplot(1, 3, 3)
    iplt.contourf(cube_p3, levels=levels_diff, cmap=cmap_diff, 
                  norm=norm_diff, extend='both')
    plt.gca().coastlines(resolution='10m')

    #cb_label_str = 'Differences' + ' [' + str(cube_p3.units) + ']'
    cb_label_str = 'Differences' + ' [%]'
    plt.colorbar(orientation='horizontal', 
                 ticks=tick_lvl_diff,
                 # format="%.1e",
                 pad=0.05, 
                 aspect=30, 
                 shrink=0.9).set_label(cb_label_str)
    suptitle = [ctl_str + ' - ' + exp_str]
    plt.title(suptitle, fontsize=10, linespacing=1.2)
   
#   Export plot to png
    if filename_prefix:
        plot_filename = output_path + \
                filename_prefix + \
                ctl_str + '_vs_' + exp_str + '_' + short_str + '.png'
    else:
        plot_filename = output_path + \
                ctl_str + '_vs_' + exp_str + '_' + short_str + '.png'
        
    print(f'Saving plots as: ' + plot_filename)
    plt.savefig(plot_filename)


#   ##############################
#               START
#   ##############################

if __name__ == "__main__":

    print("\n-------------------------------------" + \
          "\nStart of diff_plots_AQUM_vs_HadUK.py:" + \
          "\n-------------------------------------"
         )

    stall = time.time()

#   #############################
#      Load AQUM fields
#   #############################

    print("\n------------------" + \
          "\nLoading AQUM fiedls:" + \
          "\n------------------"
         )
    st = time.time()
    
    # AQUM data path and filename from parsed arguments.
    # TODO: pass path/filename as environement input arguments
    ctl_data_path = 'Regridded_outputs/'
    ctl_data_filename = 'subdomain_2022met_AQUM_regridded_time_avg_after.nc'

    exp_data_path = 'Regridded_outputs/'
    exp_data_filename = 'subdomain_2022met_AQUM_regridded_time_avg_before.nc'

    ctl_input_file = ctl_data_path + '/' + ctl_data_filename
    exp_input_file = exp_data_path + '/' + exp_data_filename

    # Derive Met Year to be processed based on filename
    pattern = '^[0-9][0-9][0-9][0-9]met'
    match = re.match(pattern, ctl_data_filename)
    met_year = 2022

    pdt1 = PartialDateTime(year=met_year, month=1, day=1)
    pdt2 = PartialDateTime(year=met_year, month=12, day=31)
    year_constraint = iris.Constraint(time=
                                      lambda cell: pdt1 <= cell.point < pdt2)

    # List of cube short names to be loaded (consider passing as argument)
    # short_name_list = None
    short_name_list = ["T_surf", "RHw", "total_precip",
                       "u_10m", "v_10m",
                       "SO2", "SO2_DryDep", "SS_WetDep",
                       "O3", "O3_DryDep", 
                       "SS_Conc", "SS_DryDep", "SS_WetDep"]

    short_name_list = ["T_surf", "RHw", 
                       "O3", "SO2", 
                       "wind_speed"]

    print(f'ctl_data_path: {ctl_data_path}')
    print(f'ctl_data_filename: {ctl_data_filename}')
    print(f'exp_data_path: {exp_data_path}')
    print(f'exp_data_filename: {exp_data_filename}')
    print(f'met_year: {met_year}')

    if short_name_list:
        #load cubes matching predefined short names
        short_constraint = iris.AttributeConstraint(
            short_name=lambda c: c in short_name_list
        )
        ctl_input_cube_list = iris.load(ctl_input_file, short_constraint)
        exp_input_cube_list = iris.load(exp_input_file, short_constraint)
    else:
        #load everything in the file
        ctl_input_cube_list = iris.load(ctl_input_file)
        exp_input_cube_list = iris.load(exp_input_file)
    
    print(f"\nctl_input_cube_list:\n{ctl_input_cube_list}")
    print(f"\nexp_input_cube_list:\n{exp_input_cube_list}")
    et = time.time()
    elapsed_time = et - st
    print('\nExecution time:', elapsed_time, 'seconds')


#   #############################
    # Load an annual mean gridded 
    # field of tas from HadUK-GRID 
    # @1km resolution
#   #############################

    print("\n------------------\nLoading HadUK-GRID:\n------------------\n")
    st = time.time()

    HadUK_data_path = '/data/users/fmalavel' + \
            '/AQUM/PROJECT/Corrosion_Scotland/HadUK-GRID/'

    # Surface temp
    var = 'tas'
    haduk_1km_annual_mean = HadUK_data_path + var + '_hadukgrid_uk_1km_ann_' + \
            str(met_year) + '01' + '-' + str(met_year) + '12.nc'
    haduk_cube_tas = iris.load_cube(haduk_1km_annual_mean, year_constraint)

    # Surface temp
    var = 'hurs'
    haduk_1km_annual_mean = HadUK_data_path + var + '_hadukgrid_uk_1km_ann_' + \
            str(met_year) + '01' + '-' + str(met_year) + '12.nc'
    haduk_cube_hurs = iris.load_cube(haduk_1km_annual_mean, year_constraint)

    # 10m surface wind
    var = 'sfcWind'
    haduk_1km_annual_mean = HadUK_data_path + var + '_hadukgrid_uk_1km_ann_' + \
            str(met_year) + '01' + '-' + str(met_year) + '12.nc'
    haduk_cube_sfcWind = iris.load_cube(haduk_1km_annual_mean, year_constraint)

    print(haduk_cube_tas)
    et = time.time()
    elapsed_time = et - st
    print('\nExecution time:', elapsed_time, 'seconds')


#   #############################
#    Decrease target domain size
#   #############################

   #if ctl_data_filename.split('_')[0] == 'subdomain':
    if re.search('subdomain', ctl_data_filename):
        print("\n-----------------------" + \
              "\nSubsampling HadUK data:" + \
              "\n-----------------------")

        st = time.time()

        haduk_cube_tas = extract_cube_area(haduk_cube_tas)
        haduk_cube_hurs = extract_cube_area(haduk_cube_hurs)
        haduk_cube_sfcWind = extract_cube_area(haduk_cube_sfcWind)

        print(f"\nhaduk_cube: {haduk_cube_tas}")

        et = time.time()
        elapsed_time = et - st
        print('\nExecution time:', elapsed_time, 'seconds')


#   #############################
#               PLOTS
#       AQUM: before - after
#   #############################
    
    print("\n------------------\nCreating Plots:\n------------------")
    st = time.time()

    ctl_str = 'AQUM_t_avg_after_regrid'
    exp_str = 'AQUM_t_avg_before_regrid'

    if re.search('subdomain', exp_data_filename):
        prefix = "subdomain_" + str(met_year) + "met_"
    else:
        prefix = str(met_year) + "met_"

    # Loop over a cubelist of aqum fields
    for short in short_name_list:

        short_constrain = iris.AttributeConstraint(short_name=short)

        cube_ctl = ctl_input_cube_list.extract_cube(short_constrain)
        cube_exp = exp_input_cube_list.extract_cube(short_constrain)
        
        cb_min=None
        cb_max=None
        colormap=None

        if short == "T_surf":
            colormap="Spectral"
        if short == "RHw":
            colormap="GnBu"
        if short in ["wind_speed"]:
            colormap="plasma"

        plot_fields_diffs(cube_ctl=cube_ctl,
                          cube_exp=cube_exp,
                          colormap=colormap,
                          cb_min=cb_min,
                          cb_max=cb_max,
                          short_str=short,
                          ctl_str=ctl_str,
                          exp_str=exp_str,
                          filename_prefix=prefix,
                          output_path=output_path
                         )

    et = time.time()
    elapsed_time = et - st
    print('\nExecution time:', elapsed_time, 'seconds')
    

#   #############################
#               PLOTS
#            AQUM - HadUK
#   #############################

    print("\n------------------\nCreating Plots:\n------------------")
    st = time.time()

    short_name_list = ["T_surf", "RHw",  "wind_speed"]

    ctl_str = 'HadUK-GRID'
    exp_str = 'AQUM_t_avg_before_regrid'

    if re.search('subdomain', exp_data_filename):
        prefix = "subdomain_" + str(met_year) + "met_"
    else:
        prefix = str(met_year) + "met_"

    # Loop over a cubelist of aqum fields
    for short in short_name_list:

        short_constrain = iris.AttributeConstraint(short_name=short)

        cube_exp = exp_input_cube_list.extract_cube(short_constrain)
        
        cb_min=None
        cb_max=None
        cb_diff_min=None
        cb_diff_max=None
        colormap=None

        if short == "T_surf":
            colormap="Spectral"
            cube_ctl = haduk_cube_tas
            cb_min=6
            cb_max=10
            cb_diff_min=-100
            cb_diff_max=100
        if short == "RHw":
            colormap="GnBu"
            cube_ctl = haduk_cube_hurs
            cb_min=80
            cb_max=90
            cb_diff_min=-100
            cb_diff_max=100
        if short in ["wind_speed"]:
            colormap="plasma"
            cube_ctl = haduk_cube_sfcWind
            cb_min=3
            cb_max=10
            cb_diff_min=-100
            cb_diff_max=100

        plot_fields_diffs(cube_ctl=cube_ctl,
                          cube_exp=cube_exp,
                          colormap=colormap,
                          cb_min=cb_min,
                          cb_max=cb_max,
                          cb_diff_min=cb_diff_min,
                          cb_diff_max=cb_diff_max,
                          short_str=short,
                          ctl_str=ctl_str,
                          exp_str=exp_str,
                          filename_prefix=prefix,
                          output_path=output_path
                         )

    et = time.time()
    elapsed_time = et - st
    print('\nExecution time:', elapsed_time, 'seconds')


#   ##############################
#                 END
#   ##############################

    print("\n------------------\nScript Summary Stats:\n------------------")

    etall = time.time()
    elapsed_time = etall - stall
    print('\nOverall script execution time:', elapsed_time, 'seconds')
