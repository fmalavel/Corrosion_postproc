import os
import iris
import iris.coord_categorisation as cat
import iris.cube
iris.FUTURE.datum_support = True


sim_dict = {
    "met2017_ems2000" : {"met_year":2017,
                         "ems_year":2000,
                         "verif_ids":"mi-bf350"},
    "met2017_ems2019" : {"met_year":2017,
                         "ems_year":2019,
                         "verif_ids":"mi-bf350"},
    "met2022_ems2000" : {"met_year":2022,
                         "ems_year":2000,
                         "verif_ids":"mi-bf351"},
    "met2022_ems2019" : {"met_year":2022,
                         "ems_year":2019,
                         "verif_ids":"mi-bf351"}
}

aqum_data_path = '/data/users/fmalavel/AQUM/PROJECT/' + \
        'Corrosion_Scotland/verification/'

#tavg_output_path = aqum_data_path
tavg_output_path = '/data/users/fmalavel/AQUM/PROJECT/' + \
        'Corrosion_Scotland/Corrosion_postproc/AQUM_time_mean_outputs/'

if not os.path.exists(tavg_output_path):
  os.mkdir(tavg_output_path)
  print(f"\nFolder {tavg_output_path} created!")
else:
  print(f"\nFolder {tavg_output_path} already exists")


def calc_monthly_and_annual_mean(met_year=None, ems_year=None, verif_id=None):

    if met_year is None:
        raise ValueError("'met_year' for the run to be processed is missing")
    if ems_year is None:
        raise ValueError("'ems_year' for the run to be processed is missing")
    if vid is None:
        raise ValueError("'vid' for the run to be processed is missing")

    ems_year = ems_year
    met_year = met_year
    
    if met_year == 2017:
        input_data_path = aqum_data_path + str(met_year) + "_" + vid + "/"
    elif met_year == 2022:
        input_data_path = aqum_data_path + str(met_year) + "_" + vid + "/"
    else:
        raise ValueError(f"Output for {met_year} met_year do not exist")
    
    cubelist_by_month_out = iris.cube.CubeList()
    cubelist_by_year_out = iris.cube.CubeList()

    yr_constraint = iris.Constraint(
        time=lambda cell: cell.point.year == met_year)
    
    cl = iris.load(input_data_path + \
                   str(met_year) + "met_" + str(ems_year) + \
                   "ems_gridded_cube_list.nc",
                   yr_constraint)
    print(cl)
    
    for cube in cl:
        print("\n-------------------------------------------------\n")
        print(f'Processing {cube.attributes["short_name"]}')
        #print(cube)
        
        if not cube.coord('grid_latitude').has_bounds():
            print(f'adding bounds to grid_latitude coord')
            cube.coord('grid_latitude').guess_bounds()
        if not cube.coord('grid_longitude').has_bounds():
            print(f'adding bounds to grid_longitude coord')
            cube.coord('grid_longitude').guess_bounds()
        if not cube.coord('time').has_bounds():
            print(f'adding bounds to time coord')
            cube.coord('time').guess_bounds()
    
        cube.remove_coord('forecast_period')
        cube.remove_coord('forecast_day')
        #cube = cube.collapsed('time', iris.analysis.MEAN)
        #print(cube)
        #print(cube.data)
        #cubelist_out.append(cube)
    
        iris.coord_categorisation.add_year(cube, 'time', name='year')
        iris.coord_categorisation.add_month(cube, 'time', name='month')
        #print(cube)
    
        cm = cube.aggregated_by(['month'], iris.analysis.MEAN)
        #print(cm)
        cy = cube.aggregated_by(['year'], iris.analysis.MEAN)
        #print(cy)
    
        cubelist_by_month_out.append(cm)
        cubelist_by_year_out.append(cy)
    
    print(cubelist_by_month_out)
    print(cubelist_by_year_out)
    
    file_mon_out = str(met_year) + "met_" + str(ems_year) + \
            "ems_gridded_cube_list_monthly_mean.nc"
    
    file_year_out = str(met_year) + "met_" + str(ems_year) + \
            "ems_gridded_cube_list_annual_mean.nc"
    
    iris.save(cubelist_by_month_out, tavg_output_path + file_mon_out)
    iris.save(cubelist_by_year_out, tavg_output_path + file_year_out)
    

if __name__ == "__main__":

    for run_dic in sim_dict:
        met = sim_dict[run_dic]['met_year']
        ems = sim_dict[run_dic]['ems_year']
        vid = sim_dict[run_dic]['verif_ids']

        calc_monthly_and_annual_mean(met_year=met, ems_year=ems, verif_id=vid)
