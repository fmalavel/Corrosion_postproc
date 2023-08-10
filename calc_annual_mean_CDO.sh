path_2022met='/data/users/fmalavel/AQUM/PROJECT/Corrosion_Scotland/verification/2022_mi-bf351'
file_2022met_2000ems='2022met_2000ems_gridded_cube_list'
file_2022met_2019ems='2022met_2019ems_gridded_cube_list'

cdo daymean ${path_2022met}'/'${file_2022met_2000ems}'.nc' ${file_2022met_2000ems}'_daily_CDO.nc'
cdo monmean ${file_2022met_2000ems}'_daily_CDO.nc' ${file_2022met_2000ems}'_monthly_CDO.nc' 

path_2017met='/data/users/fmalavel/AQUM/PROJECT/Corrosion_Scotland/verification/2017_mi-bf350'
file_2017met_2000ems='2017met_2000ems_gridded_cube_list'
file_2017met_2019ems='2017met_2019ems_gridded_cube_list'

#cdo daymean ${path_2022met}'/'${file_2022met_2019ems}'_CDO.nc' ${file_2022met_2019ems}'_daily_CDO.nc'
#cdo monmean ${file_2017met_2000ems}'_daily_CDO.nc' ${file_2017met_2000ems}'_monthly_CDO.nc' 
