# Corrosion_postproc
Postprocessing code for AQUM outputs.

## Purpose
This repository contains a collection of python routines that are used to post-process AQUM outputs that have been created using ADAQ verification tools.

The code produces netcdf files containing the AQUM reprocessed fields to be used for deriving __corrosion rate maps__.  

## Current functionalities
- Calculate monthly and daily means from AQUM hourly model outputs.
- Regrid AQUM model outputs onto HadUK-GRID 1x1km grid
- Plots the regridded AQUM and HadUK-GRID fields

## Expected inputs data

The regridding code can process hourly, monthly and yearly AQUM data. 
Below is a example print of a loaded cubelist of AQUM fields to be regridded: 
```
0: sea_salt_wet_deposition / (kg/m2/s) (time: 8761; grid_latitude: 182; grid_longitude: 146)
1: sea_salt_dry_deposition / (kg/m2/s) (time: 8761; grid_latitude: 182; grid_longitude: 146)
2: ozone_dry_deposition / (kg/m2/s)    (time: 8761; grid_latitude: 182; grid_longitude: 146)
3: sulfur_dioxide_dry_deposition / (kg/m2/s) (time: 8761; grid_latitude: 182; grid_longitude: 146)
4: relative_humidity_wrt_water / (%)   (time: 8761; grid_latitude: 182; grid_longitude: 146)
5: mass_concentration_of_ozone_in_air / (ug/m3) (time: 8761; grid_latitude: 182; grid_longitude: 146)
6: mass_concentration_of_seasalt_dry_aerosol_in_air / (ug/m3) (time: 8761; grid_latitude: 182; grid_longitude: 146)
7: mass_concentration_of_sulfur_dioxide_in_air / (ug/m3) (time: 8761; grid_latitude: 182; grid_longitude: 146)
8: precipitation_amount / (kg m-2)     (time: 8761; grid_latitude: 182; grid_longitude: 146)
9: surface_temperature / (K)           (time: 8761; grid_latitude: 182; grid_longitude: 146)
10: x_wind / (m s-1)                    (time: 8761; grid_latitude: 183; grid_longitude: 146)
11: y_wind / (m s-1)                    (time: 8761; grid_latitude: 183; grid_longitude: 146)
12: wind_to_direction / (degrees)       (time: 8761; grid_latitude: 183; grid_longitude: 146)
13: wind_speed / (m s-1)                (time: 8761; grid_latitude: 183; grid_longitude: 146)
```

This will produce time averaged cubes interpolated onto HadUK-GRID projection coordinates, e.g. for ozone:
```
mass_concentration_of_ozone_in_air / (ug/m3) (projection_y_coordinate: 1450; projection_x_coordinate: 900)
    Dimension coordinates:
        projection_y_coordinate                                      x                              -
        projection_x_coordinate                                      -                              x
    Scalar coordinates:
        forecast_day                         1.0 Days
        forecast_period                      18.0 hours, bound=(6.0, 30.0) hours
        level_height                         20.000004 m, bound=(0.0, 36.664005) m
        model_level_number                   1
        sigma                                0.9977232, bound=(1.0, 0.99582815)
        time                                 2017-07-02 11:30:00, bound=(2016-12-31 23:00:00, 2018-01-01 00:00:00)
    Cell methods:
        mean                                 time (1 hour)
        mean                                 time
    Attributes:
        Conventions                          'CF-1.7'
        STASH                                m01s34i001
        label                                '2017met_2019ems'
        short_name                           'O3'
        source                               'Data from Met Office Unified Model'
```

## Dependencies
The code relies on [iris scitool](https://scitools-iris.readthedocs.io/en/stable). 

## Ressources and information
HadUK-GRID [data description](https://www.metoffice.gov.uk/research/climate/maps-and-data/data/haduk-grid/haduk-grid).  
HadUK-GRID [data access on CEDA](https://catalogue.ceda.ac.uk/uuid/4dc8450d889a491ebb20e724debe2dfb).  
AQUM [standard jobs](https://code.metoffice.gov.uk/trac/nwpscience/wiki/ResearchSuites/AQUM).  
