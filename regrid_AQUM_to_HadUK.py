# TODO: cache regridder --> https://scitools-iris.readthedocs.io/en/latest/userguide/interpolation_and_regridding.html#caching-a-regridder
# TODO: Lazy-data, chunk the time dimension --> https://scitools-iris.readthedocs.io/en/latest/userguide/interpolation_and_regridding.html#regridding-lazy-data
# TODO: ESMF regridding https://iris-esmf-regrid.readthedocs.io/en/latest/_api_generated/esmf_regrid.schemes.html#esmf_regrid.schemes.ESMFAreaWeighted

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

import tracemalloc
import linecache
import os
import time


output_path = 'Regridded_outputs/'

# check whether output directory already exists
if not os.path.exists(output_path):
  os.mkdir(output_path)
  print(f"\nFolder {output_path} created!")
else:
  print(f"\nFolder {output_path} already exists")


def display_top(snapshot, key_type='lineno', limit=10):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("\nTop %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        print("#%s: %s:%s: %.1f KiB"
              % (index, frame.filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))


def precision_round(numbers, digits = 3):
    '''
    Parameters:
    -----------
    numbers : scalar, 1D , or 2D array(-like)
    digits: number of digits after decimal point
    
    Returns:
    --------
    out : same shape as numbers
    '''
    import numpy as np

    numbers = np.asarray(np.atleast_2d(numbers))
    out_array = np.zeros(numbers.shape) # the returning array
    
    for dim0 in range(numbers.shape[0]):
        powers = [int(F"{number:e}".split('e')[1]) for number in numbers[dim0, :]]
        out_array[dim0, :] = [round(number, -(int(power) - digits))
                         for number, power in zip(numbers[dim0, :], powers)]
        
    # returning the original shape of the `numbers` 
    if out_array.shape[0] == 1 and out_array.shape[1] == 1:
        out_array = out_array[0, 0]
    elif out_array.shape[0] == 1:
        out_array = out_array[0, :]
    
    return out_array


def guess_coord_names(cube, axes):
    """
    Guess the name of the coordinate corresponding to the required axes

    :param cube: iris Cube
    :param axes: List of axes, eg 'X','Y','Z','T'

    :returns: List of coordinate names corresponding to these axes.
              If an axes not found, then value in list is None.
              Will try to return dimension coordinates if possible.
    """

    coord_names = [None] * len(axes)
    for coord in cube.coords():
        axis = iris.util.guess_coord_axis(coord)
        for i, ax in enumerate(axes):
            if axis == ax:
                if coord_names[i] is None:
                    coord_names[i] = coord.name()

    return coord_names


def plot_regridded_fields(native_aqum,
                          regrid_aqum,
                          haduk,
                          colormap=None,
                          cb_min = None,
                          cb_max = None,
                          title_str=None,
                          plot_filename_prefix=None):
    """
    Contour plot the three cubes provided to the function.
    :native_aqum:   iris cube of AQUM outputs on AQUM native grid
    :regrid_aqum:   iris cube of AQUM outputs on HadUK-GRID grid
    :haduk:         iris cube of HadUK-GRID data on 1km grid
    :colormap:      matplotlib colormap
    :cb_min:        minimum value displayed on contour plots
    :cb_max:        maximum value displayed on contour plots
    :title_str:     contour plots and plot filename
    """

    if native_aqum is None:
        raise ValueError(f"\n No cube provided for native_aqum")
    if regrid_aqum is None:
        raise ValueError(f"\n No cube provided for regrid_aqum")
    if haduk is None:
        raise ValueError(f"\n No cube provided for haduk")
    
    if not isinstance(native_aqum, iris.cube.Cube):
        raise ValueError(f"{native_aqum} is not an Iris cube")
    
    if not isinstance(regrid_aqum, iris.cube.Cube):
        raise ValueError(f"{regrid_aqum} is not an Iris cube")
    
    if not isinstance(haduk, iris.cube.Cube):
        raise ValueError(f"{haduk} is not an Iris cube")

    """
    # Calculate time mean:
    print(f"\nTime-averaging fields for plotting")
    cube_p1 = native_aqum.collapsed('time', iris.analysis.MEAN)
    cube_p2 = regrid_aqum.collapsed('time', iris.analysis.MEAN)
    """
    cube_p1 = native_aqum.copy()
    cube_p2 = regrid_aqum.copy()
    cube_p3 = haduk.copy()

    print(f"\nCubes to be plotted:")
    print(cube_p1.summary(True))
    print(cube_p2.summary(True))
    print(cube_p3.summary(True))

    if colormap is None:
        cmap = plt.get_cmap("Spectral")
        cmap_r = cmap.reversed()
        cmap = cmap_r
    else:
        cmap = plt.get_cmap(colormap)

    if cb_min is None:
        cb_min = precision_round(np.percentile(cube_p1.data, 0.1),0)
    if cb_max is None:
        cb_max = precision_round(np.percentile(cube_p1.data, 99.9),0)
    
    print(f"\ncb_min: {cb_min}")
    print(f"cb_max: {cb_max}\n")

    levels = np.linspace(cb_min, cb_max, 21)
    tick_lvl = levels[0::5]
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)

#   Start of figure
    fig = plt.figure(figsize=(12,8), dpi=100)

#   Fig 1
    plt.subplot(1, 3, 1)
    iplt.contourf(cube_p1, levels=levels, cmap=cmap, norm=norm, extend='both')
    plt.gca().coastlines(resolution='10m')

    cb_label_str = title_str + ' [' + str(cube_p1.units) + ']'
    plt.colorbar(orientation='horizontal', 
                 ticks=tick_lvl,
                 # format="%.1e",
                 pad=0.05, 
                 aspect=30, 
                 shrink=0.9).set_label(cb_label_str)
    suptitle = 'AQUM on native grid'
    plt.title(suptitle, fontsize=10, linespacing=1.2)

#   Fig 2
    plt.subplot(1, 3, 2)
    iplt.contourf(cube_p2, levels=levels, cmap=cmap, norm=norm, extend='both')
    plt.gca().coastlines(resolution='10m')

    cb_label_str = title_str + ' [' + str(cube_p2.units) + ']'
    plt.colorbar(orientation='horizontal', 
                 ticks=tick_lvl,
                 # format="%.1e",
                 pad=0.05, 
                 aspect=30, 
                 shrink=0.9).set_label(cb_label_str)
    suptitle = 'AQUM on HadUK-GRID'
    plt.title(suptitle, fontsize=10, linespacing=1.2)

#   Fig 3
    plt.subplot(1, 3, 3)
    iplt.contourf(cube_p3, levels=levels, cmap=cmap, norm=norm, extend='both')
    plt.gca().coastlines(resolution='10m')

    cb_label_str = cube_p3.name() + ' [' + str(cube_p3.units) + ']'
    plt.colorbar(orientation='horizontal', 
                 ticks=tick_lvl,
                 # format="%.1e",
                 pad=0.05, 
                 aspect=30, 
                 shrink=0.9).set_label(cb_label_str)
    suptitle = 'HadUK-GRID'
    plt.title(suptitle, fontsize=10, linespacing=1.2)
   
#   Export plot to png
    if plot_filename_prefix:
        plot_filename = output_path + plot_filename_prefix + title_str + '_time_avg.png'
    else:
        plot_filename = output_path + title_str + '_time_avg.png'
        
    print(f'Saving plots as: ' + plot_filename)
    plt.savefig(plot_filename)


def wind_speed(cubelist, short_name="ws"):
    """
    Calculate wind speed from x and y wind components (u and v)
    which are cubes within cubelist.
    Appends new cube to existing cubelist.

    :param cubelist: iris cubelist, must contain cubes named 'x_wind' and
                     'y_wind'
    :param short_name: required short_name to be added to output cube
    :returns: cubelist with extra wind speed cube included
    """

    try:
        x_wind = cubelist.extract_cube("x_wind")
    except iris.exceptions.ConstraintMismatchError:
        warnings.warn("wind speed not calculated: x_wind missing")
        return cubelist
    try:
        y_wind = cubelist.extract_cube("y_wind")
    except iris.exceptions.ConstraintMismatchError:
        warnings.warn("wind speed not calculated: y_wind missing")
        return cubelist

    # Check they are on same grid (U&V often output on different grids)
    assert x_wind.coord("grid_longitude") == y_wind.coord("grid_longitude")
    assert x_wind.coord("grid_latitude") == y_wind.coord("grid_latitude")

    ws = (x_wind**2 + y_wind**2) ** 0.5
    ws.convert_units("m s-1")
    ws.rename("wind_speed")
    if "label" in x_wind.attributes:
        ws.attributes["label"] = x_wind.attributes["label"]
    ws.attributes["short_name"] = short_name
    cubelist.append(ws)

    return cubelist


def wind_direction(cubelist, short_name="wind_dir"):
    """
    Calculate the direction wind is going to,
    using the x and y wind components (u and v)
    which are within the cubelist.
    Also rotates the winds to a regular lat-long grid if originally on a rotated grid.
    Appends new cube to existing cubelist

    :param cubelist: iris cubelist, must contain cubes named
                     'x_wind' and'y_wind'
    :param short_name: required short_name to be added to output cube
    :returns: cubelist with extra wind direction cube included
    """

    # Getting the x and y wind components
    try:
        x_wind = cubelist.extract_cube("x_wind")
    except iris.exceptions.ConstraintMismatchError:
        raise ValueError("wind direction not calculated, x_wind missing")
    try:
        y_wind = cubelist.extract_cube("y_wind")
    except iris.exceptions.ConstraintMismatchError:
        raise ValueError("wind direction not calculated, y_wind missing")
    
    # Checking that the wind components are on the same grid
    assert x_wind.coord("grid_longitude") == y_wind.coord("grid_longitude")
    assert x_wind.coord("grid_latitude") == y_wind.coord("grid_latitude")
    wind_dir = x_wind.copy()

    # Rotating the wind components, from the rotated model, to real world
    lat_lon_coord_system = iris.coord_systems.GeogCS(
        semi_major_axis=iris.fileformats.pp.EARTH_RADIUS
    )
    u_wind, v_wind = iris.analysis.cartography.rotate_winds(
        x_wind, y_wind, lat_lon_coord_system
    )

    # Getting the wind direction
    wind_dir.data = np.rad2deg(np.arctan2(u_wind.data, v_wind.data)) % 360

    # Changing units and names
    wind_dir.rename("wind_to_direction")
    wind_dir.attributes["short_name"] = short_name
    wind_dir.units = cf_units.Unit("degrees")
    cubelist.append(wind_dir)

    return cubelist


def regrid_aqum_on_haduk(cube_list_in, cube_target, method="Linear"):
    '''
    Assume AQUM fields are on native rotated-pole grid
    HadUK-Grid is provided on ONS 1km resolution TransverseMercator grid.
    '''

    print("\nStarting regriding AQUM cubes onto HadUK-GRID:")
    ycoord_tgt, xcoord_tgt = guess_coord_names(cube_target, ["Y", "X"])
    #print( cube_target.coord(ycoord_tgt) )

    cube_list_out = iris.cube.CubeList()

    if method in ["Linear", "Nearest"]:
        if method == "Linear":
            regrid_method=iris.analysis.Linear()
        elif method == "Nearest":
            regrid_method=iris.analysis.Nearest()
    else:
        raise ValueError(f"\nRegridding method {method} not available")

    for cube in cube_list_in:
        print(f" ... {cube.name()} | being regridded using" + \
              f" 'iris.analysis.{method}'")
        cube_regridded = cube.regrid(cube_target, regrid_method)
        cube_list_out.append(cube_regridded)

    print("\nEnd of regridding!")

    return(cube_list_out)


def extract_cube_area(cube=None, 
                      MINLON=None,
                      MAXLON=None,
                      MINLAT=None,
                      MAXLAT=None):
    """
    Function that extract a region from a cube that is not on 
    a regular latlon grid.

    :cube:      an iris cube
    :MINLON:    Minimum longitude to be used constraining the cube
    :MAXLON:    Maximum longitude to be used constraining the cube
    :MINLAT:    Minimum latitude to be used constraining the cube
    :MAXLAT:    Minimum latitude to be used constraining the cube
    """

    if MINLON is None:
        MINLON = -7.5
    if MAXLON is None:
        MAXLON = -1.5
    if MINLAT is None:
        MINLAT = 55.75
    if MAXLAT is None:
        MAXLAT = 59.5

    if cube is None:
        raise ValueError(f"No cube has been provided")
    
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError(f"{cube} is not an Iris cube")

    ycoord, xcoord = None, None
    ycoord, xcoord  = guess_coord_names(cube, ["Y", "X"])
    #print(cube.coord(xcoord))
    #print(cube.coord(ycoord))

    if None in [xcoord, ycoord]:
        raise ValueError(f"Can not unroate coordinates for {cube}")

    # define reg latlon coordinate system and pair of locations
    # to be derived in the cube coordinate space.
    lat_lon_coord_system = iris.coord_systems.GeogCS(
        semi_major_axis=iris.fileformats.pp.EARTH_RADIUS
    )
    lons = np.array([MINLON, MAXLON], dtype=float)
    lats = np.array([MINLAT, MAXLAT], dtype=float)

    # Calculate pair location in cube coordinates
    if xcoord != "longitude" and ycoord != "latitude":

        if xcoord == "grid_longitude" and ycoord == "grid_latitude":
            # rotated coord ...
            pole_lon = cube.coord(xcoord).coord_system.grid_north_pole_longitude
            pole_lat = cube.coord(ycoord).coord_system.grid_north_pole_latitude
            
            # Perform rotation
            rot_lons, rot_lats = rotate_pole(lons, lats, pole_lon, pole_lat)
            rot_lons = rot_lons + 360
            
            print(f"\nReg Lat {lats} converted to " + \
                  f"rotated lat coord: {rot_lats}")
            print(f"Reg Lon {lons} converted to " + \
                  f"rotated lon coord: {rot_lons}")

            lat_constraint = iris.Constraint(
                grid_latitude=lambda cell: rot_lats[0] < cell < rot_lats[1]
                )
            lon_constraint = iris.Constraint(
                grid_longitude=lambda cell: rot_lons[0] < cell < rot_lons[1]
                )

            cube = cube.extract(lon_constraint & lat_constraint)

        elif (xcoord == "projection_x_coordinate" and 
              ycoord == "projection_y_coordinate"):
            # Other coordinate system (note this may work for x/ycoords other than
            # those considered here
            ll_crs = lat_lon_coord_system.as_cartopy_crs()
            cube_crs = cube.coord(xcoord).coord_system.as_cartopy_crs()

            # Convert to lat/lon points
            cube_lonlats = cube_crs.transform_points(ll_crs, lons, lats)
            cube_lons = cube_lonlats[:, 0]
            cube_lats = cube_lonlats[:, 1]

            print(f"\nReg Lat {lats} converted to " + \
                  f"projection_y_coordinate: {cube_lats}")
            print(f"Reg Lon {lons} converted to " + \
                  f"projection_x_coordinate: {cube_lons}")
            
            lat_constraint = iris.Constraint(
                projection_y_coordinate=lambda cell: cube_lats[0] < cell < cube_lats[1]
                )
            lon_constraint = iris.Constraint(
                projection_x_coordinate=lambda cell: cube_lons[0] < cell < cube_lons[1]
                )

            cube = cube.extract(lon_constraint & lat_constraint)

    return(cube)


def add_reg_latlon_to_cube(cube):
    """
    // Unused // 
    Add regular lat-lon auxillary coordinates to a cube.
    Function only works for input cube on rotated coord.

    :cube:  a iris cube (e.g. an AQUM cube with rotated pole coordinates)
    """

    if not isinstance(cube, iris.cube.Cube):
        raise ValueError(f"{cube} is not an Iris cube")

    ycoord, xcoord  = guess_coord_names(cube, ["Y", "X"])
    print(f'\n{xcoord}:')
    print(cube.coord(xcoord))
    print(f'\n{ycoord}:')
    print(cube.coord(ycoord))

    if None in [xcoord, ycoord]:
        raise ValueError(f"Can not unroate coordinates for {cube}")
        
    #coordinate system for lat lon grid to be added:
    # Note: iris.fileformats.pp.EARTH_RADIUS = 6371229.0  
    # lat_lon_coord_system = iris.coord_systems.GeogCS(6371229.0)
    lat_lon_coord_system = iris.coord_systems.GeogCS(
        semi_major_axis=iris.fileformats.pp.EARTH_RADIUS
    ) 

    if xcoord != "longitude" and ycoord != "latitude":

        if xcoord == "grid_longitude" and ycoord == "grid_latitude":
            #rotated coord ...

            #1. Definitions of coordinate systems.
            src_crs = cube.coord('grid_latitude').coord_system.as_cartopy_crs()
            tgt_cs = lat_lon_coord_system
            tgt_crs = tgt_cs.as_cartopy_crs()

            #2. Convert the 1d arrays to two 2d arrays
            cube_lon_pts = cube.coord(xcoord).points
            cube_lat_pts = cube.coord(ycoord).points

            cube_lon, cube_lat = np.meshgrid(cube_lon_pts, cube_lat_pts)
            #print("\ncube_lon")
            #print(cube_lon)
            #print(cube_lon.shape)
            #print("\ncube_lat")
            #print(cube_lat)
            #print(cube_lat.shape)

            #3. Calculate unrotated coordinates
            lonlats = tgt_crs.transform_points(src_crs, cube_lon, cube_lat)
            print("\nunrotated coordinates lonlats shape:")
            print(lonlats.shape)
            print("\nlonlats[...,0]")
            print(lonlats[...,0])
            print(lonlats[...,0].shape)
            print("\nlonlats[...,1]")
            print(lonlats[...,1])
            print(lonlats[...,1].shape)
                
            #4. Construct the iris coords 
            latitude = iris.coords.AuxCoord(
                lonlats[...,1], 
                standard_name='latitude', 
                units="degrees",
                coord_system=tgt_cs) 

            longitude = iris.coords.AuxCoord(
                lonlats[...,0],
                standard_name='longitude',
                units="degrees",
                coord_system=tgt_cs)
            
            print("\n")
            print(longitude)
            print(latitude)
            cube.add_aux_coord(longitude, (0,1))
            cube.add_aux_coord(latitude, (0,1))
            print(cube)

        elif (
            xcoord == "projection_x_coordinate" and \
            ycoord == "projection_y_coordinate"
        ):
            # TODO: Other coordinate system 
            pass

    return(cube)


#   ##############################
#               START
#   ##############################

if __name__ == "__main__":

    print("\n------------------------------------" + \
          "\nStart of regrid_AQUM_to_HadUK.py:" + \
          "\n------------------------------------\n"
         )

    tracemalloc.start()
    stall = time.time()

    # Select if processed cubes to be exported to netcdf files
    l_export_to_nc = True

    # Derive wind direction (requires x_wind and y_wind to be requested)
    l_calc_wind_direction = True

    # Extract a subset of the domain to speed-up interpolation
    l_extract_subdomain = True

    # Select if Time averaging of input cubes 
    # is done before or after regridding.
    l_time_avg_before_rg = True
    l_time_avg_before_rg = False

    if l_time_avg_before_rg:
        l_time_avg_after_rg = False
    else:
        l_time_avg_after_rg = True

    # Read AQUM data path and filename from parsed arguments.
    # Load default file if no argument provided.
    aqum_path = None
    aqum_filename = None

    if len(sys.argv) > 1:
        if os.path.exists(sys.argv[1]):
            aqum_path = sys.argv[1]
        else:
            raise ValueError(f"\n {sys.argv[1]} does not exist")

    if len(sys.argv) > 2:
        if os.path.exists(aqum_path + '/' + sys.argv[2]):
            aqum_filename = sys.argv[2]
        else:
            raise ValueError(f"\n {sys.argv[2]} does not exist")

    if (aqum_path and aqum_filename) is None:
        # Use precomputed monthly/annual means
        #aqum_path = '/data/users/fmalavel/AQUM/PROJECT/Corrosion_Scotland/' + \
        #        'verification'
        aqum_path = './AQUM_time_mean_outputs'
        aqum_filename = '2022met_2019ems_gridded_cube_list_annual_mean.nc'
        aqum_filename = '2022met_2019ems_gridded_cube_list_monthly_mean.nc'

    aqum_input_file = aqum_path + '/' + aqum_filename

    # Derive Met Year to be processed based on filename
    pattern = '^[0-9][0-9][0-9][0-9]met'
    match = re.match(pattern, aqum_filename)
    met_year = int(match.group()[0:4])

    pdt1 = PartialDateTime(year=met_year, month=1, day=1)
    pdt2 = PartialDateTime(year=met_year, month=12, day=31)
    year_constraint = iris.Constraint(time=
                                      lambda cell: pdt1 <= cell.point < pdt2)

    # Note on wind components: section 0 winds are on section 0 wind vectors 
    # are output on different grids.
    # Use section 3 (10m winds) or section 15 wind diagnostics which are on
    # B-GRID to derive wind speed and wind direction.

    # List of cube short names to be loaded (consider passing as argument)
    # short_name_list = None
    short_name_list = ["T_surf", "RHw", "total_precip",
                       "u_10m", "v_10m",
                       "SO2", "SO2_DryDep", "SS_WetDep",
                       "O3", "O3_DryDep", 
                       "SS_Conc", "SS_DryDep", "SS_WetDep"]

    # List of cube short_names to be regridded (consider passing as argument)
    #short_name_list = None
    rg_short_name_list = ["T_surf", "RHw", 
                          "u_10m", "v_10m",
                          "SO2", "SO2_DryDep", "SS_WetDep",
                          "O3", "O3_DryDep",
                          "SS_Conc", "SS_DryDep", "SS_WetDep"
                         ]

    #rg_short_name_list = ["T_surf", "RHw", "u_10m", "v_10m"]

    print('aqum_path:', aqum_path)
    print('aqum_filename:', aqum_filename)
    print('met_year:', met_year)


#   #############################
#      Load AQUM fields
#   #############################

    print("\n------------------\nLoading AQUM fiedls:\n------------------")
    st = time.time()

    if short_name_list:
        #load cubes matching predefined short names
        short_constraint = iris.AttributeConstraint(
            short_name=lambda c: c in short_name_list
        )
        input_cube_list = iris.load(aqum_input_file, short_constraint)
    else:
        #load everything in the file
        input_cube_list = iris.load(aqum_input_file)
    
    if l_calc_wind_direction:
        wind_direction(input_cube_list, short_name="wind_direction")
        wind_speed(input_cube_list, short_name="wind_speed")
        rg_short_name_list += ["wind_direction", "wind_speed"]

    print(f"\n{input_cube_list}")
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

    haduk_data_path = '/data/users/fmalavel/AQUM/PROJECT/Corrosion_Scotland/HadUK-GRID/'

    # Surface temp
    var = 'tas'
    haduk_1km_annual_mean = haduk_data_path + var + '_hadukgrid_uk_1km_ann_' + \
            str(met_year) + '01' + '-' + str(met_year) + '12.nc'
    haduk_cube_tas = iris.load_cube(haduk_1km_annual_mean, year_constraint)

    # Surface temp
    var = 'hurs'
    haduk_1km_annual_mean = haduk_data_path + var + '_hadukgrid_uk_1km_ann_' + \
            str(met_year) + '01' + '-' + str(met_year) + '12.nc'
    haduk_cube_hurs = iris.load_cube(haduk_1km_annual_mean, year_constraint)

    # 10m surface wind
    var = 'sfcWind'
    haduk_1km_annual_mean = haduk_data_path + var + '_hadukgrid_uk_1km_ann_' + \
            str(met_year) + '01' + '-' + str(met_year) + '12.nc'
    haduk_cube_sfcWind = iris.load_cube(haduk_1km_annual_mean, year_constraint)

    #take a copy of any of the HadGRID loaded cubes for AQUM fields interpolation 
    haduk_cube = haduk_cube_sfcWind.copy()

    print(haduk_cube)
    et = time.time()
    elapsed_time = et - st
    print('\nExecution time:', elapsed_time, 'seconds')


#   #############################
#    Decrease target domain size
#   #############################

    if l_extract_subdomain:

        print("\n------------------\n" + \
              "Decrease cubes size:" + \
              "\n------------------\n")

        st = time.time()

        haduk_cube_tas = extract_cube_area(haduk_cube_tas)
        haduk_cube_hurs = extract_cube_area(haduk_cube_hurs)
        haduk_cube_sfcWind = extract_cube_area(haduk_cube_sfcWind)
        haduk_cube = extract_cube_area(haduk_cube)

        print(f"\nhaduk_cube: {haduk_cube}")

        et = time.time()
        elapsed_time = et - st
        print('\nExecution time:', elapsed_time, 'seconds')


#   #############################
#        Calc Time Average
#        before regridding
#   #############################

    if l_time_avg_before_rg:

        print("\n------------------\n" + \
              "Average over Time coord before regridding" + \
              "\n------------------\n")
        st = time.time()

        input_cube_list_avg = iris.cube.CubeList()

        for i, cube in enumerate(input_cube_list):
            time_coord = guess_coord_names(cube, ["T"])
            if time_coord:
                cube_avg = cube.copy().collapsed('time', iris.analysis.MEAN)
                input_cube_list_avg.append(cube_avg)
        
        # remove cubes from original cubelist ...
        input_cube_list.clear()
        
        # ... repopulate with time averaged cubes
        for cube in input_cube_list_avg:
            input_cube_list.append(cube)
            
        print(f"{input_cube_list}")
        et = time.time()
        elapsed_time = et - st
        print('\nExecution time:', elapsed_time, 'seconds')


#   #############################
#               REGRID
#   #############################
           
    print("\n------------------\n" + \
          "REGRIDDING:" + \
          "\n------------------")
    st = time.time()

    # Copy AQUM cubes to be regridded into a new cube list.
    aqum_native_cubes = iris.cube.CubeList()

    if rg_short_name_list is None:
        rg_short_name_list = short_name_list

    for cube in input_cube_list:
        if cube.attributes["short_name"] in rg_short_name_list:
            if cube.units == "K":
                cube.convert_units("degC")
            aqum_native_cubes.append(cube)

    print("\nAQUM input fields:")
    print(aqum_native_cubes)

    aqum_regridded_cubes = regrid_aqum_on_haduk(aqum_native_cubes, haduk_cube)

    print("\nAQUM regridded fields:")
    print(aqum_regridded_cubes)
    et = time.time()
    elapsed_time = et - st
    print('\nExecution time:', elapsed_time, 'seconds')


#   #############################
#        Calc Time Average
#        after regridding
#   #############################

    if l_time_avg_after_rg:

        print("\n------------------\n" + \
              "Average over Time coord after regridding" + \
              "\n------------------\n")

        # 1 - AQUM cubelist on native grid
        st = time.time()
        aqum_native_cubes_avg = iris.cube.CubeList()

        for i, cube in enumerate(aqum_native_cubes):
            print(f"{i} collapsing {cube.summary(True)}")
            time_coord = guess_coord_names(cube, ["T"])
            if time_coord:
                cube_avg = cube.copy().collapsed('time', iris.analysis.MEAN)
                aqum_native_cubes_avg.append(cube_avg)

        # remove cubes from original cubelist ...
        aqum_native_cubes.clear()

        # ... repopulate with time averaged cubes
        for i, cube in enumerate(aqum_native_cubes_avg):
            aqum_native_cubes.append(cube)
        et = time.time()
        elapsed_time = et - st
        print('\nExecution time 1/2:', elapsed_time, 'seconds\n')


        # 2 - AQUM cubelist on HadUK-GRID grid
        st = time.time()
        aqum_regridded_cubes_avg = iris.cube.CubeList()

        for j, cube in enumerate(aqum_regridded_cubes):
            print(f"{j} collapsing {cube.summary(True)}")
            time_coord = guess_coord_names(cube, ["T"])
            if time_coord:
                cube_avg = cube.copy().collapsed('time', iris.analysis.MEAN)
                aqum_regridded_cubes_avg.append(cube_avg)

        # remove cubes from original cubelist ...
        aqum_regridded_cubes.clear()

        # ... repopulate with time averaged cubes
        for j, cube in enumerate(aqum_regridded_cubes_avg):
            aqum_regridded_cubes.append(cube)
        et = time.time()
        elapsed_time = et - st
        print('\nExecution time 2/2:', elapsed_time, 'seconds')


#   #############################
#               PLOTS
#   #############################

    print("\n------------------\nCreating Plots:\n------------------")
    st = time.time()

    if l_extract_subdomain:
        prefix = "subdomain_" + str(met_year) + "met_"
    else:
        prefix = str(met_year) + "met_"

    # Loop over a cubelist of aqum fields
    for cube_rg in aqum_regridded_cubes:

        short = cube_rg.attributes["short_name"]
        con_short = iris.AttributeConstraint(short_name=short)

        native_aqum = aqum_native_cubes.extract_cube(con_short)
        regrid_aqum = aqum_regridded_cubes.extract_cube(con_short)

        cb_min=None
        cb_max=None
        colormap=None

        if short == "T_surf":
            haduk_cube = haduk_cube_tas
            cb_min=4
            cb_max=12
            colormap="Spectral"
        if short == "RHw":
            haduk_cube = haduk_cube_hurs
            cb_min=70
            cb_max=90
            colormap="GnBu"
        if short in ["wind_speed"]:
            haduk_cube = haduk_cube_sfcWind
            cb_min=1
            cb_max=9
            colormap="plasma"
        if short in ["u_10m", "v_10m"]:
            haduk_cube = haduk_cube_sfcWind
            cb_min=-5
            cb_max=5
            colormap="bwr"
        if short in ["wind_direction"]:
            haduk_cube = haduk_cube_sfcWind
            cb_min=0
            cb_max=360
            colormap=None

        plot_regridded_fields(native_aqum, regrid_aqum, haduk_cube, 
                              cb_min=cb_min, 
                              cb_max=cb_max, 
                              colormap=colormap,
                              title_str = short,
                              plot_filename_prefix = prefix
                             )
    et = time.time()
    elapsed_time = et - st
    print('\nExecution time:', elapsed_time, 'seconds')


#   #############################
#   Saving regridded Fields as nc
#   #############################

    if l_export_to_nc:

        print("\n------------------" + \
              "\nExporting to netcdf:" + \
              "\n------------------\n")
        st = time.time()
    
        if l_time_avg_before_rg:
            nc_output_nat = str(met_year) + 'met_AQUM_native_' + \
                    'time_avg_before.nc'
            nc_output_reg = str(met_year) + 'met_AQUM_regridded_' + \
                    'time_avg_before.nc'
        elif l_time_avg_after_rg:
            nc_output_nat = str(met_year) + 'met_AQUM_native_' + \
                    'time_avg_after.nc'
            nc_output_reg = str(met_year) + 'met_AQUM_regridded_' + \
                    'time_avg_after.nc'

        if l_extract_subdomain:
            nc_output_nat = "subdomain_" + nc_output_nat
            nc_output_reg = "subdomain_" + nc_output_reg
        else:
            nc_output_nat = "fulldomain_" + nc_output_nat
            nc_output_reg = "fulldomain_" + nc_output_reg

        iris.save(aqum_native_cubes, output_path + nc_output_nat)
        iris.save(aqum_regridded_cubes, output_path + nc_output_reg)

        print(f"Native cubes exported to {nc_output_nat} ")
        print(f"Reggrided cube exported to {nc_output_reg} ")

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

    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')

    print("\n[ Top 10 ]")
    for stat in top_stats[:10]:
        print(stat)
