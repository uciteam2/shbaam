#!/usr/bin/env python2

import rtree
import math
import shapely.prepared
import shapely.geometry
import fiona
import sys
import netCDF4
import datetime
import csv
import subprocess
import os
import numpy as np

"""
Computes total terrestrial water storage anomaly timeseries.

params:
    longitudes {list} list containing longitude values from the original input netCDF4 file
    latitudes {list} list containing latitude values from the original input netCDF4 file
    swe_average {list} list containing SWE Averages calculated over the total surface area for a particular grid cell
    surface_area {list} list containing the surface area of each grid cell
    input_netCDF4 {file} file object which is the representation of the input netCDF4 file

return:
    {list} a list containing the total swe for a particular time (month in this case)
"""


def water_storage_timeseries(longitudes, latitudes, swe_averages, surface_areas, input_netCDF4):
    print('Compute total terrestrial water storage anomaly timeseries')

    # compute total surface area
    total_surface_area = 0
    for grid_cell in range(len(swe_averages)):
        total_surface_area += surface_areas[grid_cell]

    # init timeseries
    timeseries = []
    
    # grab all time dimensions for the timeseries, there should be one for each month we are measuring
    times = len(input_netCDF4.dimensions['time'])
    for time in range(times):
        total_swe_anomaly = 0

        # iterate through each grid cell, use swe_average since it has an average for each grid cell
        for grid_cell in range(len(swe_averages)):

            # init values to be used in calculations
            lon = longitudes[grid_cell]
            lat = latitudes[grid_cell]
            surface_area = surface_areas[grid_cell]
            swe_average = swe_averages[grid_cell]

            # in original code, he did something with the mask, but I don't think we have to do this? line 363-366 {?}

            # add to total surface area

            # get original SWE for the grid cell at the time for this loop
            swe_total = input_netCDF4.variables['SWE'][time, lat, lon]

            # calculate the anomaly for this grid cell
            anomaly = swe_total - swe_average

            # multiply by the surface area to eliminate m^2
            anomaly *= surface_area

            # add to the running total for this period in time
            total_swe_anomaly += anomaly

            # want millimeters averaged over the entire area
            # compute mass of volume of each grid cell

            # add them up and divide over total area

        # get average
        timeseries.append(total_swe_anomaly / total_surface_area)

    return timeseries


"""
Creates a series of datetime strings for each of the time dimensions in the input netCDF4 file. Will use later when creating the output CSV file
NOTE: this relies on my conc file to work accurately, because it gets access to the conc time variable and takes the date from there
params:
    input_netCDF4 {file} file object which is the representation of the input netCDF4 file

return:
    timestrings {list} a list of datetime strings with each date representing a different month in the netCDF4 file
"""


def create_timestrings(input_netCDF4):
    print('Determining datestrings...')
    # grab initial date from conc file. this uses the units from the conc file
    timestring = ' '.join(input_netCDF4.variables['time'].units.split(' ')[2:])
    origin_date = datetime.datetime.strptime(timestring, '%Y-%m-%d %H:%M:%S')

    # init return array
    timestrings = [origin_date]

    for month_delta in range(1, len(input_netCDF4.variables['time'])):
        # adds a new month to the timestrings, passing in the previous month to add to

        timestrings.append(add_month(timestrings[month_delta - 1]))

    return [date.strftime('%m/%d/%Y') for date in timestrings]


"""
Helper function for creating datestrings. Increments the month by one. This can't be done using timedelta
because a month is not a uniform measure.

params:
    date {datetime} a datetime object to increment by a month

return:
    {datetime} an object which is a month ahead of the input datetime
"""


def add_month(date):
    new_year = date.year
    new_month = date.month + 1
    if new_month > 12:
        new_month = 1
        new_year += 1
    new_day = date.day
    return datetime.datetime(new_year, new_month, new_day)


"""
Creates a csv with the SWE average for each month
params:
    timestrings {list} a list of strings in MM/DD/YYYY format. One for each time dimension in the input netCDF4
    swe_sums {list} a list of swe_sums derived from the water_storage_timeseries function
    output_file {str} file string of the output file. It should have already been tested before input to this function
"""


def create_csv(timestrings, swe_sums, output_csv):
    print('Creating CSV and writing average SWE to file')

    # only write if we have the same number of swe sums and timestrings
    if(len(timestrings) == len(swe_sums)):

        # Open csv at the location that was given in the args
        with open(output_csv, 'wb') as csvfile:

            # create a csv writer
            csvwriter = csv.writer(csvfile, dialect='excel')

            # write the headers
            csvwriter.writerow(['Month', 'SWE Sum'])

            # loop through the results and write to the csv file
            for i in range(len(timestrings)):
                csvwriter.writerow([timestrings[i], swe_sums[i]])

    else:
        # errors out if they aren't the same length
        print('ERROR: Timestrings length is not the same as SWE sums')
        print('Timestrings length: ', len(timestrings))
        print('SWE Sums length: ', len(swe_sums))


"""
Retrieves the fillvalue from the original input netCDF4 file. If there is no fill value specified,
then it returns the default

params:
    input_netCDF4 {netCDF4.Dataseet} Dataset object from the netCDF4 library

returns:
    the fillvalue to be used when creating the netCDF4 output
"""

def get_fillvalue(input_netCDF4):
    fillvalue = netCDF4.default_fillvals['f4']
    if 'RUNSF' in input_netCDF4.variables:
        runsf = input_netCDF4.variables['RUNSF']
        if '_FillValue' in runsf.ncattrs():
            fillvalue = runsf._FillValue
            print('The fill value for RUNSF is: ', str(fillvalue))
        else:
            fillvalue = None

    return fillvalue


"""
Creates an output netCDF4 file with...

params:
    input_filepath {str} a filepath given as one of the arguments to run the script. The location of the input nc4 file
    output_filepath {str} a filepath given as one of the arguments to run the script. The location of the output nc4 file
    fillvalue {str} the fillvalue to use when creating the masked arrays. Retrieve using the get_fillvalue function defined earlier
    longitudes {list} longitudes taken from the original input nc4 file
    latitudes {list}  latitudes taken from the original input nc4 file
    swe_average {list} list containing SWE Averages calculated over the total surface area for a particular grid cell
    timeseries {list} a list containing the total swe for a particular time (month in this case)
"""


def create_output_netCDF4(input_netCDF4, output_filepath, fillvalue, gld_lon, gld_lat, intersect_lon, intersect_lat, swe_averages, timeseries):
    print('creating output netCDF4 file...')

    # create the output nc4 file from the filepath supplied in the args
    output_netCDF4 = netCDF4.Dataset(output_filepath, 'w', format='NETCDF3_CLASSIC')
    # input_netCDF4 = netCDF4.Dataset(input_filepath, 'r')
    

    # create dimensions
    create_dimensions(output_netCDF4, gld_lon, gld_lat)

    # create the variables then return them to use later in creating variable attributes
    variables = create_variables(output_netCDF4, fillvalue)

    # create global attributes for the nc4 file
    create_global_attributes(output_netCDF4)

    # create variable attributes from the original nc4 file
    copy_variable_attributes(input_netCDF4, output_netCDF4)

    # modify the CRS variable attributes (no idea why this is here)
    modify_crs_attribute(output_netCDF4, variables)

    # populate static data
    print('populating static data')

    lon = variables['lon']
    lon[:] = input_netCDF4.variables['lon'][:]
    # variables['lon'][:] = input_netCDF4.variables['lon'][:]
    lat = variables['lat']
    lat[:] = input_netCDF4.variables['lat'][:]

    # variables['lat'][:] = input_netCDF4.variables['lat'][:]

    populate_dynamic_data(output_netCDF4, input_netCDF4, intersect_lon, intersect_lat, swe_averages, timeseries)

    # close files
    output_netCDF4.close()
    input_netCDF4.close()


"""
creates the dimensions for the output netCDF4 file.

params:
    output_netCDF4 {netCDF4.Dataset} the dataset which was created
    longitudes {list} longitudes taken from the original input nc4 file
    latitudes {list}  latitudes taken from the original input nc4 file
"""


def create_dimensions(output_netCDF4, longitudes, latitudes):
    print('creating dimensions...')
    time = output_netCDF4.createDimension('time', None)
    lat = output_netCDF4.createDimension('lat', len(latitudes))
    lon = output_netCDF4.createDimension('lon', len(longitudes))
    nv = output_netCDF4.createDimension('nv', 2)  # no idea why this is here


"""
Creates the variables for the output nc4 file
params:
    output_netCDF4 {netCDF4.Dataset} the dataset which was created
    fillvalue {str} the fillvalue to be used in the swe variable

returns:
    {dict} the new variables objects of the output dataset in {'name': variable} format
"""


def create_variables(output_netCDF4, fillvalue):
    print('creating variables')
    time = output_netCDF4.createVariable('time', 'i4', ('time',))
    time_bands = output_netCDF4.createVariable(
        'time_bands', 'i4', ('time', 'nv'))
    lat = output_netCDF4.createVariable('lat', 'f4', ('lat',))
    lon = output_netCDF4.createVariable('lon', 'f4', ('lon',))
    swe = output_netCDF4.createVariable('swe', 'f4', ('time', 'lat', 'lon'), fill_value=fillvalue)

    # not sure what crs is
    crs = output_netCDF4.createVariable('crs', 'i4')

    return {
        'time': time,
        'time_bands': time_bands,
        'lat': lat,
        'lon': lon,
        'swe': swe,
        'crs': crs
    }


"""
Creates the global variables for the output nc4 file
params:
    output_netCDF4 {netCDF4.Dataset} the dataset which was created
"""


def create_global_attributes(output_netCDF4):
    # Get current time without microseconds
    dt = datetime.datetime.utcnow()
    dt = dt.replace(microsecond=0)

    # grab the current version of shbaam
    version = subprocess.Popen('bash ../version.sh', stdout=subprocess.PIPE, shell=True).communicate()
    version = version[0].rstrip()

    # create global attributes
    output_netCDF4.Conventions = 'CF-1.6'
    output_netCDF4.title = ''
    output_netCDF4.institution = ''
    output_netCDF4.source = 'SHBAAM: ' + version + \
        'GRACE: '
    output_netCDF4.history = 'date created: ' + dt.isoformat() + '+00:00'
    output_netCDF4.references = 'https://github.com/c-h-david/shbaam/'
    output_netCDF4.comment = ''
    output_netCDF4.featureType = 'timeSeries'


"""
Copies the attributes from the input netCDF4 file provided in the program args over to
the newly created output file.

params:
    input_filepath {str} the filepath to the input nc4 file
    destination {netCDF4.Dataset} the dataset which was created
"""


def copy_variable_attributes(source, destination):
    pass
    #TODO MATCH CEDRIC CODE
    # this should be revised to only grab feasable attributes


"""
These are for the WGS84 spheroid

params:
    output_netCDF4 {netCDF4.Dataset} the dataset which was created
    output_variables {dict} the variables objects of the output dataset in {'name': variable} format
"""


def modify_crs_attribute(output_netCDF4, output_variables):
    output_variables['swe'].grid_mappings = 'crs'
    output_variables['crs'].grid_mapping_name = 'latitude_longitude'
    output_variables['crs'].semi_major_axis = '6378137'
    output_variables['crs'].inverse_flattening = '298.257223563'


"""
populates the output file with the relavent data
"""


def populate_dynamic_data(output_netCDF4, input_netCDF4, longitudes, latitudes, swe_averages, timeseries):
    print('populate dynamic data...')

    for grid_cell in range(len(swe_averages)):
        lon = longitudes[grid_cell]
        lat = latitudes[grid_cell]
        swe_average = swe_averages[grid_cell]

        for time in range(len(timeseries)):
            try:
                output_netCDF4.variables['swe'][time, lat, lon] = input_netCDF4.variables['SWE'][time, lat, lon] - swe_average
            except Exception as e:
                print(e)

    # why do we do this?
    output_netCDF4.variables['time'][:] = input_netCDF4.variables['time'][:]
# =================================================Brian==========================================================
# ================================================================================================================


# =================================================Michelle======================================================
# ================================================================================================================
'''
Perform calculations on the grid cells: find long-term means and surface areas for each grid cell
	- total_num_cells:	the number of grid cells we are interested in, used for looping
	- grid_lats:		array of latitudes for each grid cell
	- grid_lons: 		array of longitudes (corresponding to the lats) for each grid cell
	- times: 		the total number of "times" from the time dimension of the netCDF file
	- cdf_file: 		the cdf_file itself - used to access the SWE values
	- lat_interval:		the size of the latitude interval from the netCDF file
	- lon_interval:		the size of the longitude interval from the netCDF file

Returns: Tuple containing 2 arrays: (time_averages, surface_areas)
'''


def grid_calculations(total_num_cells, grid_lats, grid_lons, actual_lats, times, cdf_file, lat_interval, lon_interval):
	# construct empty arrays to fill with values
	time_averages = [0] * total_num_cells
	surface_areas = [0] * total_num_cells

	# iterate through the lats and lons for each grid cell
	for grid in range(total_num_cells):
		# retrieve lat and lon
		lat = grid_lats[grid]  # the index of latitude for the current cell
		lon = grid_lons[grid]  # the index of longitude for the current cell

		# get the surface area for this cell and add it to the surface_areas array
		SA = 6371000 * math.radians(lat_interval) * 6371000 * math.radians(lon_interval) * math.cos(math.radians(actual_lats[lat]))  # make a global var for 6371000
		surface_areas[grid] = SA

		# iterate through the grid cell at each time in the netCDf file's time dimension
		for time in range(times):
			# create a running total of SWE values in the time_averages array
		    time_averages[grid] += cdf_file.variables['SWE'][time, lat, lon]

	# divide the values in the time_averages array by times to get the average
	time_averages = [x/times for x in time_averages]

	# return 2 arrays: surface areas and time_averages
	return (time_averages, surface_areas)

# =================================================/Michelle======================================================
# ================================================================================================================

# =================================================Caleb==========================================================
# ================================================================================================================

# Purpose:


# *******************************************************************************
# Import Python modules
# *******************************************************************************


def check_command_line_arg():
    # Checks the length of arguements and if input files exist
    IS_arg = len(sys.argv)
    if IS_arg != 6:
        print('ERROR - 5 and only 5 arguments can be used')
        raise SystemExit(22)

    for shb_file in sys.argv[1:3]:
        try:
            with open(shb_file) as file:
                pass
        except IOError as e:
                print('ERROR - Unable to open ' + shb_file)
                raise SystemExit(22)

        print('[+] Command Line Entered Properly')


def readPolygonShpFile(polyFile):
    print('Read polygon shapefile')
    polygon_file = fiona.open(polyFile, 'r')
    polygon_features = len(polygon_file)
    print(' - The number of polygone features is: ' + str(polygon_features))
    print('[+] Returned Polygon File Properly')
    return polygon_file


def createShapeFile(gld_dim_lat_length, gld_dim_lon_length, gld_lon_dimension_array, gld_lat_length_dimension_array, polygonShapeFile, output_pnt_shp):
    shapeFile_Driver = polygonShapeFile.driver
    shapeFile_Driver_Copy = shapeFile_Driver

    shapeFile_crs = polygonShapeFile.crs
    shapeFilePoint_crs = shapeFile_crs.copy()

	# Kept the same names, just in case
    shapeFile_schema = {'geometry': 'Point','properties': {'lon': 'int:4','lat': 'int:4'}}

    with fiona.open(output_pnt_shp, 'w', driver=shapeFile_Driver_Copy,crs=shapeFilePoint_crs, schema=shapeFile_schema) as pointFile:
        for lon in range(gld_dim_lon_length):
            gld_lon = gld_lon_dimension_array[lon]
            if (gld_lon > 180):
            # Shifts GLD range [0:360] to [-180:180]
                gld_lon -= 360
            for lat in range(gld_dim_lat_length):
                gld_lat = gld_lat_length_dimension_array[lat]
                shapeFilePoint_Prepared = {'lon': lon,
                            'lat': lat}
                shapeFilePoint_geometry = shapely.geometry.mapping(shapely.geometry.Point((gld_lon, gld_lat)))
                pointFile.write({
                    'properties': shapeFilePoint_Prepared,
                    'geometry': shapeFilePoint_geometry,
                    })

    print('[+] New ShapeFile Created')


def createSpatialIndex(pointFile):
    print('Create spatial index for the bounds of each point feature')
    index = rtree.index.Index()

    for feature in pointFile:
        feature_id = int(feature['id'])
    	# the first argument of index.insert has to be 'int', not 'long' or 'str'
        shape = shapely.geometry.shape(feature['geometry'])
        index.insert(feature_id, shape.bounds)
    	# creates an index between the feature ID and the bounds of that feature

    print(' - Spatial index created')
    return index


def find_intersection(polygon, index, points):
    intersect_tot = 0
    intersect_lon = []
    intersect_lat = []

    for area in polygon:
        shape_geo = shapely.geometry.shape(area['geometry'])
        shape_prep = shapely.prepared.prep(shape_geo)
	    # a 'prepared' geometry allows for faster processing after
        for point_id in [int(x) for x in list(index.intersection(shape_geo.bounds))]:
            shape_feature = points[point_id]
            shape_file = shapely.geometry.shape(shape_feature['geometry'])
            if shape_prep.contains(shape_file):
                # CHANGE NAME !!!????
                JS_dom_lon = shape_feature['properties']['lon']
                JS_dom_lat = shape_feature['properties']['lat']
                intersect_lon.append(JS_dom_lon)
                intersect_lat.append(JS_dom_lat)
                intersect_tot += 1

    print(' - The number of grid cells found is: '+str(intersect_tot))
    return (intersect_tot, intersect_lon, intersect_lat)


if __name__ == '__main__':
    check_command_line_arg()

    input_gld_nc4 = sys.argv[1]  # shb_grc_ncf; Concatenated File
    input_pol_shp = sys.argv[2]  # shb_pol_ncf; Polygon File
    output_pnt_shp = sys.argv[3]  # shb_pnt_shp; Point File
    output_swe_csv = sys.argv[4]  # shb_wsa_csv;
    output_swe_ncf = sys.argv[5]  # shb_wsa_ncf

    print('Read GLD netCDF file')
    f = netCDF4.Dataset(input_gld_nc4, 'r')

    # Get Dimension Sizes
    number_of_lon = len(f.dimensions['lon'])  # IS_grc_lon
    print(' - The number of longitudes is: '+str(number_of_lon))
    number_of_lat = len(f.dimensions['lat'])  # IS_grc_lat
    print(' - The number of latitudes is: '+str(number_of_lat))
    num_of_time_steps = len(f.dimensions['time'])  # IS_grc_time
    print(' - The number of time steps is: '+str(num_of_time_steps))

    # Value of Dimension Arrays
    gld_lon = f.variables['lon']  # ZV_grc_lon
    gld_lat = f.variables['lat']  # ZV_grc_lat
    gld_time = f.variables['time']  # ZV_grc_time

    # Get Interval Sizes
    gld_lon_interval_size = abs(gld_lon[1] - gld_lon[0])
    print(' - The interval size for longitudes is: ' + str(gld_lon_interval_size))
    gld_lat_interval_size = abs(gld_lat[1]-gld_lat[0])
    print(' - The interval size for latitudes is: ' + str(gld_lat_interval_size))
    if len(gld_time) > 1:
	    gld_time_interval_size = abs(gld_time[1] - gld_time[0])
    else:
	    gld_time_interval_size=0

    # Get Fill Values
    # gld_fill=netCDF4.default.fillvals['f4']
    gld_fill = get_fillvalue(f)

    if 'RUNSF' in f.variables:
	    gld_fill=var._FillValue
	    print(' - The fill value for RUNSF is: '+str(gld_fill))
    else:
	    gld_fill=None

    print('[+] Variables are set up properly')  # Can Delete


    polyShapeFile = readPolygonShpFile(input_pol_shp)  # shb_pol_lay
    createShapeFile(number_of_lat, number_of_lon, gld_lon, gld_lat, polyShapeFile, output_pnt_shp)

    point_features=fiona.open(output_pnt_shp, 'r')  # shb_pnt_lay
    index=createSpatialIndex(point_features)
    intersect_tot, intersect_lon, intersect_lat = find_intersection(polyShapeFile, index, point_features)

    time_averages, surface_areas = grid_calculations(intersect_tot, intersect_lat, intersect_lon, gld_lat, num_of_time_steps, f, gld_lat_interval_size, gld_lon_interval_size)

    swe_time_series = water_storage_timeseries(intersect_lon, intersect_lat, time_averages, surface_areas, f)

    print('SWE timeseries average: {}'.format(np.average(swe_time_series)))
    print('SWE timeseries min: {}'.format(np.min(swe_time_series)))
    print('SWE timeseries max: {}'.format(np.max(swe_time_series)))

    timestrings = create_timestrings(f)

    create_csv(timestrings, swe_time_series, output_swe_csv)
    fillvalue = get_fillvalue(f)

    create_output_netCDF4(f, output_swe_ncf, fillvalue, gld_lon, gld_lat, intersect_lon, intersect_lat, time_averages, swe_time_series)

    print('[+] Script Completed')
