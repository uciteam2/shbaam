#!/usr/bin/env python

import sys
import os.path
import subprocess
import netCDF4
import numpy
import datetime
import fiona
import shapely.geometry
import shapely.prepared
import rtree
import math
import csv
import dateutil.relativedelta

def print_filenames(filenames):
    print('Input filenames: ')
    print(' - swe        netcdf    ' + filenames['swe_netcdf'])
    print(' - polygon    shapefile ' + filenames['polygon_shapefile'])

    print('Output filenames: ')
    print(' - point      shapefile ' + filenames['point_shapefile'])
    print(' - timeseries csv       ' + filenames['timeseries_csv'])
    print(' - map        netcdf    ' + filenames['map_netcdf'])

def check_arguments():
    if len(sys.argv) != 6:
        print('Error: must be 5 arguments')
        raise SystemExit(22)   

def check_input_filenames(input_filenames):
    for input_file in input_filenames:
        try:
            with open(input_file) as f:
                pass
        except IOError as e:
            print('ERROR - Unable to open \'' + input_file + '\'')
            raise SystemExit(22) 

def get_time_step(time_array):
    if len(time_array) > 1:
        return abs(time_array[1] - time_array[0])
    else:
        return 0

def close_swe_netcdf(swe):
    swe['file'].close()

def read_swe_netcdf(swe_netcdf):
    print('Read GLDAS netCDF file')

    f = netCDF4.Dataset(swe_netcdf, 'r')

    swe = {}
    swe['filename']        = swe_netcdf
    swe['file']            = f
    swe['swe']             = f.variables['SWE']
    swe['longitude_size']  = len(f.dimensions['lon'])
    swe['latitude_size']   = len(f.dimensions['lat'])
    swe['time_size']       = len(f.dimensions['time'])
    swe['longitude_array'] = f.variables['lon']
    swe['latitude_array']  = f.variables['lat']
    swe['time_array']      = f.variables['time']
    swe['time_step']       = get_time_step(swe['time_array'])
    swe['fill_value']      = netCDF4.default_fillvals['f4']
    swe['longitude_step']  = abs(swe['longitude_array'][1]                     \
                               - swe['longitude_array'][0])
    swe['latitude_step']   = abs(swe['latitude_array'][1]                      \
                               - swe['latitude_array'][0])


    print(' - The number of longitudes : '         + str(swe['longitude_size']))
    print(' - The number of latitudes : '          + str(swe['latitude_size']))
    print(' - The number of time steps : '         + str(swe['time_size']))
    print(' - The interval size for longitudes : ' + str(swe['longitude_step']))
    print(' - The interval size for latitudes : '  + str(swe['latitude_step']))
    print(' - The interval size for time : '       + str(swe['time_step']))
    print(' - The fill value : '                   + str(swe['fill_value']))

    return swe

def create_point_shapefile(swe, polygon, point_shapefile):
    print('Create a point shapefile with all the GLDAS grid cells')
    if os.path.exists(point_shapefile):
        print(' - Use cached point shapefile.')
        return

    longitude_array = swe['longitude_array']
    latitude_array  = swe['latitude_array']
    polygon_driver  = polygon.driver
    point_driver    = polygon_driver
    polygon_crs     = polygon.crs
    point_crs       = polygon_crs.copy()

    point_schema = {'geometry': 'Point',                                       \
                    'properties': {'lon_index':  'int:4',                      \
                                   'lat_index':  'int:4'}}
    with fiona.open(point_shapefile, 'w', driver = point_driver,
                                          crs    = point_crs,
                                          schema = point_schema) as point: 
        for swe_longitude_index in range(len(longitude_array)):
            longitude = longitude_array[swe_longitude_index]
            for swe_latitude_index in range(len(latitude_array)):
               latitude = latitude_array[swe_latitude_index]
               point_prepared = {'lon_index':  swe_longitude_index,            \
                                 'lat_index':  swe_latitude_index}
               point_geometry = shapely.geometry.mapping(                      \
                                shapely.geometry.Point((longitude,latitude)))
               point.write({'properties': point_prepared,                      \
                            'geometry'  : point_geometry})

    print(' - Point shapefile created')

def create_spatial_index(point):
    print('Create spatial index for the bounds of each point feature')

    index = rtree.index.Index()
    for point_feature in point:
        point_feature_id = int(point_feature['id'])
        # The first argument of index.insert has to be int, not long or str
        point_shape = shapely.geometry.shape(point_feature['geometry'])
        index.insert(point_feature_id, point_shape.bounds)
        # Creates an index between the feature ID and the bounds of that feature

    return index

def find_intersection(polygon, point, index, swe):
    print('Find GLDAS grid cells that intersect with polygon')
    target_lon_indexes = []
    target_lat_indexes = []
  
    for polygon_feature in polygon:
        polygon_shape = shapely.geometry.shape(polygon_feature['geometry'])
        polygon_prepared = shapely.prepared.prep(polygon_shape)
        # A 'prepared' geometry allows for faster processing after
        for point_feature_id in [int(x) for x in list(                         \
                                 index.intersection(polygon_shape.bounds))]:
            point_feature = point[point_feature_id]
            point_shape   = shapely.geometry.shape(point_feature['geometry'])
            if polygon_prepared.contains(point_shape):
                longitude_index = point_feature['properties']['lon_index']
                latitude_index  = point_feature['properties']['lat_index']
                target_lon_indexes.append(longitude_index)
                target_lat_indexes.append(latitude_index)

    print(' - The number of grid cells found : ' + str(len(target_lon_indexes)))
    return {'lon': target_lon_indexes, 'lat': target_lat_indexes}

def compute_target_region(swe, polygon_shapefile, point_shapefile):
    polygon = fiona.open(polygon_shapefile, 'r')
    print(' - The number of polygon features : ' + str(len(polygon)))    

    create_point_shapefile(swe, polygon, point_shapefile) 
    point  = fiona.open(point_shapefile, 'r')
    index  = create_spatial_index(point)
    target = find_intersection(polygon, point, index, swe)

    polygon.close()
    point.close()

    return target

# For each intersecting grid cell, we obtain
# sum of each month's swe change / months
# For a grid cell, long term mean is the average of change.
def compute_long_term_means(target_indexes, swe):
    print('Find long-term mean for each intersecting GLDAS grid cell')
    long_term_means = []
    for longitude_index, latitude_index in                                     \
        zip(target_indexes['lon'], target_indexes['lat']):

        mean = 0
        for time in range(swe['time_size']):
            mean += swe['swe'][time,                                           \
                               latitude_index,                                 \
                               longitude_index] 
        long_term_means.append(mean)
    long_term_means = [x/swe['time_size'] for x in long_term_means]

    return long_term_means

def compute_surface_areas(target_indexes, swe):
    print('Compute surface area of each grid cell')
    surface_areas = [0] * len(target_indexes['lon'])
    
    # in meters
    earth_radius = 6371000 
    for i in range(len(target_indexes['lon'])):
        latitude_index   = target_indexes['lat'][i]
        latitude         = swe['latitude_array'][latitude_index]
        surface_areas[i] = earth_radius * math.radians(swe['latitude_step'])   \
                                        * earth_radius                         \
                                        * math.radians(swe['longitude_step'])  \
                                        * math.cos(math.radians(latitude))

    return surface_areas

def compute_total_surface_area(surface_areas):
    print('Compute total surface area')

    total_surface_area = sum(surface_areas)

    print(' - the area (m2) : ' + str(total_surface_area))
    return total_surface_area

def compute_timeseries(swe, target_indexes, surface_areas,                     \
                       long_term_means, total_surface_area):
    print('Compute total snow water equivalent anomaly timeseries')
    timeseries=[]

    # The original units for snow water equivalent is millimeter, so we need to
    # convert it to meter
    millimeters_in_one_meter = 1000
    for time in range(swe['time_size']):
        cumulative_value = 0
        for longitude_index, latitude_index, surface_area, long_term_mean in   \
            zip(target_indexes['lon'], target_indexes['lat'],                  \
                surface_areas, long_term_means):
            value = (swe['swe'][time, latitude_index, longitude_index]         \
                   - long_term_mean) / millimeters_in_one_meter                \
                    * surface_area
            cumulative_value += value
        timeseries.append(millimeters_in_one_meter *                           \
                          cumulative_value / total_surface_area)

    return timeseries

def compute_time_strings(start_time_string, time_size):
    print(' - Determine time strings')
    time_strings = []
    start_time = datetime.datetime.strptime(start_time_string,                 \
                                            '%Y-%m-%d %H:%M:%S')
    for month_passed in range(time_size):
        time_delta = dateutil.relativedelta.relativedelta(months=month_passed)
        time = start_time + time_delta
        # Only preserve month, day, year
        time_string = time.strftime('%m/%d/%Y')

        time_strings.append(time_string)

    return time_strings

def write_timeseries_csv(timeseries_csv, timeseries, time_size):
    print('Write timeseries_csv')

    time_strings = compute_time_strings('2002-04-01 00:00:00', time_size)

    with open(timeseries_csv, 'wb') as csvfile:
        csvwriter = csv.writer(csvfile, dialect='excel')
        for i in range(len(timeseries)):
            line = [time_strings[i], timeseries[i]] 
            csvwriter.writerow(line)

    print(' - timeseries_csv written')

def write_map_netcdf(map_netcdf, swe, target_indexes, long_term_means):
    print('Write map_netcdf')
    
    print(' - Create netCDF file')
    h = netCDF4.Dataset(map_netcdf, 'w', format='NETCDF3_CLASSIC')

    print(' - Create dimension')
    time_dimension = h.createDimension('time', None)
    lat_dimension  = h.createDimension('lat',  swe['latitude_size'])
    lon_dimension  = h.createDimension('lon',  swe['longitude_size'])
    
    print(' - Create variable')
    time_variable = h.createVariable('time', 'i4', ('time',))
    lat_variable  = h.createVariable('lat',  'f4', ('lat',))
    lon_variable  = h.createVariable('lon',  'f4', ('lon',))
    swe_variable  = h.createVariable('SWE',  'f4', ('time', 'lat', 'lon',),    \
                                     fill_value = swe['fill_value'])
    
    print(' - Populdate global attributes')
    current_time = datetime.datetime.utcnow().replace(microsecond=0)
    version = subprocess.Popen('bash ../version.sh',                           \
                               stdout=subprocess.PIPE,shell=True).communicate()
    version = version[0].rstrip()

    h.Conventions='CF-1.6'
    h.title = ''
    h.institution = ''
    h.source = 'SHBAAM: ' + version + ', GLDAS: '                              \
                + os.path.basename(swe['filename'])
    h.history = 'date created: ' + current_time.isoformat() + '+00:00'
    h.references = 'https://github.com/c-h-david/shbaam/'
    h.comment = ''
    h.featureType = 'timeSeries'

    print(' - Copy existing variable attributes')
    source_time_variable        = swe['file'].variables['time']
    time_variable.standard_name = source_time_variable.standard_name
    time_variable.units         = source_time_variable.units
    time_variable.calendar      = source_time_variable.calendar
   
    source_lat_variable         = swe['file'].variables['lat']
    lat_variable.standard_name  = source_lat_variable.standard_name
    lat_variable.units          = source_lat_variable.units
    lat_variable.axis           = source_lat_variable.axis

    source_lon_variable         = swe['file'].variables['lon']
    lon_variable.standard_name  = source_lon_variable.standard_name
    lon_variable.units          = source_lon_variable.units
    lon_variable.axis           = source_lon_variable.axis

    source_swe_variable         = swe['swe']
    swe_variable.long_name      = source_swe_variable.long_name
    swe_variable.units          = u'mm'
    swe_variable.code           = source_swe_variable.code
    swe_variable.table          = source_swe_variable.table
    swe_variable.missing_value  = source_swe_variable.missing_value
    

    print(' - Populate static data')
    lat_variable[:]  = swe['latitude_array'][:]
    lon_variable[:]  = swe['longitude_array'][:]
    time_variable[:] = swe['time_array'][:]

    print(' - Populate dynamic data')
    for lon_index, lat_index, long_term_mean in                                \
        zip(target_indexes['lon'], target_indexes['lat'], long_term_means):
        for time_index in range(swe['time_size']):
            swe_variable[time_index, lat_index, lon_index] =                   \
            swe['swe'][time_index, lat_index, lon_index] - long_term_mean     

    h.close()
    print(' - Finish writing map_netcdf file')

def print_computations(timeseries):
    print('Check some computations')
    print(' - Average of time series: '+ str(numpy.average(timeseries)))
    print(' - Maximum of time series: '+ str(numpy.max(timeseries)))
    print(' - Minimum of time series: '+ str(numpy.min(timeseries)))

def read_filenames():
    check_arguments()

    filenames = {}
    filenames['swe_netcdf']          = sys.argv[1]
    filenames['polygon_shapefile']   = sys.argv[2]
    filenames['point_shapefile']     = sys.argv[3]
    filenames['timeseries_csv']      = sys.argv[4]
    filenames['map_netcdf']          = sys.argv[5] 

    check_input_filenames([filenames['swe_netcdf'],                            \
                           filenames['polygon_shapefile']])
    print_filenames(filenames)

    return filenames

def main():
    filenames          = read_filenames()

    # Read swe from GLDAS file
    swe                = read_swe_netcdf(filenames['swe_netcdf'])

    target_indexes     = compute_target_region(swe,                            \
                                               filenames['polygon_shapefile'], \
                                               filenames['point_shapefile'])    

    long_term_means    = compute_long_term_means(target_indexes, swe)
    surface_areas      = compute_surface_areas(target_indexes, swe)
    total_surface_area = compute_total_surface_area(surface_areas)
    timeseries         = compute_timeseries(swe, target_indexes, surface_areas,\
                                            long_term_means, total_surface_area)

    write_timeseries_csv(filenames['timeseries_csv'], timeseries,              \
                         swe['time_size'])
    write_map_netcdf(filenames['map_netcdf'], swe, target_indexes,             \
                     long_term_means)

    print_computations(timeseries)

    close_swe_netcdf(swe)    

if __name__ == '__main__':
    main()
