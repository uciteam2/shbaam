#!/usr/bin/env python

import sys
import getopt
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


advanced               = False
filename_tags_provided = False
add_all_variable_names = False

def get_time_step(time_array):
    """Calculates the difference of the the time indexes \
        Returns the int time_step"""
    if len(time_array) > 1:
        return abs(time_array[1] - time_array[0])
    else:
        return 0


def close_netcdf(data):
    """Closes the input netCDF file"""
    data['file'].close()


def get_fill_value(f):
    """Takes the netCDF file \
        Returns the int fill_value"""
    fill_value = netCDF4.default_fillvals['f4']
    if 'RUNSF' in f.variables:
        variable = f.variables['RUNSF']
        if '_FillValue' in variable.ncattrs():
            fill_value = variable._FillValue
        else:
            fill_value = None

    return fill_value


def read_netcdf(netcdf):
    """opens the input netCDF file                          \
        Returns the Dict, containing the necessary values for \
            calculations and later file creation"""
    print('Read GLDAS netCDF file')

    f = netCDF4.Dataset(netcdf, 'r')

    data = {}

    data['data'] = []
    data['filename'] = netcdf
    data['file'] = f
    data['longitude_size'] = len(f.dimensions['lon'])
    data['latitude_size'] = len(f.dimensions['lat'])
    data['time_size'] = len(f.dimensions['time'])
    data['longitude_array'] = f.variables['lon']
    data['latitude_array'] = f.variables['lat']
    data['time_array'] = f.variables['time']
    data['time_step'] = get_time_step(data['time_array'])
    data['fill_value'] = get_fill_value(f)
    data['longitude_step'] = abs(data['longitude_array'][1] \
                                 - data['longitude_array'][0])
    data['latitude_step'] = abs(data['latitude_array'][1] \
                                - data['latitude_array'][0])

    print(' - Number of longitudes : ' + str(data['longitude_size']))
    print(' - Number of latitudes : ' + str(data['latitude_size']))
    print(' - Number of time steps : ' + str(data['time_size']))
    print(' - Interval size for longitudes : ' + str(data['longitude_step']))
    print(' - Interval size for latitudes : ' + str(data['latitude_step']))
    print(' - Interval size for time : ' + str(data['time_step']))
    print(' - Fill value : ' + str(data['fill_value']))

    return data


def create_point_shapefile(data, polygon, point_shapefile):
    """Creates a output point shapefile, containing the geometry \
        and properties of the Point                             \
        If the file already exists, the cached one will be used instead"""
    print('Create a point shapefile with all the GLDAS grid cells')
    if os.path.exists(point_shapefile):
        print(' - Use cached point shapefile.')
        return

    longitude_array = data['longitude_array']
    latitude_array = data['latitude_array']
    polygon_driver = polygon.driver
    point_driver = polygon_driver
    polygon_crs = polygon.crs
    point_crs = polygon_crs.copy()

    point_schema = {'geometry': 'Point', \
                    'properties': {'lon_index': 'int:4', \
                                   'lat_index': 'int:4'}}
    with fiona.open(point_shapefile, 'w', driver=point_driver,
                    crs=point_crs,
                    schema=point_schema) as point:
        for data_longitude_index in range(len(longitude_array)):
            longitude = longitude_array[data_longitude_index]
            if longitude > 180:
                longitude -= 360
            for data_latitude_index in range(len(latitude_array)):
                latitude = latitude_array[data_latitude_index]
                point_prepared = {'lon_index': data_longitude_index, \
                                  'lat_index': data_latitude_index}
                point_geometry = shapely.geometry.mapping( \
                    shapely.geometry.Point((longitude, latitude)))
                point.write({'properties': point_prepared, \
                             'geometry': point_geometry})

    print(' - Point shapefile created')


def create_spatial_index(point):
    """Return an Index object containing               \
        the int id and 4-tuple bound of the point shape"""
    print('Create spatial index for the bounds of each point feature')

    index = rtree.index.Index()
    for point_feature in point:
        point_feature_id = int(point_feature['id'])
        # The first argument of index.insert has to be int, not long or str
        point_shape = shapely.geometry.shape(point_feature['geometry'])
        index.insert(point_feature_id, point_shape.bounds)
        # Creates an index between the feature ID and the bounds of that feature

    return index


def find_intersection(polygon, point, index, data):
    """Iterates through the features of the polygon                             \
            and the index object's list of intersection bounds                  \
        Returns the Dict, containing the list of longitude and latitude indexes \
            that intersect with the polygon"""
    print('Find GLDAS grid cells that intersect with polygon')
    target_lon_indexes = []
    target_lat_indexes = []

    for polygon_feature in polygon:
        polygon_shape = shapely.geometry.shape(polygon_feature['geometry'])
        polygon_prepared = shapely.prepared.prep(polygon_shape)
        # A 'prepared' geometry allows for faster processing after
        for point_feature_id in [int(x) for x in list( \
                index.intersection(polygon_shape.bounds))]:
            point_feature = point[point_feature_id]
            point_shape = shapely.geometry.shape(point_feature['geometry'])
            if polygon_prepared.contains(point_shape):
                longitude_index = point_feature['properties']['lon_index']
                latitude_index = point_feature['properties']['lat_index']
                target_lon_indexes.append(longitude_index)
                target_lat_indexes.append(latitude_index)

    print(' - The number of grid cells found : ' + str(len(target_lon_indexes)))
    return {'lon': target_lon_indexes, 'lat': target_lat_indexes}


def compute_target_region(data, polygon_shapefile, point_shapefile):
    """Opens the input netCDF input file,                                         \
            used to create the output point shapefile                             \
        Returns the Dict, containing the list of longitude and latitude indexes   \
            that intersect with the polygon, computed by find_intersction"""
    polygon = fiona.open(polygon_shapefile, 'r')
    print(' - The number of polygon features : ' + str(len(polygon)))

    create_point_shapefile(data, polygon, point_shapefile)
    point = fiona.open(point_shapefile, 'r')
    index = create_spatial_index(point)
    target = find_intersection(polygon, point, index, data)

    polygon.close()
    point.close()

    return target


def compute_long_term_means(target_indexes, data):
    """Returns the list of means by iterating through the data, \
    provided by the input netCDF file"""
    print('Find long-term mean for each intersecting GLDAS grid cell')
    long_term_means = []
    for longitude_index, latitude_index in \
            zip(target_indexes['lon'], target_indexes['lat']):

        mean = 0
        for time in range(data['time_size']):
            mean += data['cur_data'][time, \
                                     latitude_index, \
                                     longitude_index]
        long_term_means.append(mean)
    long_term_means = [x / data['time_size'] for x in long_term_means]

    return long_term_means


def compute_surface_areas(target_indexes, data):
    """Returns the list of surface areas,                       \
    computed by using the data, provided by the input netCDF file"""

    print('Compute surface area of each grid cell')
    surface_areas = [0] * len(target_indexes['lon'])

    # in meters
    earth_radius = 6371000
    for i in range(len(target_indexes['lon'])):
        latitude_index = target_indexes['lat'][i]
        latitude = data['latitude_array'][latitude_index]
        surface_areas[i] = earth_radius * math.radians(data['latitude_step']) \
                           * earth_radius \
                           * math.radians(data['longitude_step']) \
                           * math.cos(math.radians(latitude))

    return surface_areas


def compute_total_surface_area(surface_areas):
    """Returns the int total Surface area by        \
        adding the list of given surface areas"""
    print('Compute total surface area')

    total_surface_area = sum(surface_areas)

    print(' - the area (m2) : ' + str(total_surface_area))
    return total_surface_area


def compute_timeseries(data, target_indexes, surface_areas, \
                       long_term_means, total_surface_area):
    """Returns the list of timeseries,                      \
        by dividing cumulative value and the total surface area"""
    print('Compute total snow water equivalent anomaly timeseries')
    timeseries = []

    # The original units for snow water equivalent is millimeter, so we need to
    # convert it to meter
    millimeters_in_one_meter = 1000
    for time in range(data['time_size']):
        cumulative_value = 0
        for longitude_index, latitude_index, surface_area, long_term_mean in \
                zip(target_indexes['lon'], target_indexes['lat'], \
                    surface_areas, long_term_means):
            value = (data['cur_data'][time, latitude_index, longitude_index] \
                     - long_term_mean) / millimeters_in_one_meter \
                    * surface_area
            cumulative_value += value
        timeseries.append(millimeters_in_one_meter * \
                          cumulative_value / total_surface_area)

    return timeseries


def compute_time_strings(start_time_string, time_size):
    """Returns the list of time strings, determined by months"""
    print(' - Determine time strings')
    time_strings = []
    start_time = datetime.datetime.strptime(start_time_string, \
                                            '%Y-%m-%d %H:%M:%S')
    for month_passed in range(time_size):
        time_delta = dateutil.relativedelta.relativedelta(months=month_passed)
        time = start_time + time_delta
        # Only preserve month, day, year
        time_string = time.strftime('%m/%d/%Y')

        time_strings.append(time_string)

    return time_strings


def write_timeseries_csv(timeseries_csv, timeseries, time_size, start_time):
    """Takes the list of time strings and           \
    creates the output CSV file of the time strings"""
    print('Write timeseries_csv: ' + timeseries_csv)

    time_strings = compute_time_strings(start_time, time_size)

    with open(timeseries_csv, 'wb') as csvfile:
        csvwriter = csv.writer(csvfile, dialect='excel')
        for i in range(len(timeseries)):
            line = [time_strings[i], timeseries[i]]
            csvwriter.writerow(line)

    print(' - timeseries_csv written')


def write_map_netcdf(map_netcdf, data, target_indexes, command_info):
    """Creates the ouput nc file, containing dimensions, variables, \
    global attributes, variable attributes, and static data"""
    print('Write map_netcdf: ' + map_netcdf)

    print(' - Create netCDF file')
    h = netCDF4.Dataset(map_netcdf, 'w', format='NETCDF3_CLASSIC')

    print(' - Create dimension')
    time_dimension = h.createDimension('time', None)
    lat_dimension = h.createDimension('lat', data['latitude_size'])
    lon_dimension = h.createDimension('lon', data['longitude_size'])

    print(' - Create variable')
    time_variable = h.createVariable('time', 'f8', ('time',))
    lat_variable = h.createVariable('lat', 'f8', ('lat',))
    lon_variable = h.createVariable('lon', 'f8', ('lon',))

    print(' - Populdate global attributes')
    current_time = datetime.datetime.utcnow().replace(microsecond=0)
    version = subprocess.Popen('bash ../version.sh', \
                               stdout=subprocess.PIPE, shell=True).communicate()
    version = version[0].rstrip()

    h.Conventions = 'CF-1.6'
    h.title = ''
    h.institution = ''
    h.source = 'SHBAAM: ' + version + ', ' + command_info['source'] + ': ' \
               + os.path.basename(data['filename'])
    h.history = 'date created: ' + current_time.isoformat() + '+00:00'
    h.references = 'https://github.com/c-h-david/shbaam/'
    h.comment = ''
    h.featureType = 'timeSeries'

    print(' - Copy existing variable attributes')
    source_time_variable = data['file'].variables['time']
    time_variable.standard_name = source_time_variable.standard_name
    time_variable.units = source_time_variable.units
    time_variable.calendar = source_time_variable.calendar

    source_lat_variable = data['file'].variables['lat']
    lat_variable.standard_name = source_lat_variable.standard_name
    lat_variable.units = source_lat_variable.units
    lat_variable.axis = source_lat_variable.axis

    source_lon_variable = data['file'].variables['lon']
    lon_variable.standard_name = source_lon_variable.standard_name
    lon_variable.units = source_lon_variable.units
    lon_variable.axis = source_lon_variable.axis

    print(' - Populate static data')
    lat_variable[:] = data['latitude_array'][:]
    lon_variable[:] = data['longitude_array'][:]
    time_variable[:] = data['time_array'][:]

    for data_info in data['data']:
        variable_name = data_info[0]
        values = data_info[1]
        long_term_means = data_info[2]

        data_variable = h.createVariable(variable_name, \
                                         'f4', ('time', 'lat', 'lon',), \
                                         fill_value=data['fill_value'])

        source_data_variable = values
        data_variable.long_name = source_data_variable.long_name
        data_variable.units = u'mm'
        data_variable.code = source_data_variable.code
        data_variable.table = source_data_variable.table
        data_variable.missing_value = source_data_variable.missing_value

    print(' - Populate dynamic data')
    for lon_index, lat_index, long_term_mean in \
            zip(target_indexes['lon'], target_indexes['lat'], long_term_means):
        for time_index in range(data['time_size']):
            data_variable[time_index, lat_index, lon_index] = \
                values[time_index, lat_index, lon_index] - long_term_mean

    h.close()
    print(' - Finish writing map_netcdf file')


def print_computations(timeseries):
    """Prints the time series results"""
    print('Check some computations')
    print(' - Average of time series: ' + str(numpy.average(timeseries)))
    print(' - Maximum of time series: ' + str(numpy.max(timeseries)))
    print(' - Minimum of time series: ' + str(numpy.min(timeseries)))


def form_filename(command_info, file_type, variable_name=None):
    """Returns the parsed filename and its file extensions"""
    if file_type == 'shp':
        filename = '.'.join(filter(None, [command_info['source'], \
                                          command_info['model'], 'pnt_tst.shp']))
    elif file_type == 'csv':
        filename = '_'.join(filter(None, \
                                   ['timeseries', variable_name, \
                                    command_info['location'], 'tst.csv']))
    elif file_type == 'nc':
        filename = '_'.join(filter(None, ['map_swea', \
                                          command_info['location'], 'tst.nc']))

    return command_info['output_folder'] + '/' + filename


def read_command_info(filenames, command_info):
    """Checks if the file exists and can print the usage of the file    \
    based on whether the command is in advanced mode"""
    if advanced and len(filenames) != 3:
        print_usage()
        print('ERROR: Exactly 5 pathnames must be provided')
        print('       You must provide a netcdf file, ' + \
              'a polygon shapefile, and a output folder')
        raise SystemExit(22)
    if not advanced and len(filenames) != 5:
        print_usage()
        print('ERROR: Exactly 5 pathnames must be provided.')
        print('       You must provide a netcdf file, a polygon shapefile,' + \
              'pathes for a point shapefile, a timeseries file, a map netcdf file')
        raise SystemExit(22)

    command_info['netcdf'] = filenames[0]
    command_info['polygon_shapefile'] = filenames[1]

    if advanced:
        command_info['output_folder'] = filenames[2].strip('/')
        command_info['point_shapefile'] = form_filename(command_info, 'shp')
        command_info['map_netcdf'] = form_filename(command_info, 'nc')
    else:
        command_info['point_shapefile'] = filenames[2]
        command_info['timeseries_csv'] = filenames[3]
        command_info['map_netcdf'] = filenames[4]

    check_if_input_files_exist([command_info['netcdf'], \
                                command_info['polygon_shapefile']])

    return command_info


def print_usage():
    """Prints in the terminal of the instructions of the commands"""
    print('Usage: ')
    print('  ./shbaam_sewa.py [input_netcdf_file] [input_polygon_shapefile]' + \
          '[output_temporary_shapefile] [output_timeseries_csv] [output_map_netcdf]')
    print('  ./shbaam_swea.py -d [data1],[data2] ' + \
          '[netcdf_file] [shapefile] [output_folder]')
    print('')

    print('Options:')
    print('  -a,  Select all the variable names in the netcdf file to compute.')
    print('  -d,  Select the variable names to compute, SWE, Canint, etc. ' + \
          'Seperate the names with a comma. The default variable name is SWE.')
    print('  -p,  Add tags to the output filenames. Provide the tags in ' + \
          'the following order, [DATA SOURCE],[DATA MODEL],[TARGET LOCATION]')
    print('  -t,  Provide start time string in the format of YYYY-MM-DD ' + \
          'HH:MM:SS' + '. The default is 2002-04-01 00:00:00')

    print('')
    print('Examples:')
    print('  ./shbaam_swea.py GLDAS_netcdf.nc Nepal.shp ' + \
          'GLDAS.VIC.shp timeseries_swea_Nepal.csv map_swea_Nepal.nc')
    print('  ./shbaam_swea.py -d SWE,Canint -p GLDAS,VIC,Nepal' + \
          ' -t \'2002-04-01 00:00:00\' netcdf_file.nc4 Nepal.shp ./output_folder')


def read_command_line():
    """Checks command line and parses the options                   \
        Returns the dict, containing the Source, Model, Location,   \
        start_time, and variable_names"""
    global advanced
    global filename_tags_provided
    global add_all_variable_names

    try:
        options, arguments = getopt.getopt(sys.argv[1:], 'hd:p:t:a')
    except getopt.GetoptError:
        print_usage()

    command_info = {'source': '', 'model': '', 'location': '', \
                    'start_time': '2002-04-01 00:00:00', 'variable_names': []}
    for option, argument in options:
        if option == '-a':
            advanced = True
            add_all_variable_names = True
        elif option == '-h':
            print_usage()
            raise SystemExit(0)
        elif option == '-d':
            advanced = True
            command_info['variable_names'] = argument.split(',')
        elif option == '-p':
            filename_tags_provided = True
            tags = argument.split(',')
            if len(tags) != 3:
                print_usage()
                print('ERROR: You must provide the -p option in the order ' + \
                      'of SOURCE,MODEL,LOCATION')
                raise SystemExit(22)
            command_info['source'] = tags[0]
            command_info['model'] = tags[1]
            command_info['location'] = tags[2]
        elif option == '-t':
            command_info['start_time'] = argument

    # Set SWE as the default data to be processed
    if not advanced:
        command_info['variable_names'] = ['SWE']

    read_command_info(arguments, command_info)

    return command_info


def set_up(data, variable_name, command_info):
    """Registers each variable before computing     \
    curr_data is the currently processed variable"""
    if advanced:
        command_info['timeseries_csv'] = form_filename(command_info, \
                                                       'csv', variable_name)
    data['cur_data'] = data['file'].variables[variable_name]


def print_command_info(command_info):
    """Prints the information of the name and type of  input/output files"""
    print('Input files: ')
    print(' - data       netcdf    ' + command_info['netcdf'])
    print(' - polygon    shapefile ' + command_info['polygon_shapefile'])

    if advanced:
        print('Output folder:')
        print(' - ' + command_info['output_folder'])
    else:
        print('Output files:')
        print(' - point      shapefile ' + command_info['point_shapefile'])
        print(' - timeseries csv       ' + command_info['timeseries_csv'])
        print(' - map        netcdf    ' + command_info['map_netcdf'])


def check_if_input_files_exist(input_command_info):
    """Checks if the input files exist"""
    for input_file in input_command_info:
        try:
            with open(input_file) as f:
                pass
        except IOError as e:
            print('ERROR - Unable to open \'' + input_file + '\'')
            raise SystemExit(22)


def get_computable_variables(command_info, data):
    """Returns the list of proper variable names"""
    return [variable_name for variable_name, \
                              variable in data['file'].variables.iteritems() \
            if variable_name not in ['time', 'lat', 'lon']]


def validate_variable_name(command_info, data):
    """Validates if the variable_name used is        \
        in accordance computable/proper variable names"""
    all_variables = get_computable_variables(command_info, data)
    if advanced:
        for variable_name in command_info['variable_names']:
            if variable_name not in all_variables:
                print('ERROR: variable name ' + str(variable_name) + \
                      ' is not in the netcdf file')
                raise SystemExit(22)


def post_command_info_setup(command_info, data):
    """Ensures the initial variable names in command_info contains nothing  \
    and proceeds to add the list of proper variable names from the function,\
    get_computable_variables"""
    if add_all_variable_names:
        assert len(command_info['variable_names']) == 0
        command_info['variable_names'] = get_computable_variables(command_info, data)


def main():
    command_info = read_command_line()
    print_command_info(command_info)

    data = read_netcdf(command_info['netcdf'])
    post_command_info_setup(command_info, data)
    validate_variable_name(command_info, data)
    target_indexes = compute_target_region(data, \
                                           command_info['polygon_shapefile'], \
                                           command_info['point_shapefile'])
    surface_areas = compute_surface_areas(target_indexes, data)
    total_surface_area = compute_total_surface_area(surface_areas)
    for variable_name in command_info['variable_names']:
        print('------------------ ' + variable_name + ' ------------------')
        set_up(data, variable_name, command_info)

        long_term_means = compute_long_term_means(target_indexes, data)
        timeseries = compute_timeseries(data, target_indexes, surface_areas, \
                                        long_term_means, total_surface_area)

        # Write to output
        write_timeseries_csv(command_info['timeseries_csv'], timeseries, \
                             data['time_size'], command_info['start_time'])

        print_computations(timeseries)
        data['data'].append((variable_name, data['cur_data'], long_term_means))

    print('------------------ End of Computation ------------------')
    write_map_netcdf(command_info['map_netcdf'], data, target_indexes, command_info)

    close_netcdf(data)


if __name__ == '__main__':
    main()
