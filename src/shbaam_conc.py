#!/usr/bin/env python

import sys
import os.path
import netCDF4
import numpy

def print_usage():
    print('Usage: ./shbaam_conc.py [input file]... [ouput file]')

# Check if there are should be at least 2 input files, 1 output file.
# Which is a total of at least 4 arguments.
def check_program_parameters():
    if len(sys.argv) <= 4:
        print_usage()
        print('Error: There should be at least 2 input files, 1 output file')
        raise SystemExit(22)

# All files except for the last one are input files
def get_input_output_files(files):
    return (files[:-1], files[-1])

def print_information(input_files, output_file):
    print('input  files: ')
    for input_file in input_files:
      print('  ' + str(input_file))
    print('output files: ')
    print('  ' + str(output_file))
    print('')

def abort_if_input_files_unreadable(input_files):
    for input_file in input_files:
        try:
            with open(input_file) as f:
                pass
        except IOError as e:
            print('Error: Unable to open ' + input_file)
            print(e)
            raise SystemExit(22)

def abort_if_output_file_exists(output_file):
    if os.path.exists(output_file):
        print('Aborts: output file \'' + output_file + '\' already exists.')
        raise SystemExit(22)


def create_attributes(output_netcdf4, input_netcdf4):
    for attribute_name in input_netcdf4.ncattrs():
        if attribute_name == u'history':
            continue
        output_netcdf4.setncattr(attribute_name, input_netcdf4.getncattr(attribute_name))

def create_dimensions(output_netcdf4, input_netcdf4):
    for dimension_name, dimension in input_netcdf4.dimensions.iteritems():
        size = None if dimension.isunlimited() else len(dimension)
        output_netcdf4.createDimension(dimension_name, size)

def create_variables(output_netcdf4, input_netcdf4, input_files):
    for variable_name, variable in input_netcdf4.variables.iteritems():
        output_variable = output_netcdf4.createVariable(variable_name, variable.datatype, variable.dimensions)

        for attribute_name in variable.ncattrs():
            if variable_name == 'time' and attribute_name == 'units':
                output_variable.setncattr(attribute_name, u'months since 2002-04-01 00:00:00')
            else:
                output_variable.setncattr(attribute_name, variable.getncattr(attribute_name))


    for variable_name, variable in input_netcdf4.variables.iteritems():
        if variable_name == 'time':
            output_netcdf4.variables['time'][:] = numpy.arange(0, len(input_files), 1)
        if variable_name == 'lat':
            output_netcdf4.variables['lat'][:] = input_netcdf4.variables['lat'][:]   
        if variable_name == 'lon':
            output_netcdf4.variables['lon'][:] = input_netcdf4.variables['lon'][:]
        if variable_name == 'SWE':
            all_swe = []
            for i in range(len(input_files)):
                input_file = netCDF4.Dataset(input_files[i], 'r')
                # Here assumes all input files has 1 unique time step
                all_swe.append(input_file.variables['SWE'][0]) 
                input_file.close()
            output_netcdf4.variables['SWE'][:] = all_swe[:]
        if variable_name == 'Canint':
            all_canint = []
            for i in range(len(input_files)):
                input_file = netCDF4.Dataset(input_files[i], 'r')
                # Here assumes all input files has 1 unique time step
                all_canint.append(input_file.variables['Canint'][0])
                input_file.close()
            output_netcdf4.variables['Canint'][:] = all_canint[:]

def handle_concatenation(input_files, output_file):
    output_netcdf4 = netCDF4.Dataset(output_file, 'w', format='NETCDF4')
    input_netcdf4 = netCDF4.Dataset(input_files[0], 'r')

    create_attributes(output_netcdf4, input_netcdf4)
    create_dimensions(output_netcdf4, input_netcdf4)
    create_variables(output_netcdf4, input_netcdf4, input_files)

    output_netcdf4.close()
    input_netcdf4.close()
   


def main():
    check_program_parameters()
    input_files, output_file = get_input_output_files(sys.argv[1:])
    print_information(input_files, output_file)
    abort_if_input_files_unreadable(input_files)    
    abort_if_output_file_exists(output_file)

    handle_concatenation(input_files, output_file)
    
    print('Success')


if __name__ == '__main__':
    main()
