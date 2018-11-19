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

def abort_if_input_files_unreadable(input_files):
    for input_file in input_files:
        try:
            with open(input_file) as f:
                pass
        except IOError as e:
            print('Error: Unable to open ' + input_file)
            print(e)
            raise SystemExit(22)

def copy_global_attributes(output_netcdf4, input_netcdf4):
    print('Copy global attributes')
    for attribute_name in input_netcdf4.ncattrs():
        if attribute_name == u'history':
            continue
        output_netcdf4.setncattr(attribute_name, input_netcdf4.getncattr(attribute_name))

def copy_dimensions(output_netcdf4, input_netcdf4):
    print('Copy dimensions')
    for dimension_name, dimension in input_netcdf4.dimensions.iteritems():
        size = None if dimension.isunlimited() else len(dimension)
        output_netcdf4.createDimension(dimension_name, size)

def copy_variables(output_netcdf4, input_netcdf4, input_files):
    print('Copy variables')
    for variable_name, variable in input_netcdf4.variables.iteritems():
        output_variable = output_netcdf4.createVariable(variable_name, variable.datatype, variable.dimensions)

        for attribute_name in variable.ncattrs():
           output_variable.setncattr(attribute_name, variable.getncattr(attribute_name))

    output_netcdf4.variables['lat'][:]  = input_netcdf4.variables['lat'][:]   
    output_netcdf4.variables['lon'][:]  = input_netcdf4.variables['lon'][:]

    for i in range(len(input_files)):
        input_file = netCDF4.Dataset(input_files[i], 'r')

        output_netcdf4.variables['time'][i] = 0 #numpy.arange(0, len(input_files), 1)
        output_netcdf4.variables['SWE'][i]    = input_file.variables['SWE'][0]
        output_netcdf4.variables['Canint'][i] = input_file.variables['Canint'][0]

        input_file.close()

def handle_concatenation(input_files, output_file):
    output_netcdf4 = netCDF4.Dataset(output_file, 'w', format='NETCDF4')
    input_netcdf4 = netCDF4.Dataset(input_files[0], 'r')

    copy_global_attributes(output_netcdf4, input_netcdf4)
    copy_dimensions(output_netcdf4, input_netcdf4)
    copy_variables(output_netcdf4, input_netcdf4, input_files)

    output_netcdf4.close()
    input_netcdf4.close()

def main():
    check_program_parameters()
    input_files, output_file = get_input_output_files(sys.argv[1:])
    print_information(input_files, output_file)
    abort_if_input_files_unreadable(input_files)    

    handle_concatenation(input_files, output_file)
    
    print('Files concatenated sucessfully')


if __name__ == '__main__':
    main()
