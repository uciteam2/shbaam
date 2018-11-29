#!/usr/bin/env python

import sys
import netCDF4
import unittest

sys.path.append('../src')


from shbaam_conc import *

class Test(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        folder = '../input/GLDAS/GLDAS_VIC10_M/2002/'

        self.input_files = [folder + 'GLDAS_VIC10_M.A200204.001.grb.SUB.nc4', folder + 'GLDAS_VIC10_M.A200205.001.grb.SUB.nc4', \
                            folder + 'GLDAS_VIC10_M.A200206.001.grb.SUB.nc4', folder + 'GLDAS_VIC10_M.A200207.001.grb.SUB.nc4', \
                            folder + 'GLDAS_VIC10_M.A200208.001.grb.SUB.nc4', folder + 'GLDAS_VIC10_M.A200209.001.grb.SUB.nc4', \
                            folder + 'GLDAS_VIC10_M.A200210.001.grb.SUB.nc4', folder + 'GLDAS_VIC10_M.A200211.001.grb.SUB.nc4', \
                            folder + 'GLDAS_VIC10_M.A200212.001.grb.SUB.nc4']
        self.output_file =  folder + 'GLDAS_VIC10_M.A200204_200212.nc4'

        handle_concatenation(self.input_files, self.output_file)

        self.output_netcdf4 = netCDF4.Dataset(self.output_file, 'r')
        self.input_netcdf4 = netCDF4.Dataset(self.input_files[0], 'r')


    @classmethod
    def tearDownClass(self):
        self.output_netcdf4.close()
        self.input_netcdf4.close()

    # Test cases below
    
    def test_get_input_output_files(self):
        files = ['input1', 'input2', 'input3', 'output']
        input_files, output_file = get_input_output_files(files)

        self.assertEqual(input_files, ['input1', 'input2', 'input3'])
        self.assertEqual(output_file, 'output')

    def test_create_attributes(self):
        for attribute_name in self.input_netcdf4.ncattrs():
            self.assertEqual(self.output_netcdf4.getncattr(attribute_name), self.input_netcdf4.getncattr(attribute_name))

    def test_create_variables(self):
        for variable_name in self.input_netcdf4.variables:
            # Check variables' attributes
            for attribute_name in self.input_netcdf4.variables[variable_name].ncattrs():
                # TODO(Jiayou Shi): need to check this, but it is now empty
                if (variable_name == 'time' and attribute_name == 'units'):
                    self.assertEquals(self.output_netcdf4.variables[variable_name].getncattr(attribute_name), u'months since 2002-04-01 00:00:00')
                else:
                    self.assertTrue(numpy.array_equal(self.output_netcdf4.variables[variable_name].getncattr(attribute_name), \
                                                      self.input_netcdf4.variables[variable_name].getncattr(attribute_name)))
        
            # Check variables' values
            if variable_name == 'time':
                self.assertTrue(numpy.array_equal(self.output_netcdf4.variables['time'], \
                                                  numpy.arange(0, len(self.input_files), 1)))
            if variable_name == 'lat':
                self.assertTrue(numpy.array_equal(self.output_netcdf4.variables['lat'], \
                                                  self.input_netcdf4.variables['lat']))
            if variable_name == 'lon':
                self.assertTrue(numpy.array_equal(self.output_netcdf4.variables['lon'], \
                                                  self.input_netcdf4.variables['lon'])) 
            # Loop through output's SWE, it should be an array of each file' SWE from first to the last
            # And it should match the order from 04 to 12
            if variable_name == 'SWE':
                self.assertEqual(len(self.output_netcdf4.variables['SWE']), len(self.input_files))
                for time in range(len(self.input_files)):
                    for lat in range(len(self.output_netcdf4.variables['lat'])):
                        input_file = netCDF4.Dataset(self.input_files[time], 'r')
                        self.assertTrue(numpy.array_equal(self.output_netcdf4.variables['SWE'][time][lat], \
                                                          input_file.variables['SWE'][0][lat]))
                        input_file.close()


            if variable_name == 'Canint':
                self.assertEqual(len(self.output_netcdf4.variables['Canint']), len(self.input_files))
                for time in range(len(self.input_files)):
                    for lat in range(len(self.output_netcdf4.variables['Canint'])):
                        input_file = netCDF4.Dataset(self.input_files[time], 'r')
                        self.assertTrue(numpy.array_equal(self.output_netcdf4.variables['Canint'][time][lat], \
                                                          input_file.variables['Canint'][0][lat]))
                        input_file.close()

    def test_create_dimensions(self):
        for dimension_name in self.input_netcdf4.dimensions:
            if (dimension_name == 'time'):
                continue

            self.assertEqual(self.input_netcdf4.dimensions[dimension_name].name, \
                             self.output_netcdf4.dimensions[dimension_name].name)
            self.assertEqual(self.input_netcdf4.dimensions[dimension_name].size, \
                             self.output_netcdf4.dimensions[dimension_name].size)
            self.assertEqual(self.input_netcdf4.dimensions[dimension_name].isunlimited(), \
                             self.output_netcdf4.dimensions[dimension_name].isunlimited())

if __name__ == '__main__':
    unittest.main()
