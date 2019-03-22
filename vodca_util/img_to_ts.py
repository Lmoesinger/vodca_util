import numpy as np
import netCDF4 as nc
import pandas as pd
from smecv_grid.grid import SMECV_Grid_v042
import os
from vodca_util.functions_to_apply import get_monthly_stats, get_yearly_stats


class ApplyFun(object):
    def __init__(self, in_path, out_dict):
        '''
        initializes the object
        :param in_path: string, path to folder created by unzipping the vodca data
        :param out_dict: dictionary, each entry is a pair of "filename: function handle".
        '''

        self.path = in_path
        self.grid = SMECV_Grid_v042()
        self.fillvalue = np.nan
        self.out_dict = out_dict

    def convert_all(self, nchunks):
        '''
        Calculates the yearly statistics globally and saves them to a netcdf file
        :param nchunks: number of chunks. E.g. with nchunks=10 only about 1/10 of the data is in the memory at the same
        time. Try to keep nchunks as low as possible, as each time it has to open all files again.
        :return:
        '''

        gpis = self.grid.activegpis.filled()
        chunks = np.array_split(gpis, nchunks)
        chunk_statistics_list = []
        for chunk in chunks:
            stats_dict = {}
            df = self.read_data(chunk)
            for fun in self.out_dict:
                stats_dict[fun] = out_dict[fun](df)
            chunk_statistics_list.append(stats_dict)
            del df

        for fun in self.out_dict:
            stats = [chunk[fun] for chunk in chunk_statistics_list]
            variables = list(stats[0])
            variable_dict = {}
            for variable in variables:
                variable_dict[variable] = pd.concat([stat[variable] for stat in stats], axis=1)
            self.write_to_nc(variable_dict, fun)

    def read_data(self, gpis):
        lonlats = self.grid.gpi2lonlat(gpis)
        lons = lonlats[0]
        lats = lonlats[1]
        cols = ((np.array(lons) - 0.25 / 2) / 0.25 + 720).astype(int)
        rows = (-(np.array(lats) + 0.25 / 2) / 0.25 + 360).astype(int)

        year_dirs = os.listdir(self.path)
        vod_list = []
        time_list = []
        for year_dir in year_dirs:
            print(year_dir)
            year_path = os.path.join(self.path, year_dir)
            fnames = os.listdir(year_path)
            for fname in fnames:
                df = nc.Dataset(os.path.join(year_path, fname))
                vod = df['vod'][:].filled(np.nan)
                vod = vod.squeeze()[rows, cols]
                time = df['time'][:].filled(np.nan)[0]
                vod_list.append(vod)
                time_list.append(time)

        df = pd.DataFrame(vod_list, pd.to_datetime('1970-01-01') + pd.to_timedelta(time_list, unit='d'), columns=gpis)
        df.sort_index(inplace=True)
        return df

    def write_to_nc(self, variable_dict, fname):
        '''
        Writes the calculated statistics to a netcdf file
        :param variable_dict: dictionary, each entry is a yearly statistic, e.g. mean, variance, etc.
        :param fname: file name
        :return: nothing
        '''
        dates = np.array(variable_dict[list(variable_dict)[0]].index)
        lonlats = self.grid.gpi2lonlat(variable_dict[list(variable_dict)[0]].columns)
        lons = lonlats[0]
        lats = lonlats[1]
        cols = ((np.array(lons) - 0.25 / 2) / 0.25 + 720).astype(int)
        rows = (-(np.array(lats) + 0.25 / 2) / 0.25 + 360).astype(int)
        ds = self.create_nc_file(dates, fname)

        for nc_var_name in variable_dict.keys():
            temp = np.empty((len(ds['time']), len(ds['lat']), len(ds['lon'])))
            temp[:] = self.fillvalue
            variable = variable_dict[nc_var_name]
            temp[:, rows, cols] = variable.values
            nc_var = ds.createVariable(nc_var_name, temp.dtype, ('time', 'lat', 'lon'), fill_value=self.fillvalue)
            nc_var[:] = temp
        ds.close()

    def create_nc_file(self, dates, fname):
        '''
        Sets up the output netcdf file. Creates the file, creates and sets the values the dimensions
        :param dates: np.array, containing the years where observations are available
        :param fname: string, file name
        :return: ds, the initialized netcdf4 object
        '''
        ds = nc.Dataset(fname, mode='w')

        lon_1d = np.unique(self.grid.arrlon.filled())
        londim = ds.createDimension('lon', lon_1d.__len__())
        lonvar = ds.createVariable('lon', lon_1d.dtype, ('lon',))
        lonvar.units = 'degrees_east'
        lonvar.standard_name = 'longitude'
        lonvar.valid_range = [-180.0, 180.0]
        lonvar[:] = lon_1d

        lat_1d = np.unique(self.grid.arrlat.filled())[::-1]
        latdim = ds.createDimension('lat', lat_1d.__len__())
        latvar = ds.createVariable('lat', lat_1d.dtype, ('lat',))
        latvar.units = 'degrees_north'
        latvar.standard_name = 'latitude'
        latvar.valid_range = [-90.0, 90.0]
        latvar[:] = lat_1d

        timedim = ds.createDimension("time", len(dates))
        timevar = ds.createVariable('time', dates.dtype, ('time',))
        # timevar.units = 'year'
        ds.calendar = 'standard'
        timevar.calendar = 'standard'
        timevar.long_name = 'time'
        timevar[:] = dates
        return ds


if __name__ == "__main__":
    """
    Adjust values as needed!
    """
    # path to folder containint the vodca data.
    # Assumes that no additional files have been added to it after unzipping it
    in_path = '/data-write/RADAR/vod_merged/released/X-band/'
    # path to folder where output files will be created
    out_path = '/data-write/RADAR/vod_merged/released/'

    # dictionary, each name is the filename while the entry is the function to be applied
    out_dict = {os.path.join(out_path, 'vodca_v01-0_yearly_stats_X-band.nc'): get_yearly_stats,
                os.path.join(out_path, 'vodca_v01-0_monthly_stats_X-band.nc'): get_monthly_stats}
    # Number of chunks.
    #   low -> increased processing speed but high memory requirement
    #   high -> lower speed but also lower memory required
    # Therefore you want to take the lowest number your memory can support.
    nchunks = 2

    # initialize converter
    converter = ApplyFun(in_path, out_dict)
    # run it
    converter.convert_all(nchunks)

