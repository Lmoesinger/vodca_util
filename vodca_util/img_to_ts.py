import numpy as np
import netCDF4 as nc
import pandas as pd
from smecv_grid.grid import SMECV_Grid_v042
import datetime as dt
import os


class Img2Ts(object):
    def __init__(self, in_path, out_fname):
        self.path = in_path
        self.grid = SMECV_Grid_v042()
        self.fillvalue = np.nan
        self.out_fname = out_fname

    def convert_all(self, nchunks):

        gpis = self.grid.activegpis.filled()
        chunks = np.array_split(gpis, nchunks)
        #
        # lon_chunks = np.array_split(self.grid.activearrlon.filled(), nchunks)
        # lat_chunks = np.array_split(self.grid.activearrlat.filled(), nchunks)
        statistics_list = []
        for chunk in chunks:
            statistics = self.convert_gpis(chunk)
            statistics_list.append(statistics)
        variable_dict = {}
        for variable in statistics:
            variable_dict[variable] = pd.concat([stats[variable] for stats in statistics_list], axis=1)

        self.write_to_nc(variable_dict)
        #
        # for cell in self.grid.get_cells().filled():
        #     self.convert_gpis(self.grid.grid_points_for_cell(cell)[0].filled())

    def convert_gpis(self, gpis):
        lonlats = self.grid.gpi2lonlat(gpis)
        lons = lonlats[0]
        lats = lonlats[1]
        cols = ((np.array(lons) - 0.25 / 2) / 0.25 + 720).astype(int)
        rows = (-(np.array(lats) + 0.25 / 2) / 0.25 + 360).astype(int)

        year_dirs = os.listdir(self.path)
        # year_dirs.sort()
        vod_list = []
        time_list = []
        for year_dir in year_dirs:
            print(year_dir)
            year_path = os.path.join(self.path, year_dir)
            fnames = os.listdir(year_path)
            # fnames.sort(key=self._sorting)
            for fname in fnames:
                # print(fname)
                df = nc.Dataset(os.path.join(year_path, fname))
                vod = df['vod'][:].filled(np.nan)
                vod = vod.squeeze()[rows, cols]
                time = df['time'][:].filled(np.nan)[0]
                vod_list.append(vod)
                time_list.append(time)

        df = pd.DataFrame(vod_list, pd.to_datetime('1970-01-01') + pd.to_timedelta(time_list, unit='d'), columns=gpis)
        df.sort_index(inplace=True)
        df_grouped = df.groupby(df.index.year)
        return {'mean': df_grouped.mean(),
                'min': df_grouped.min(),
                'max': df_grouped.max(),
                'n_obs': df_grouped.count(),
                'variance': df_grouped.var()}


    def write_to_nc(self, variable_dict):
        dates = np.array(variable_dict[list(variable_dict)[0]].index)
        lonlats = self.grid.gpi2lonlat(variable_dict[list(variable_dict)[0]].columns)
        lons = lonlats[0]
        lats = lonlats[1]
        cols = ((np.array(lons) - 0.25 / 2) / 0.25 + 720).astype(int)
        rows = (-(np.array(lats) + 0.25 / 2) / 0.25 + 360).astype(int)
        ds = self.create_nc_file(dates)

        for nc_var_name in variable_dict.keys():
            temp = np.empty((len(ds['time']), len(ds['lat']), len(ds['lon'])))
            temp[:] = self.fillvalue
            variable = variable_dict[nc_var_name]
            temp[:, rows, cols] = variable.values
            nc_var = ds.createVariable(nc_var_name, variable.values.dtype, ('time', 'lat', 'lon'), fill_value=self.fillvalue)
            nc_var[:] = temp
        ds.close()

    def create_nc_file(self, dates):
        ds = nc.Dataset(self.out_fname, mode='w')

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
        timevar.units = 'year'
        ds.calendar = 'standard'
        timevar.calendar = 'standard'
        timevar.long_name = 'time'
        timevar[:] = dates
        return ds


    def _sorting(self, L):
        splitup = L.split('_')
        return pd.to_datetime(splitup[-1][:-3])




if __name__ == "__main__":
    in_path = '/data-write/RADAR/vod_merged/released/C-Band/'
    out_fname = '/data-write/RADAR/vod_merged/released/vodca_v01-0_yearly_stats_C-band.nc'
    converter = Img2Ts(in_path, out_fname)
    converter.convert_all(2)