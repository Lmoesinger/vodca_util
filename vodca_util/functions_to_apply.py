import numpy as np
import pandas as pd


def get_yearly_stats(df):
    '''
    Given a list of gpis, read the vod data of them and calculate the yearly statistics for them.
    :param gpis: list-like, grid point ids
    :return: statistics of each gpi
    '''
    df_grouped = df.groupby(df.index.year)
    n_obs = df_grouped.count()
    n_obs[np.isnan(n_obs)] = 0
    return {'mean': df_grouped.mean(),
            'min': df_grouped.min(),
            'max': df_grouped.max(),
            'n_obs': n_obs,
            'variance': df_grouped.var()}


def get_monthly_stats(df):
    df_grouped = df.groupby([df.index.year, df.index.month])
    monthly_mean = df_grouped.mean()
    monthly_mean.index = [int(str(date[0]) + str(date[1]).zfill(2)) for date in monthly_mean.index]
    # monthly_mean.set_index(monthly_mean.index.map(lambda t: pd.datetime(*t, 1)), inplace=True)
    return {'mean': monthly_mean}
