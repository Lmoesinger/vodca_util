import numpy as np
import netCDF4 as nc
import pandas as pd
from smecv_grid.grid import SMECV_Grid_v042
import os


class Img2Ts(object):
    def __init__(self, path):
        self.path = path
        self.grid = SMECV_Grid_v042()

    def resample(self):


        df = nc.Dataset(self.path)