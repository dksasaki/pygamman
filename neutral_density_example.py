
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from gamman import gamman_module as gn

ds = xr.open_dataset('/home/otel/Desktop/deletar/pyroms_tests/202202_pbs_ceresIV_2010-05-01_00_00_00_2010-05-08_00_00_00.nc')

xm, ym = np.meshgrid(ds.longitude.values, ds.latitude.values)

s = ds.so[0,0].values.ravel()
t = ds.thetao[0,0].values.ravel()
p = 0*xm.ravel()
lon = xm.ravel()
lat = ym.ravel()

imask =np.where(~np.isnan(t))

a = gn.gamma_n(s[imask],t[imask],p[imask],lon[imask],lat[imask], [imask[0].size])

gamma = np.zeros_like(s)
gamma[imask] = a[0]
gamma = gamma.reshape(xm.shape)