import numpy as np
from astropy.io import fits

def sens_error(file, cam):
    with fits.open(file) as f:
        s = f[0].data
        w = s[0,:]
        s = s[1:,:]
        print(s.shape)
        s[s==0] = np.nan
        avgsens = np.nanpercentile(s, [1,5,16,50,84,95,99], axis=0)
        print(avgsens.shape)
        print(w[0:10])
        np.savetxt(f'sens_percentiles-{cam}.csv', np.hstack((w[:,None],avgsens.T)), delimiter=',')

sens_error('sens-b-v1.1.fits', 'b')
sens_error('sens-r-v1.1.fits', 'r')
sens_error('sens-z-v1.1.fits', 'z')
