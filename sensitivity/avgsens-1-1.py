#! /usr/bin/env python

from astropy.io import fits
from astropy.table import Table
import numpy as np
from scipy.integrate import simpson
import os.path as path
import os
import sys
import getopt

LVM_ROOT='/uufs/chpc.utah.edu/common/home/sdss50/'
#DRPALL = os.getenv('LVM_SPECTRO_REDUX')+'/1.1.0/drpall-1.1.0.fits'
DRPALL = 'drpall_txt.csv'
t = Table.read(DRPALL)

select = ((t['mjd']>60255) * (t['exptime']==900) * (t['tileid']>11111) * (t['fluxcal']=='STD'))
t = t[select]
t = t[2:-1]

print(len(t))

def avg_sens(p, res, line):
    if path.exists(p) is True:
        with fits.open(p) as frame:
            h = frame['PRIMARY'].header
            if h['FLUXCAL'] == 'STD':
                if line==0:
                    res[0,:] = frame['WAVE'].data
                    print(res[0,0:10])
                    line += 1
                s = Table(frame['FLUXCAL_STD'].data)
                try:
                    res[line,:] = s['mean']
                except Exception:
                    print('WARNING no sens info: '+p)
            else:
                print('WARNING fluxcal==F in '+p)
    else:
        print('WARNING no file:'+p)


def mk_sens():
    res_b = np.zeros((len(t)+1, 4401))
    res_r = np.zeros((len(t)+1, 3591))
    res_z = np.zeros((len(t)+1, 4561))
    line = 0
    for obs in t:
        d = obs['location']
        expnum = obs['expnum']
        avg_sens(d+'/'+f'lvmFFrame-b-{expnum:08d}.fits', res_b, line)
        avg_sens(d+'/'+f'lvmFFrame-r-{expnum:08d}.fits', res_r, line)
        avg_sens(d+'/'+f'lvmFFrame-z-{expnum:08d}.fits', res_z, line)
        line += 1
        if(line%100 == 0):
            print(line)
    
    fits.writeto('sens-b-v1.1.fits', res_b, overwrite=True)
    fits.writeto('sens-r-v1.1.fits', res_r, overwrite=True)
    fits.writeto('sens-z-v1.1.fits', res_z, overwrite=True)


def mean_sens(cam, file, filter=False):
    with fits.open(file) as f:
        s = f[0].data
    w = s[0,:]
    s = s[1:,:]
    #return w, s
    n,_ = s.shape
    if cam == "b":
        filt = np.exp(-0.5 * ((w - 4500) / 250) ** 2)
    elif cam == "r":
        filt = np.exp(-0.5 * ((w - 6500) / 250) ** 2)
    else:
        filt = np.exp(-0.5 * ((w - 8500) / 250) ** 2)
    I2 = simpson(y=filt, x=w)
    for i in range(n):
        I1 = simpson(y=s[i,:] * filt, x=w)
        if np.isfinite(I1) and I1>0:
            s[i,:] /= (I1/I2)
        else:
            s[i,:] = np.nan
    # filter out bad calibrations, origin still unclear
    if filter is True:
        r = np.nanstd(s[:,1500:2500], axis=1)
        m = np.nanmedian(r)
        s[np.where(r>1.1*m),:] = np.nan
    return w, s

def mk_mean_sens():
    w, s = mean_sens('b', 'sens-b-v1.1.fits', filter=True)
    m = np.nanmedian(s, axis=0)
    np.savetxt("mean-sens-b-v1.1.csv", np.c_[w,m], fmt='%.8g', delimiter=",")
    w, s = mean_sens('r', 'sens-r-v1.1.fits', filter=True)
    m = np.nanmedian(s, axis=0)
    np.savetxt("mean-sens-r-v1.1.csv", np.c_[w,m], fmt='%.8g', delimiter=",")
    w, s = mean_sens('z', 'sens-z-v1.1.fits', filter=True)
    m = np.nanmedian(s, axis=0)
    np.savetxt("mean-sens-z-v1.1.csv", np.c_[w,m], fmt='%.8g', delimiter=",")

    
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hsm", ["help", "output="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    for o, a in opts:
        if o in ("-m", "--mean-sens-only"):
            print("Extracting mean sensitivity curves from FITS files ...")
            mk_mean_sens()
            sys.exit()
        elif o in ("-h", "--help"):
            sys.exit()
    print("Creating sensitivity FITS files ...")
    mk_sens()
    print("Extracting mean sensitivity curves from FITS files ...")
    mk_mean_sens()
                
if __name__ == "__main__":
    main()
