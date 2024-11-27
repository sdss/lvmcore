#! /usr/bin/env python

from astropy.io import fits
import os.path as path

def process(line, outfile):
    dir = path.dirname(line)
    fl = path.basename(line)
    exposure = fl.split('-')[-1].split('.')[0]
    exposure = int(exposure)
    with fits.open(line) as hdu:
        h = hdu[0].header
        exptime = h['EXPTIME']
        flux = h['FLUXCAL']
        tileid = h['TILE_ID']
        mjd = h['MJD']
    #print(f'{dir},{mjd},{exposure},{tileid},{exptime},{flux}\n')
    outfile.write(f'{dir},{mjd},{exposure},{tileid},{exptime},{flux}\n')

i=0
with open('cframes.txt') as f:
    with open('drpall_txt.csv', 'w') as outf:
        for line in f:
            process(line[:-1], outf)
            if i%100 == 0:
                print(i)
            i = i+1




