import numpy as np
from astropy.time import Time
from astropy.io import fits

with open('lvmFFrame.txt') as file:
    lines = [line.rstrip() for line in file]

lines_b = [l for l in lines if l[-15]=='b']
lines_r = [l for l in lines if l[-15]=='r']
lines_z = [l for l in lines if l[-15]=='z']

for frame in lines_z:
    with fits.open(frame) as f:
        h = f[0].header
        try:
            tile = h['TILE_ID']
            if tile<=1111:
                continue
            cal = h['FLUXCAL']
        except KeyError:
            continue
        if cal!='NONE':
            obstime = Time(h['OBSTIME'], format='isot', scale='utc').mjd
            exposure = h['EXPOSURE']
            std_sens = h.get('STDSENMZ', 0)
            sci_sens = h.get('SCISENMZ', 0)
            print(exposure, obstime, tile, std_sens, sci_sens)
