"""
The script uses reprojected nH map to filter out gap regions in the map
and save it as a mask fits file
"""
import numpy as np
import module_plotting as mp
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from matplotlib import cm
import time

outpath = '../data/intermediate_products/'

hdu = fits.open(outpath + '3_grid_nh_repro.fits')
img = hdu[0].data

img_mask = np.copy(img)

treshold = 1.0 # treshold assumed from the nH image
img_mask[img_mask < treshold] = np.nan
img_mask[img_mask >= treshold] = 1


header = hdu[0].header
header['history'] = '--------------------------------------------------'
header['history'] = 'Used script:    5_mask_nH.py'
header['history'] = 'Author:         Vladimir Domcek'
header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
header['history'] = 'Main changes: '
header['history'] = 'Mask bellow treshold: ' + str(treshold) 

hdu[0].data = img_mask
hdu.writeto(outpath + '5_mask_bin3_nH.fits', overwrite=True)
