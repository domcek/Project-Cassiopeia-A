"""
The script reprojects 10" resolution C-band image (provided by Dr. Delaney) 
into a common WCS grid. Chandra observation WCS grid is used as a template
"""

import time
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
from reproject import reproject_exact

outpath = '../data/intermediate_products/'
Chandra_infile = '1_chandra_bin3_4.2-6.2_flux_jy.fits'

hdu1 = fits.open((outpath + Chandra_infile))
hdu2_C10 = fits.open('../data/original/C10.FITS')
# BEAM TO JY
# 17 comes from ratio between beam surface area (7.28'') and 2.5" image pixel surface area (0.43")
# 16 comes as correction from smoothing from 2.5" to 10". Difference of 4x turns into 16 in surface area
beam_to_jy = 17*16 # only works for this specific case

array_C10, footprint_C10 = reproject_exact((hdu2_C10[0].data[0, 0, :, :]/(beam_to_jy), WCS(hdu2_C10[0].header).sub(2)),
                                               hdu1[0].header)

xray_pixel_size_in_arcsec = hdu1[0].header['CDELT2'] * 3600
radio_pixel_size_in_arcsec = hdu2_C10[0].header['CDELT2'] * 3600
ratio = xray_pixel_size_in_arcsec ** 2 / radio_pixel_size_in_arcsec ** 2

print('Ratio:', ratio)
# correcting for the pixel change
print('Pixel size in arcsec', hdu1[0].header['CDELT2'] * 3600)
array_C10 *= ratio

new_secondary_hdu = fits.PrimaryHDU()
# copy spitzer header into new secondary
new_secondary_hdu.header = hdu2_C10[0].header
# copy chandra header projection part into new primary
hdu2_C10[0].header = hdu1[0].header[-80:-54]
hdu2_C10[0].data = array_C10

# create new HDU list
new_hdul = fits.HDUList()
new_hdul.append(hdu2_C10[0])
new_hdul.append(new_secondary_hdu)
new_hdul[0].verify('fix')
# new_hdul[1].verify('fix')

header = new_hdul[0].header
header['history'] = '--------------------------------------------------'
header['history'] = 'Used script:    3_reproject_C10.py'
header['history'] = 'Author:         Vladimir Domcek'
header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
header['history'] = 'Main changes: '
header['history'] = 'Delaney radio data reprojected to chandra coord system and binsize (1.5'')'
header['history'] = 'Previous radio header moved into 2nd extension of header'
header['history'] = 'Reprojected file: radio/C10.FITS'
header['history'] = 'Reproject according to file: ' + Chandra_infile
header['history'] = 'Units Jy/px'

new_hdul.writeto(outpath + '3_C10_bin3_xheader.fits', overwrite=True)

