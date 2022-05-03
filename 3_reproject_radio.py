"""
The script reprojects 2.5" resolution High resolution radio image (provided by Dr. Delaney) 
into a common WCS grid. The output of 2_radio_calib.py is used for flux calibration of the image.
Chandra observation WCS grid is used as a template
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
from reproject import reproject_exact
import time

outpath = '../data/intermediate_products/'
Chandra_infile = '1_chandra_bin3_4.2-6.2_flux_jy.fits'

hdu1 = fits.open((outpath + Chandra_infile))
hdu2_hires = fits.open('../data/original/HIRES.RADIO.FITS')

array_hires, footprint_hires = reproject_exact((hdu2_hires[0].data[0, 0, :, :], WCS(hdu2_hires[0].header).sub(2)), hdu1[0].header)

xray_pixel_size_in_arcsec = hdu1[0].header['CDELT2'] * 3600
radio_pixel_size_in_arcsec = hdu2_hires[0].header['CDELT2'] * 3600
ratio = xray_pixel_size_in_arcsec**2 / radio_pixel_size_in_arcsec**2

print('Ratio:', ratio)
# correcting for the pixel change
print('Pixel size in arcsec', hdu1[0].header['CDELT2'] * 3600)
array_hires *= ratio


new_secondary_hdu = fits.PrimaryHDU()
# copy spitzer header into new secondary
new_secondary_hdu.header = hdu2_hires[0].header
# copy chandra header projection part into new primary
hdu2_hires[0].header = hdu1[0].header[-80:-54]
hdu2_hires[0].data = array_hires

# create new HDU list
new_hdul = fits.HDUList()
new_hdul.append(hdu2_hires[0])
new_hdul.append(new_secondary_hdu)
new_hdul[0].verify('fix')
# new_hdul[1].verify('fix')

header = new_hdul[0].header
header['history'] = '--------------------------------------------------'
header['history'] = 'Used script:    3_reproject_radio.py'
header['history'] = 'Author:         Vladimir Domcek'
header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
header['history'] = 'Main changes: '
header['history'] = 'Delaney radio data reprojected to chandra coord system and binsize (1.5'')'
header['history'] = 'Previous radio header moved into 2nd extension of header'
header['history'] = 'Reprojected file: radio/HIRES.RADIO.FITS'
header['history'] = 'Reproject according to file: ' + Chandra_infile
header['history'] = 'Not fluxed'

new_hdul.writeto(outpath + '3_Hires_bin3_xheader.fits', overwrite=True)

# FLUX CALIBRATION
# array_hires[array_hires <= 5e-2] = np.nan  # cutting out noise to 0

# new_hdul[0].data = array_hires
# new_hdul.writeto(outpath + '3_Hires_bin3_xheader.fits', overwrite=True)

sum_image = 11895.7  # measured in ds9
# green_cat = (2720*4.72**(-0.77))
perly = 697.6  # result of 2_radio_calibration
array_hires *= perly / sum_image

new_hdul[0].data = array_hires
header = new_hdul[0].header[:-2]
header['history'] = 'Flux set to match Perly 2014 measurements'

new_hdul.writeto(outpath + '3_Hires_bin3_xheader_fluxed_2009_nomask.fits', overwrite=True)
