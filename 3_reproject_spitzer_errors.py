"""
The script reprojects 2.5" resolution of Spitzer error map (from Spitzer archive) 
into a common WCS grid. Chandra observation WCS grid is used as a template.
"""

from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from reproject import reproject_exact
from astropy.io import fits
import time
import numpy as np

Spitzer_inpath = '../data/original/'
Spitzer_infile = 'SPITZER_I1_34836224_0000_1_A31354259_munc.fits' # level 2 product of spitzer archive

outpath = '../data/intermediate_products/'
Chandra_infile = '1_chandra_bin3_4.2-6.2_flux_jy.fits'

hdu1 = fits.open((outpath + Chandra_infile))
hdu2_ch1 = fits.open((Spitzer_inpath + Spitzer_infile))
hdu2_ch1[0].data = (hdu2_ch1[0].data / (42545 * 1 / (hdu2_ch1[0].header['PXSCAL2'])**2))**2
# second squared because of the variance
print(hdu2_ch1[0].header['PXSCAL2'])
print((42545 * 1 / (hdu2_ch1[0].header['PXSCAL2'])**2))  # conversion from MJy/sr


array_ch1, footprint_ch1 = reproject_exact(hdu2_ch1[0], hdu1[0].header)
xray_pixel_size_in_arcsec = hdu1[0].header['CDELT2'] * 3600
spitzer_pixel_size_in_arcsec = hdu2_ch1[0].header['PXSCAL2']
ratio = xray_pixel_size_in_arcsec**2 / spitzer_pixel_size_in_arcsec**2
new_array_ch1 = array_ch1 * ratio

# new_array_ch1[np.isnan(new_array_ch1)] = 0  # changing nan values into 0

# Header update

# create second header
new_secondary_hdu = fits.PrimaryHDU()
# copy spitzer header into new secondary
new_secondary_hdu.header = hdu2_ch1[0].header
# copy chandra header projection part into new primary
hdu2_ch1[0].header = hdu1[0].header[-80:-54]
hdu2_ch1[0].data = np.sqrt(new_array_ch1)

# create new HDU list
new_hdul = fits.HDUList()
new_hdul.append(hdu2_ch1[0])
new_hdul.append(new_secondary_hdu)
new_hdul[0].verify('fix')
new_hdul[1].verify('fix')

header = new_hdul[0].header
header['history'] = '--------------------------------------------------'
header['history'] = 'Used script:    3_reproject_spitzer_errors.py'
header['history'] = 'Author:         Vladimir Domcek'
header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
header['history'] = 'Main changes: '
header['history'] = 'Spitzer data errors reprojected to chandra coord. system and binsize (1.5'')'
header['history'] = 'Previous Spitzer header moved into 2nd extension of header'
header['history'] = 'Reprojected Spitzer file: ' + Spitzer_infile
header['history'] = 'Reproject according to file: ' + Chandra_infile

new_hdul.writeto(outpath + '3_Spitzer_bin3_ch1_err_xheader_jy.fits', overwrite=True)

