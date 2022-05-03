"""
The script reprojects 10" resolution of Radio spectral index map 
(provided by Dr. Delaney) into a common WCS grid. 
Chandra observation WCS grid is used as a template.
"""

import time
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
from reproject import reproject_exact

outpath = '../data/final_products/alpha_radio/'
Chandra_infile = '1_chandra_bin3_4.2-6.2_flux_jy.fits'

hdu1 = fits.open('../data/intermediate_products/' + Chandra_infile)
hdu2_SPIXBL10 = fits.open('../data/intermediate_products/SPIXBL10CORR.FITS')

array_SPIXBL10, footprint_SPIXBL10 = reproject_exact((hdu2_SPIXBL10[0].data[0, 0, :, :], WCS(hdu2_SPIXBL10[0].header).sub(2)),
                                               hdu1[0].header)

new_secondary_hdu = fits.PrimaryHDU()
# copy spitzer header into new secondary
new_secondary_hdu.header = hdu2_SPIXBL10[0].header
# copy chandra header projection part into new primary
hdu2_SPIXBL10[0].header = hdu1[0].header[-80:-54]
hdu2_SPIXBL10[0].data = array_SPIXBL10

# create new HDU list
new_hdul = fits.HDUList()
new_hdul.append(hdu2_SPIXBL10[0])
new_hdul.append(new_secondary_hdu)
new_hdul[0].verify('fix')
# new_hdul[1].verify('fix')

header = new_hdul[0].header
header['history'] = '--------------------------------------------------'
header['history'] = 'Used script:    3_reproject_SPIXBL10cor.py'
header['history'] = 'Author:         Vladimir Domcek'
header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
header['history'] = 'Main changes: '
header['history'] = 'Delaney radio data reprojected to chandra coord system and binsize (1.5'')'
header['history'] = 'Previous radio header moved into 2nd extension of header'
header['history'] = 'Reprojected file: radio/SPIXBL10.FITS'
header['history'] = 'Reproject according to file: ' + Chandra_infile
header['history'] = 'Units Jy/beam'

new_hdul.writeto(outpath + '3_SPIXBL10CORR_bin3_xheader.fits', overwrite=True)

