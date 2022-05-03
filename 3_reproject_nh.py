"""
The script reprojects nH image (provided by Dr. Una Hwang) into a common WCS grid.
Chandra observation WCS grid is used as a template
"""

from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from reproject import reproject_exact
from astropy.io import fits
import time


outpath = '../data/intermediate_products/'
Chandra_infile = '1_chandra_bin3_4.2-6.2_flux_jy.fits'

# REPROJECTING NH
hdu1 = fits.open((outpath + Chandra_infile))
hdu2_nH = fits.open(('../data/original/grid_nh.img')) # Hwang and Laming 2012

array_nH, footprint_nH = reproject_exact(hdu2_nH[0], hdu1[0].header)

# RATIO of pixel does not apply here because I want the same numbers in different grid
# array_H -= 578

new_secondary_hdu = fits.PrimaryHDU()
# copy spitzer header into new secondary
new_secondary_hdu.header = hdu2_nH[0].header
# copy chandra header projection part into new primary
hdu2_nH[0].header = hdu1[0].header[-80:-54]
hdu2_nH[0].data = array_nH

# create new HDU list
new_hdul = fits.HDUList()
new_hdul.append(hdu2_nH[0])
new_hdul.append(new_secondary_hdu)
new_hdul[0].verify('fix')
new_hdul[1].verify('fix')

header = new_hdul[0].header
header['history'] = '--------------------------------------------------'
header['history'] = 'Used script:    3_reproject_nh.py'
header['history'] = 'Author:         Vladimir Domcek'
header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
header['history'] = 'Main changes: '
header['history'] = 'NH map reprojected to chandra coord. system and binsize (1.5'')'
header['history'] = 'Previous NH header moved into 2nd extension of header'
header['history'] = 'Reprojected file: hwang_casa_2012/grid_nh.img'
header['history'] = 'Reproject according to file: ' + Chandra_infile

new_hdul.writeto(outpath + '3_grid_nh_repro.fits', overwrite=True)



