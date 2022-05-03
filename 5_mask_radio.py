"""
The script uses reprojected radio map to filter out regions with less reliable flux
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

hdu_radio_4_72Ghz = fits.open(outpath + '3_Hires_bin3_xheader_fluxed_2009_nomask.fits')
img_radio_4_72Ghz = hdu_radio_4_72Ghz[0].data

radio_lowflux_mask = np.copy(img_radio_4_72Ghz)

treshold = 5e-3  # treshold assumed from the radio image
radio_lowflux_mask[radio_lowflux_mask < treshold] = np.nan
radio_lowflux_mask[radio_lowflux_mask >= treshold] = 1


plt.clf()
wcs_hdr = WCS(hdu_radio_4_72Ghz[0].header)
fig, ax = mp.plot_casa(figsize=[8, 6], coords=True, wcs=wcs_hdr)
cmap = cm.get_cmap('binary_r')
cmap.set_bad(color='black')
im = ax.imshow(radio_lowflux_mask, origin='lower',cmap=cmap, vmin=0, vmax=1)  # , norm=LogNorm(vmin=2e-4, vmax=7e-3))
# cbar = mp.set_colorbar(fig, im, title=r'mask')
ax.set_title(r'Radio lowflux mask')
plt.savefig(outpath + '5_mask_bin3_radio_hires_lowflux.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

header = hdu_radio_4_72Ghz[0].header
header['history'] = '--------------------------------------------------'
header['history'] = 'Used script:    5_mask_radio.py'
header['history'] = 'Author:         Vladimir Domcek'
header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
header['history'] = 'Main changes: '
header['history'] = 'Mask bellow treshold: ' + str(treshold) + 'Jy/px'

hdu_radio_4_72Ghz[0].data =radio_lowflux_mask
hdu_radio_4_72Ghz.writeto(outpath + '5_mask_bin3_radio_hires_lowflux.fits', overwrite=True)
