"""
The script loads the flux calibrated 2.5" resolution radio data and creates mask 
for regions where the flux is more then 10% of the average rms obtained from
four corners of the radio image. The output is a mask fits file.
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


rms = 7e-4 # measured in 4 corners of Hires image and averaging out

radio_fraction_map =  rms / img_radio_4_72Ghz 

header = hdu_radio_4_72Ghz[0].header
header['history'] = '--------------------------------------------------'
header['history'] = 'Used script:    4_5_mask_radio_fraction_map.py'
header['history'] = 'Author:         Vladimir Domcek'
header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
header['history'] = 'Main changes: '
header['history'] = 'Fraction map (rms in 4 corners / radio flux) '
header['history'] = 'RMS used = ' + str(rms)

hdu_radio_4_72Ghz[0].data = radio_fraction_map
hdu_radio_4_72Ghz.writeto(outpath + '4_Hires_fraction_map.fits', overwrite=True)



radio_fraction_mask = np.copy(radio_fraction_map)

treshold = 0.1  # treshold 10%
radio_fraction_mask[radio_fraction_mask > treshold] = np.nan
radio_fraction_mask[radio_fraction_mask < 0] = np.nan
radio_fraction_mask[radio_fraction_mask <= treshold] = 1


plt.clf()
wcs_hdr = WCS(hdu_radio_4_72Ghz[0].header)
fig, ax = mp.plot_casa(figsize=[8, 6], coords=True, wcs=wcs_hdr)
cmap = cm.get_cmap('binary_r')
cmap.set_bad(color='black')
im = ax.imshow(radio_fraction_mask, origin='lower',cmap=cmap, vmin=0, vmax=1)  # , norm=LogNorm(vmin=2e-4, vmax=7e-3))
# cbar = mp.set_colorbar(fig, im, title=r'mask')
ax.set_title(r'Radio thermal fraction mask')
plt.savefig(outpath + '5_mask_bin3_radio_thermal_fraction.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)

header['history'] = 'Mask everything above 10% of rms'

hdu_radio_4_72Ghz[0].data = radio_fraction_mask
hdu_radio_4_72Ghz.writeto(outpath + '5_mask_bin3_radio_thermal_fraction.fits', overwrite=True)
