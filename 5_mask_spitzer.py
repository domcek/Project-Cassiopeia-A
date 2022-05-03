"""
The script uses utilises DAOStarFinder to identify point sources in the IR image 
The 'Flux' variable is used as a proxy for the necessary point source mask size
The output is a mask fits file that masks the location of point sources.
"""

from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from photutils import DAOStarFinder
import matplotlib.pyplot as plt
from photutils import CircularAperture
from matplotlib.colors import LogNorm
import numpy as np
import time
import module_plotting as mp
from astropy.wcs import WCS
from matplotlib import cm



outpath = '../data/intermediate_products/'

hdu = fits.open(outpath + '3_Spitzer_bin3_ch1_xheader_jy.fits')
data = hdu[0].data
mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
print((mean, median, std))



star_id = 560
daofind = DAOStarFinder(fwhm=3.0, threshold=3. * std)
sources = daofind(data - median)
print(sources[star_id])

big_sources = np.copy(sources)
big_sources = big_sources[big_sources['flux'] > 100]

positions = (big_sources['xcentroid'], big_sources['ycentroid'])
apertures = CircularAperture(positions, r=3.)

positions_starid = (sources[star_id]['xcentroid'], sources[star_id]['ycentroid'])
apertures_starid = CircularAperture(positions_starid, r=3.)

positions_500 = (sources[500]['xcentroid'], sources[500]['ycentroid'])
apertures_500 = CircularAperture(positions_500, r=3.)

positions_530 = (sources[530]['xcentroid'], sources[530]['ycentroid'])
apertures_530 = CircularAperture(positions_530, r=3.)


plt.clf()
wcs_hdr = WCS(hdu[0].header)
fig, ax = mp.plot_casa(figsize=[8, 8], coords=True, wcs=wcs_hdr)
ax.set_xlim(0, 500)
ax.set_ylim(0, 500)
cmap = cm.get_cmap('binary_r')
cmap.set_bad(color='black')
im = ax.imshow(data, cmap='Greys', origin='lower', norm=LogNorm(vmin=1e-5, vmax=2e-4))
# cbar = mp.set_colorbar(fig, im, title=r'Flux density [Jy]')
apertures.plot(color='white', lw=2.5, alpha=1)
apertures_starid.plot(color='red', lw=2.5, alpha=0.5)
apertures_530.plot(color='red', lw=2.5, alpha=0.5)
apertures_500.plot(color='red', lw=2.5, alpha=0.5)
# plt.savefig(outpath + '5_Spitzer_bin3_ch1_masked_stars_daosf.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)


# picked stars (for std 5.)
print(sources[500])  # small
print(sources[560])  # medium
print(sources[530])  # large

daosf_mask = np.copy(data)
for k in range(len(sources)):
    x = int(np.round(sources[k]['xcentroid']))
    y = int(np.round(sources[k]['ycentroid']))
    circ_size = 2.5
    if sources[k]['flux'] >= 70:
        circ_size = 4.5
    if sources[k]['flux'] < 70 and sources[k]['sharpness'] > 10:
        circ_size = 3.5
    if sources[k]['flux'] <= 10 and sources[k]['sharpness'] > 3:
        circ_size = 2.5
    if sources[k]['flux'] <= 3:
        circ_size = 2
    if sources[k]['flux'] >= 500:
        circ_size = 12

    for i in range(-10, 11):
        for j in range(-10, 11):
            # print(y+i,x+j)
            if np.sqrt(i**2 + j**2) <= circ_size and ((y + i < len(data) and x + j < len(data[0]))):
                daosf_mask[y + i][x + j] = 0  # yes it is reversed correctly
daosf_mask[daosf_mask != 0] = 1

# Special masking for strong pixels
for x in range(234,240): # undetected star
    for y in range(184,189):
        daosf_mask[x,y] = 0
# Strong pixels around masked stars
daosf_mask[245, 334] = 0
daosf_mask[213, 313] = 0
daosf_mask[201, 295] = 0

daosf_mask[daosf_mask == 0] = np.nan

plt.clf()
wcs_hdr = WCS(hdu[0].header)
fig, ax = mp.plot_casa(figsize=[8, 6], coords=True, wcs=wcs_hdr)
cmap = cm.get_cmap('binary_r')
cmap.set_bad(color='black')
im = ax.imshow(data * daosf_mask, cmap='Greys', origin='lower', norm=LogNorm(vmin=1e-5, vmax=2e-4))
cbar = mp.set_colorbar(fig, im, title=r'Flux density [Jy]')
ax.set_title(r'Spitzer masked data')
plt.savefig(outpath + '5_mask_bin3_spitzer_stars_daosf_masked.pdf', bbox_inches='tight')
plt.show()
plt.close(fig)


plt.clf()
fig, ax = plt.subplots(figsize=(6,4))
ax.hist(sources['flux'], bins=500)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('Distribution of the point source fluxes')
plt.savefig(outpath + '5_daosf_point_source_flux_distribution.pdf', bbox_inches='tight')
plt.show()


header = hdu[0].header
header['history'] = '--------------------------------------------------'
header['history'] = 'Used script:    5_mask_spitzer.py'
header['history'] = 'Author:         Vladimir Domcek'
header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
header['history'] = 'Main changes: '
header['history'] = 'Mask stars with gaussians using DAOStarsFinder'

hdu[0].data = daosf_mask
hdu.writeto(outpath + '5_mask_bin3_spitzer_stars_daosf.fits', overwrite=True)