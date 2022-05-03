"""
The script loads radio and radio-to-infrared spectral index maps in lower spatial resolution (10") 
and calculates the difference map that is saved as output. 
"""

import numpy as np
from astropy.io import fits
import time
from astropy.wcs import WCS


# simple power law
def pl_index(flux1, flux2, vnu1, vnu2):
    return (np.log10(flux2) - np.log10(flux1)) / (np.log10(vnu2) - np.log10(vnu1))


def write_into_fits(inputdata, comment, outputfile):
    wcs_proj_file = fits.open(maskpath + '1_chandra_bin3_4.2-6.2_flux_jy.fits')

    new_secondary_hdu = fits.PrimaryHDU()
    new_secondary_hdu.header = wcs_proj_file[0].header
    new_secondary_hdu.data = inputdata

    # create new HDU list
    new_hdul = fits.HDUList()
    new_hdul.append(new_secondary_hdu)
    new_hdul[0].verify('fix')

    header = new_hdul[0].header
    header['history'] = '--------------------------------------------------'
    header['history'] = 'Used script:    monte_carlo_diff.py'
    header['history'] = 'Author:         Vladimir Domcek'
    header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
    header['history'] = comment

    print(outputfile)
    new_hdul.writeto(outputfile, overwrite=True)
    return 0

maskpath = '../data/intermediate_products/'
alphapath = '../data/final_products/'
outpath = '../data/final_products/alpha_diff/'

num_sim = 10000

print('Loading masks')
hdu_radio_hires_lowflux_mask = fits.open('../data/intermediate_products/5_mask_bin3_radio_hires_lowflux.fits')
radio_hires_lowflux_mask = hdu_radio_hires_lowflux_mask[0].data

hdu_spitzer_stars_mask = fits.open('../data/intermediate_products/5_mask_bin3_spitzer_stars_daosf.fits')
spitzer_stars_mask = hdu_spitzer_stars_mask[0].data

hdu_spitzer_stars_mask_smooth = fits.open('../data/intermediate_products/5_mask_bin3_spitzer_stars_smooth.fits')
spitzer_stars_mask_smooth = hdu_spitzer_stars_mask_smooth[0].data

hdu_extinction_mask = fits.open('../data/intermediate_products/5_mask_bin3_nH.fits')
extinction_mask = hdu_extinction_mask[0].data

hdu_frac_mask = fits.open('../data/intermediate_products/5_mask_bin3_radio_thermal_fraction.fits')
frac_mask = hdu_frac_mask[0].data

hdu_BL10mask = fits.open('../data/intermediate_products/5_mask_bin3_alpha_R.fits')
BL10mask = hdu_BL10mask[0].data

the_mask = radio_hires_lowflux_mask * extinction_mask * BL10mask * spitzer_stars_mask_smooth


print('Loading low-res spectral index data')
alpha_RIR = fits.open(alphapath + 'alpha_RIR_low-res/alpha_R-IR_median_lowres_masked.fits')[0].data * the_mask
alpha_RIR_err = fits.open(alphapath + 'alpha_RIR_low-res/alpha_R-IR_err_lowres_masked.fits')[0].data * the_mask

alpha_R = fits.open(alphapath + 'alpha_radio/3_SPIXBL10CORR_bin3_xheader.fits')[0].data * the_mask
alpha_R_err = fits.open(alphapath + 'alpha_radio/3_SPIXERR10_bin3_xheader.fits')[0].data * the_mask

shape = (np.shape(alpha_RIR)[0], np.shape(alpha_RIR)[1], num_sim)
diff_matrix = np.zeros(shape, dtype='float32')

diff_matrix_median, diff_matrix_err = np.zeros_like(alpha_RIR, dtype='float32'), np.zeros_like(alpha_RIR, dtype='float32')

print('Running cycles')
for x in range(len(BL10mask[0]) - 1):
    for y in range(len(BL10mask) - 1):
        if np.isnan(alpha_RIR[y, x]):
            continue
        if np.isnan(alpha_R[y, x]):
            continue

        alpha1 = np.random.normal(alpha_R[y, x], alpha_R_err[y, x], num_sim)
        alpha2 = np.random.normal(alpha_RIR[y, x], alpha_RIR_err[y, x], num_sim)

        diff = alpha2 - alpha1

        #diff_matrix[y][x] = diff
        diff_matrix_median[y][x] = np.nanmedian(diff)
        diff_matrix_err[y][x] = np.nanstd(diff)

print('Saving alpha difference matrix')
np.save(outpath + 'mc_alpha_diff.npy', diff)

print('Saving median maps')
inputdata = diff_matrix_median * the_mask
comment = 'Diff between radio and radio-to-infrared spectral index'
outputfile = outpath + 'diff_alpha_median_masked.fits'
write_into_fits(inputdata, comment, outputfile)

inputdata = diff_matrix_err * the_mask
comment = 'Diff between radio and radio-to-infrared spectral index error values'
outputfile = outpath + 'diff_alpha_err_masked.fits'
write_into_fits(inputdata, comment, outputfile)
