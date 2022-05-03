"""
The script loads radio, infrared and extinction data in higher spatial resoulution (2.5") 
and using Monte Carlo method goes through calibration procedure as described 
in Section 2 (and Fig.2) of the accompanied paper. 
The outputs are the calibrated radio, infrared images 
and the radio-to-infrared spectral index map
"""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import time


# N_H to Av conversion of Guver & Ozel (2009)
def nh_av(nh):
    av = nh / (np.random.normal(2.21, 0.09, nh.size) * 1e21)
    return av


# extinction curve equation of Indebetouw et al. (2004)
def extinc_curve(lam0, sim):
    return np.random.normal(0.61, 0.04, sim) - np.random.normal(2.22, 0.17, sim) * np.log10(lam0) + \
        np.random.normal(1.21, 0.23, sim) * (np.log10(lam0))**2


# extinction correction
def ext_cor(flux, extinction):
    F0 = 277.5
    m = -2.5 * np.log10(flux / F0)
    m_c = m - extinction
    return F0 * 10**(-0.4 * m_c)


# simple power law
def pl_index(flux1, flux2, vnu1, vnu2):
    return (np.log10(flux2) - np.log10(flux1)) / (np.log10(vnu2) - np.log10(vnu1))


def write_into_fits(inputdata, comment, outputfile):
    # wcs_proj_file = fits.open(inpath + '1_chandra_projection_header.fits')
    wcs_proj_file = fits.open(inpath + '1_chandra_bin3_4.2-6.2_flux_jy.fits')

    new_secondary_hdu = fits.PrimaryHDU()
    new_secondary_hdu.header = wcs_proj_file[0].header
    new_secondary_hdu.data = inputdata

    # create new HDU list
    new_hdul = fits.HDUList()
    new_hdul.append(new_secondary_hdu)
    new_hdul[0].verify('fix')

    header = new_hdul[0].header
    header['history'] = '--------------------------------------------------'
    header['history'] = 'Used script:    monte_carlo_high-res.py'
    header['history'] = 'Author:         Vladimir Domcek'
    header['history'] = 'Date and time:  ' + time.strftime('%X %x %Z')
    header['history'] = comment

    print(outputfile)
    new_hdul.writeto(outputfile, overwrite=True)
    return 0


inpath = '../data/intermediate_products/'
outpath = '../data/final_products/alpha_RIR_high-res/'
num_sim = 10000
vnu1 = 4.72e9
vnu2 = (3.6 * u.micron).to(u.Hz, equivalencies=u.spectral()).value
nh_frac_error = 0.1
radio_frac_error = 0.1
bkg_val, bkg_err = 1e-5, 1.5e-6
# bkg_val, bkg_err = 1e-5, 3e-6 # wide bkg test


# MASKS
print('Loading masks')
hdu_radio_hires_lowflux_mask = fits.open(inpath + '5_mask_bin3_radio_hires_lowflux.fits')
radio_hires_lowflux_mask = hdu_radio_hires_lowflux_mask[0].data

hdu_spitzer_stars_mask = fits.open(inpath + '5_mask_bin3_spitzer_stars_daosf.fits')
spitzer_stars_mask = hdu_spitzer_stars_mask[0].data

hdu_extinction_mask = fits.open(inpath + '5_mask_bin3_nH.fits')
extinction_mask = hdu_extinction_mask[0].data

the_mask = radio_hires_lowflux_mask * spitzer_stars_mask * extinction_mask

# LOAD DATA
print('Loading data')
nh_map = fits.open(inpath + '3_grid_nh_repro.fits')[0].data
Spitzer_flux_map = fits.open(inpath + '3_Spitzer_bin3_ch1_xheader_jy.fits')[0].data * the_mask
Spitzer_flux_err = fits.open(inpath + '3_Spitzer_bin3_ch1_err_xheader_jy.fits')[0].data * the_mask
radio_map = fits.open(inpath + '3_Hires_bin3_xheader_fluxed_2009_nomask.fits')[0].data * the_mask

# LOAD CONSTANT DISTRIBUTIONS
print('Simulating constant value distributions')
A36_Ak_extc = extinc_curve(3.6, num_sim)  # 3.6 micron
A36_Av_extc = 10 ** A36_Ak_extc * (1. / 8.8)  # Scaling to optical

Spitzer_bkg_distribution = np.random.normal(bkg_val, bkg_err, num_sim)

shape = (np.shape(radio_map)[0], np.shape(radio_map)[1], num_sim)

alpha_matrix = np.zeros(shape, dtype='float32')
radio_flux_matrix = np.zeros(shape, dtype='float32')
spitzer_flux_matrix = np.zeros(shape, dtype='float32')

spitzer_flux_median, spitzer_flux_err = np.zeros_like(radio_map, dtype='float32'), np.zeros_like(radio_map, dtype='float32')
radio_flux_median, radio_flux_err = np.zeros_like(radio_map, dtype='float32'), np.zeros_like(radio_map, dtype='float32')
alpha_matrix_median, alpha_matrix_err = np.zeros_like(radio_map, dtype='float32'), np.zeros_like(radio_map, dtype='float32')

print('Running cycles')
for x in range(len(radio_map[0]) - 1):
    print('x', x)
    for y in range(len(radio_map) - 1):

        if np.isnan(nh_map[y, x]):
            continue
        if np.isnan(radio_map[y, x]):
            continue
        # making distribution from nh data and fractional error
        nh = np.random.normal(nh_map[y, x] * 1e22, nh_frac_error * nh_map[y, x] * 1e22, num_sim)
        av_from_nh = nh_av(nh)

        # extinction map
        A_36 = A36_Av_extc * av_from_nh

        # extinction_A36[y][x] = np.nanmedian(A_36)
        # extinction_A36_err[y][x] = np.nanstd(A_36)

        # making distribution from Spitzer data and err files
        Spitzer_flux_simulation = np.random.normal(Spitzer_flux_map[y, x], Spitzer_flux_err[y, x], num_sim)

        # BKG correction
        Spitzer_fluxbkg_cor = Spitzer_flux_simulation - Spitzer_bkg_distribution
        Spitzer_fluxbkg_cor[Spitzer_fluxbkg_cor < 0] = np.nan

        # extinction correction
        Spitzer_flux_cor = ext_cor(Spitzer_fluxbkg_cor, A_36)

        # making distribution from radio data and fractional error
        radio_flux2009 = np.random.normal(radio_map[y, x], radio_frac_error * radio_map[y, x], num_sim)

        # Spectral index
        alpha = pl_index(Spitzer_flux_cor, radio_flux2009, vnu2, vnu1)

        spitzer_flux_matrix[y][x] = Spitzer_flux_cor
        spitzer_flux_median[y][x] = np.nanmedian(Spitzer_flux_cor)
        spitzer_flux_err[y][x] = np.nanstd(Spitzer_flux_cor)

        radio_flux_matrix[y][x] = radio_flux2009
        radio_flux_median[y][x] = np.nanmedian(radio_flux2009)
        radio_flux_err[y][x] = np.nanstd(radio_flux2009)

        alpha_matrix[y][x] = alpha
        alpha_matrix_median[y][x] = np.nanmedian(alpha)
        alpha_matrix_err[y][x] = np.nanstd(alpha)

print('Saving alpha matrix')
np.save(outpath + 'mc_alpha_R-IR_distr_matrix.npy', alpha_matrix)
# print('Saving Spitzer flux matrix')
# np.save(outpath + 'mc_spitzer_flux_distr_matrix.npy', spitzer_flux_matrix)
# print('Saving Radio flux matrix')
# np.save(outpath + 'mc_radio_flux_distr_matrix.npy', radio_flux_matrix)


print('Saving median maps')
# ALPHA MEDIAN
inputdata = alpha_matrix_median
comment = 'Spectral index map median values'
outputfile = outpath + 'alpha_R-IR_median_masked.fits'
write_into_fits(inputdata, comment, outputfile)

# ALPHA ERR
inputdata = alpha_matrix_err
comment = 'Spectral index map error values'
outputfile = outpath + 'alpha_R-IR_err_masked.fits'
write_into_fits(inputdata, comment, outputfile)

# RADIO MEDIAN
inputdata = radio_flux_median
comment = 'Radio Flux median values'
outputfile = outpath + 'radio_flux_median_masked.fits'
write_into_fits(inputdata, comment, outputfile)

# SPITZER MEDIAN
inputdata = spitzer_flux_median
comment = 'Spitzer Flux corrected median values'
outputfile = outpath + 'spitzer_flux_median_masked.fits'
write_into_fits(inputdata, comment, outputfile)


# RADIO ERROR
inputdata = radio_flux_err
comment = 'Radio Flux error values'
outputfile = outpath + 'radio_flux_err_masked.fits'
write_into_fits(inputdata, comment, outputfile)

# SPITZER ERROR
inputdata = spitzer_flux_err
comment = 'Spitzer Flux corrected error values'
outputfile = outpath + 'spitzer_flux_err_masked.fits'
write_into_fits(inputdata, comment, outputfile)

print('Done')
