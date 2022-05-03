# Flux calculation for the high resolution radio data
ipython 2_radio_calib.py 

# Reprojection of datasets to the common coordinate system
ipython 3_reproject_C10.py
ipython 3_reproject_L10.py
ipython 3_reproject_nh.py
ipython 3_reproject_radio.py
ipython 3_reproject_spitzer.py
ipython 3_reproject_spitzer_errors.py
ipython 3_reproject_spitzer_subpyraf_sm.py
ipython 3_reproject_SPIXBL10cor.py
ipython 3_reproject_SPIXERR10.py

# Creating masks
ipython 4_5_mask_radio_fraciton_map.py
ipython 5_mask_alpha_R.py
ipython 5_mask_nH.py
ipython 5_mask_radio.py
ipython 5_mask_spitzer.py
# additional masks are produced by zenodo_casa_figures.ipynb
# 5_mask_bin3_spitzer_stars_smooth.fits was created manually based on spitzer uncertainty map

# Monte Carlo process
ipython monte_carlo_high-res.py
ipython monte_carlo_low-res_smooth.py
ipython monte_carlo_diff.py

# Figures, table 1 values and additional region masks 
# are produced in zenodo_casa_figures.ipynb and zenodo_casa_SED.ipynb