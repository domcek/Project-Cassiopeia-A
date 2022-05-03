# Data Reproduction Package for “Mapping the spectral index of Cassiopeia A: evidence for flattening from radio to infrared“

V. Domček*, J. Vink, J.V. Hernández Santisteban, T. DeLaney and P. Zhou

*Corresponding author: vdomcek@gmail.com, OrcID: 0000-0002-7770-4538

Published in Monthly Notices of the Royal Astronomical Society: https://academic.oup.com/mnras/article/502/1/1026/6126073

Available on ArXiv: https://arxiv.org/abs/2005.12677

## Software

- Operating system: Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-64-generic x86_64)
- Software used in this work is publicly available at:
    - Anaconda / Python v3.7.4, including Astropy v3.2.2, NumPy v1.17.2, Matplotlib v3.3.1, SciPy v1.3.1, reproject v0.6, photutils v0.7.2, jupyterlab v1.1.4
    - Anaconda / Python v2.7.17 with pyraf v2.1.15 and its dependencies (https://pypi.org/project/pyraf/)

## Data

All data files required to run the project are accessible in the Zenodo repository: https://doi.org/10.5281/zenodo.4478615

## Scripts

Instruction for the order to run the project is provided in bash file “./Run.sh”. Provided all python packages are installed at correct version and zenodo folder structure is preserved, the whole project can be run by executing only “./Run.sh”.

Basic details on the function of the script is documented in the beginning of each script.

## Figures and tables

All figures of the paper are an output of Jupyter lab notebook “zenodo_casa_figures.ipynb”. Figure 9 is produced by Jupyter lab notebook “zenodo_casa_SED.ipynb”, which also contains values for Table.1.

“./zenodo_casa_figures.ipynb”

“./zenodo_casa_SED.ipynb”
