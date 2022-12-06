# required imports
from scipy.interpolate import interp1d
from datetime import timedelta
from astropy.io import fits
from tqdm import tqdm
import pandas as pd
import numpy as np
import time


def create_nonlinear_grid(x_min, x_max):
    """
    Function to create nonlinear spaced grid in accordance with telescope resultion
    x_min: (float) staring vale of the grid
    x_max: (float) upper limit for the final grid
    returns: (array) nonlinear spaced grid
    """
    grid = [x_min]
    newest_val = 0
    while newest_val < x_max:
        newest_val = grid[-1] * (1 + (70 / 3e5))
        grid.append(newest_val)
    return grid


def generate_filenames(df):
    """
    Function generate the fits filename based on the identification data
    """
    plate = str(df['PLATE']).zfill(4)
    mjd = str(df['MJD']).zfill(5)
    fiber = str(df['FIBER']).zfill(4)
    id = f'spec-{plate}-{mjd}-{fiber}.fits'
    return id


def create_simulations(wavelengths, min_wave, max_wave, fluxes, sigmas, n):
    """
    Function to take in spectra data and return many simulations of stacked spectra
    wavelengths (array): array of the wavelengths of each flux measurement
    min_wave (float): minimum wavelength to stack from
    max_wave (float): minimum wavelength to stack to
    fluxes (array): series of spectra data with the wavelength as index
    sigma (array): standard deviation values for simulations
    n (int): number of simulations to return
    Returns: DataFrame with many columns each corresponding to a simulation of stacked data
    """

    x_index = create_nonlinear_grid(min_wave, max_wave)  # widest possible range of wavelengths

    # create storage for simulations during each loop
    simulations = {}
    for i in range(n):
        means = np.full(sigmas.shape, 0)  # means are all going to be 0

        # make first simulation the true data
        if i == 0:
            noise = np.full(sigmas.shape, 0)
        else:
            noise = np.random.normal(means, sigmas)  # generate random numbers

        clipped_noise = np.clip(noise, -10, 10)  # clip to be +/- 10 at most
        perturbed_fluxes = fluxes + clipped_noise

        # now we clean and interpolate the data as before
        norm_flux = perturbed_fluxes / np.median(perturbed_fluxes)

        # interpolate the data
        model = interp1d(wavelengths, norm_flux, 'cubic', bounds_error=False)

        # add the interpolation to the larger DataFrame
        simulations[f'sim {i + 1}'] = model(x_index)

    # combine all the spectra into a single DataFrame
    sim_df = pd.DataFrame.from_dict(simulations)

    # format the DataFrame slightly better
    sim_df.index = x_index
    sim_df.index.names = ['WAVE']

    return sim_df


def prepare_all_spectra(ref_df, folder_dir, files, min_wave, max_wave, n_simulations):
    """
    Function to load fits files from a folder, before cleaning and returning a DataFrame
    ref_df (DataFrame)): the reference fits file loaded as a DataFrame
    folder_dir (str): directory to where all the fits files are stored
    files (list -> str): all files that are to be loaded
    min_wave (float): minimum wavelength to stack from
    max_wave (float): minimum wavelength to stack to
    n_simulations (int): number of simulations to create
    Returns: DataFrame of all spectral data
    """
    # start a timer
    start = time.time()

    # instantiate DataFrame to hold cleaned data
    dict_of_sims = {}

    # run loop to load and interpolate each galaxy's spectra
    for fits_file in tqdm(files, desc='loading and prepping data'):
        try:
            # clean and prepare the data
            galaxy_data = fits.open(f"{folder_dir}/{fits_file}")
            Table = pd.DataFrame(galaxy_data[1].data)
            Table['Z'] = ref_df.loc[fits_file, 'Z']
            Table['LAMBDA'] = np.power(10, Table['LOGLAM']) / (1 + Table['Z'])

            # set up normalized sigma and normalized flux
            Table['SIGMA'] = np.sqrt(1 / Table['IVAR'])

            # create slightly perturbed simulations of the spectral data
            sims = create_simulations(wavelengths=Table['LAMBDA'].values,
                                      min_wave=min_wave,
                                      max_wave=max_wave,
                                      fluxes=Table['FLUX'].values,
                                      sigmas=Table['SIGMA'].values,
                                      n=n_simulations)

            dict_of_sims[fits_file] = sims

        # if for some reason an error occurs skip it
        except:
            continue

    # stop the timer
    end = time.time()
    print(f"\nloaded all spectra data in {timedelta(seconds=end - start)}!")
    return dict_of_sims


def stack_spectra(ref_dir: str, folder_dir: str, files: list[str], min_wave: float, max_wave: float, n_simulations=1) -> 'pd.DataFrame':
    """
    Function to stack spectra and perturb spectra for future error analysis
    :param ref_dir: path to the reference file including plate and redshift information (includes file name)
    :param folder_dir: path to the folder containing spectra fits files (includes folder name)
    :param files: list of file names that should be included in the stacking
    :param min_wave: leftmost wavelength to stack from
    :param max_wave: rightmost wavelength to stack to
    :param n_simulations: the number of simulations or columns of returned stacked spectra
    :return: DataFrame with index of wavelengths and columns of stacked spectra (first is normal, others are perturbed)
    """

    # load the reference table and convert important columns to DataFrame
    print('loading reference table ...\n')
    ref_Table = fits.open(ref_dir)[1].data
    ref_df = pd.DataFrame(data={'PLATE': ref_Table['PLATE'].astype(str),
                                'MJD': ref_Table['MJD'].astype(str),
                                'FIBER': ref_Table['FIBER'].astype(str),
                                'Z': ref_Table['Z'].astype(float)})
    ref_df['FILENAME'] = ref_df.apply(generate_filenames, axis=1)
    ref_df.set_index('FILENAME', inplace=True)

    # load files and stack their spectra
    data = prepare_all_spectra(ref_df, folder_dir, files, n_simulations=n_simulations)

    # loop through each dictionary key and take the median value across them all
    x_axis = create_nonlinear_grid(min_wave, max_wave)
    cols = len(x_axis)
    rows = len(data.keys())

    # calculate the nan-median for each group of sims
    sim_dict = {}
    for i in range(0, n_simulations):

        # instantiate an array of 0's
        arr = np.full((cols, rows), np.nan)

        # fill each row with a galaxies ith sim
        for n, gal in enumerate(data.keys()):
            arr[:, n] = data[gal][f'sim {i + 1}'].values

        stack = np.nanmedian(arr, axis=1)
        sim_dict[f'stack_sim_{i + 1}'] = stack

    simulations = pd.DataFrame.from_dict(sim_dict)
    simulations.index = x_axis
    simulations.index.names = ['WAVE']
    return simulations


# set initial parameters and paths for stacking code
ref_dir = ...
folder_dir = ...
files = ...
n_simulations = 50

# run the code and return a DataFrame for further analysis
df = stack_spectra(ref_dir=ref_dir, folder_dir=folder_dir, files=files, min_wave=3600, max_wave=10400, n_simulations=1)
