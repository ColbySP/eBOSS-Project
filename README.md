# eBOSS-Project

The code is assocaited with the eBOSS research project underway at Bates College.

While potentially not the most simply or pretty implimentation takes in a few parameters and returns
a pandas DataFrame of stacked spectrum. For error propagation we perturb the initial data slightly
using the variance measure provided in each spectra's fits file. The code allows for the user to
input the number of simulations desired and returns that many columns for further uncertainty
calculations. The first column is always the unperturbed data, while others are perturbed using
gaussian noise. Other inputs are the reference directory for the galaxies to get a measure of
redshift. A directory where all the specta fits files are kept, a list of the file names to include
in the stacked procedure, and lastly the min and max wavelenght values to constrain stacking to.

If you have any issues with the code, email me at czsp32@gmail.com and I will be able to clarify.
