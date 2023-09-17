# nebular-extract-analyse
This repository contains code for extracting spectroscopic data from CSV files, correcting the fluxes and errors, calculating temperatures and densities using various line ratios, and calculating ionic and elemental abundances using various ICFs. The code outputs the results in a readable and familiar format.

The code uses the PyNeb package to perform all of the calculations and analysis. PyNeb is a Python package for analyzing astronomical spectroscopic data. It provides a variety of tools for extracting features from spectra, fitting spectral models, and computing chemical abundances.

To use the code, simply clone the repository. Then, you can run the code by passing the path to your CSV file and the name of the file that will be read by PyNeb as command-line arguments.

For example, to run the code on a CSV file called "my_data.csv" and output the extracted data to a file called "my_results.dat",which is then used by the code to compute all the calculations automatically.

You would use the following command: 'python3 dataE_A.py my_data.csv new_file.dat'

The code will output a new CSV file containing the extracted and corrected fluxes, in addtion to csv files containing ionic and elemental abundances, also it will output tables with computed tempratures and densities using the relevant line ratios.

This code is useful for astronomers who need to extract, correct, and analyze spectroscopic data from CSV files. It can also be used by students and researchers who are interested in learning more about spectroscopic data analysis.
