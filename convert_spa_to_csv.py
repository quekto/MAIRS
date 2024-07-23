from LoadSpectrum import read_spa
import pathlib
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import sys, time

def read_spa(filepath):
    '''
    Input
    Read a file `*.spa`
    ----------
    Output
    Return spectra, wavelenght (nm), titles
    '''
    with open(filepath, 'rb') as f:
        f.seek(564)
        Spectrum_Pts = np.fromfile(f, np.int32, 1)[0]
        f.seek(30)
        SpectraTitles = np.fromfile(f, np.uint8, 255)
        SpectraTitles = ''.join([chr(x) for x in SpectraTitles if x != 0])

        f.seek(576)
        Max_Wavenum = np.fromfile(f, np.single, 1)[0]
        Min_Wavenum = np.fromfile(f, np.single, 1)[0]
        # print(Min_Wavenum, Max_Wavenum, Spectrum_Pts)
        Wavenumbers = np.flip(np.linspace(
            Min_Wavenum, Max_Wavenum, Spectrum_Pts))

        f.seek(288)

        Flag = 0
        while Flag != 3:
            Flag = np.fromfile(f, np.uint16, 1)

        DataPosition = np.fromfile(f, np.uint16, 1)
        f.seek(DataPosition[0])

        Spectra = np.fromfile(f, np.single, Spectrum_Pts)
    return Spectra, Wavenumbers, SpectraTitles

basepath = '.'
# Get paths of .spa, .SPA files.
paths = [str(x) for x in list(pathlib.Path(basepath).rglob('*.spa'))] + \
    [str(x) for x in list(pathlib.Path(basepath).rglob('*.SPA'))]
print('Files detected: {}'.format(len(paths)))

length = len(paths)
counter = 0
for path in paths:
  counter += 1
  
  spectra_tmp, wavenumber_tmp, title_tmp = read_spa(path)
  df = pd.DataFrame({"wavenumber": wavenumber_tmp, "spectra": spectra_tmp})

  csv_path = os.path.splitext(path)[0] + '.csv'

  df.to_csv(csv_path, index=False)

  if counter % 100 == 0:
    sys.stdout.write(f"\r{counter} / {length} files read.")
    sys.stdout.flush()
    # time.sleep(0.5)
sys.stdout.write(f"\rAll {length} files read!")

# # Plot the data
# plt.figure(figsize=(8, 6))
# plt.plot(wavelength_tmp, spectra_tmp, marker='o', linestyle='-', color='b', label='Data')

# # Adding title and labels
# plt.title('X vs Y Plot')
# plt.xlabel('X Axis')
# plt.ylabel('Y Axis')

# # Adding grid
# plt.grid(True)

# # Adding legend
# plt.legend()

# # Show the plot
# plt.show()
