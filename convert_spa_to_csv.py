# ==============
# Version 2.0.0
# Created by Ian Fitch Mochida at University of Tokyo
# 2024-10-18
# ==============

import argparse
import os
import pathlib
import sys
import time
from pathlib import Path, PurePath

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from LoadSpectrum import read_spa

cnt: int = 0

def parse_new() -> pathlib.PosixPath: 
    # Create the parser
    parser = argparse.ArgumentParser(description="Process some directory.")
    # Add the positional argument for the target directory
    parser.add_argument('directory', nargs='?', default='.', help='Target directory (default: current directory)')
    # Parse the command-line arguments
    args = parser.parse_args()
    # Get the directory argument
    basepath = Path(args.directory) # Parent directory

    file_count =len([str(x) for x in list(pathlib.Path(basepath).rglob('*.spa'))] + \
        [str(x) for x in list(pathlib.Path(basepath).rglob('*.SPA'))])
    print(f'Number of spa files detected: {file_count}')

    return basepath

def write_to_csv(path: pathlib.PosixPath, flag:int):
    global cnt
    spectra_tmp, wavenumber_tmp, title_tmp = read_spa(path)
    df = pd.DataFrame({"wavenumber": wavenumber_tmp, "spectra": spectra_tmp})

    # Set output path
    # Remove extension with `splitext` and change it to csv.
    csv_path =  PurePath(path).parent.joinpath(PurePath(path).stem + '.csv')

    # Convert df -> csv
    df.to_csv(csv_path, index=False)

    cnt += 1
    if flag == 2: # Concat IP and OP when both detected
        csv_path_IPOP:str = str(path)[:-8] + 'IPOPow.csv'

        spa_path_IP:str   = str(path)[:-8] + 'IPow.spa'
        spa_path_OP:str   = str(path)[:-8] + 'OPow.spa'
        
        # IP
        spectra_tmp, wavenumber_tmp, title_tmp = read_spa(spa_path_IP)
        df_IP = pd.DataFrame({"wavenumber": wavenumber_tmp, "spectra": spectra_tmp})
        # OP
        spectra_tmp, wavenumber_tmp, title_tmp = read_spa(spa_path_OP)
        df_OP = pd.DataFrame({"wavenumber": wavenumber_tmp, "spectra": spectra_tmp})

        df_IPOP = pd.concat([df_IP, df_OP], axis=1)
        df_IPOP.to_csv(csv_path_IPOP, index=False)

        cnt += 1

def recursive(path: pathlib.PosixPath) -> None:
    ipop_flag: int = 0
    for child in path.iterdir():
        if child.is_dir(): # First go down if still directory
            recursive(child)
        elif child.suffix not in ['.spa', '.SPA']: # Skip file other than spa
            continue
        else:
            if child.stem[-4:] in ['IPow', 'OPow']: # Check if the file is IPow or OPow
                ipop_flag += 1
            write_to_csv(child, ipop_flag)
            if ipop_flag == 2:
                ipop_flag = 0

    
    sys.stdout.write(f'\r{cnt} files processed.'.ljust(30))
    sys.stdout.flush()

def main():
    global cnt
    basepath: pathlib.PosixPath = parse_new() # 大元のディレクトリ
    recursive(basepath) # 再帰的に下っていき、spaファイルにたどり着いたら変換

    sys.stdout.write(f'\r{cnt} files created!'.ljust(30))
    sys.stdout.flush()

    return

if __name__ == "__main__":
    main()
    print() # Cleaning up stdout


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
# # Adding legend
# plt.legend()

# # Show the plot
# plt.show()
