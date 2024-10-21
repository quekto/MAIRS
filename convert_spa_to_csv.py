# ==============
# Version 2.1.0
# Created by Ian Fitch Mochida at University of Tokyo
# 2024-10-18 19:31 - 2.1.0 : Modified the data structure from accessing directories every time from holding a tree of the initial directory structure.
# 2024-10-18 19:31 - 2.1.1 : Counter modification.
# 2024-10-21 - 2.1.2 : Fixed bug where output within LoadSpectrum of wavenumber was not wavenumber.
# ==============

import argparse
import os
import pathlib
import sys
import time
from pathlib import Path, PurePath

import pandas as pd

from LoadSpectrum import read_spa

cnt: int = 0
file_count: int = 0

def parse_new() -> pathlib.PosixPath: 
    global file_count
    # Create the parser
    parser = argparse.ArgumentParser(description="Process some directory.")
    # Add the positional argument for the target directory
    parser.add_argument('directory', nargs='?', default='.', help='Target directory (default: current directory)')
    # Parse the command-line arguments
    args = parser.parse_args()
    # Get the directory argument
    basepath:pathlib.PosixPath = Path(args.directory) # Parent directory

    file_count =len([str(x) for x in list(pathlib.Path(basepath).rglob('*.spa'))] + \
        [str(x) for x in list(pathlib.Path(basepath).rglob('*.SPA'))])
    print(f'Number of spa files detected: {file_count}')

    return basepath

def write_to_csv(path: str):
    global cnt, file_count
    spectra_tmp, wavenumber_tmp, title_tmp = read_spa(path)
    df = pd.DataFrame({"wavenumber": wavenumber_tmp, "spectra": spectra_tmp})

    # Set output path
    # Remove extension with `splitext` and change it to csv.
    csv_path =  PurePath(path).parent.joinpath(PurePath(path).stem + '.csv')

    # Convert df -> csv
    df.to_csv(csv_path, index=False)

    cnt += 1
    if cnt % 10 == 0:
        sys.stdout.write(f'\r{cnt}/{file_count} files processed.'.ljust(30))
        sys.stdout.flush()

    if PurePath(path).stem[-4:] == "IPow":
        csv_path_IPOP:str = path[:-8] + 'IPOPow.csv'

        spa_path_IP:str   = path[:-8] + 'IPow.spa'
        spa_path_OP:str   = path[:-8] + 'OPow.spa'
        
        # IP
        spectra_tmp, wavenumber_tmp, title_tmp = read_spa(spa_path_IP)
        df_IP = pd.DataFrame({"wavenumber": wavenumber_tmp, "spectra": spectra_tmp})
        # OP
        spectra_tmp, wavenumber_tmp, title_tmp = read_spa(spa_path_OP)
        df_OP = pd.DataFrame({"wavenumber": wavenumber_tmp, "spectra": spectra_tmp})

        df_IPOP = pd.concat([df_IP, df_OP], axis=1)
        df_IPOP.to_csv(csv_path_IPOP, index=False)

        cnt += 1
        if cnt % 10 == 0:
            sys.stdout.write(f'\r{cnt}/{file_count} files processed.'.ljust(30))
            sys.stdout.flush()

def recursive(tree: dict, parent_dir:str) -> None:
    ipop_flag: int = 0
    for child in tree:
        dir = parent_dir + '/' + child
        if tree[child]: # If its a directory, True.
            recursive(tree[child], dir)
        elif PurePath(child).suffix not in ['.spa', '.SPA']: # Skip file other than spa
            continue
        else:
            write_to_csv(dir)

def build_directory_tree(root_dir) -> dict:
    tree = {}
    for entry in os.listdir(root_dir):
        path = os.path.join(root_dir, entry)
        if os.path.isdir(path):
            tree[entry] = build_directory_tree(path)  # Recursively build tree for subdirectories
        else:
            tree[entry] = None  # Files are leaf nodes
    return tree

def main():
    global cnt
    basepath: pathlib.PosixPath = parse_new() # 大元のディレクトリ
    tree = build_directory_tree(basepath)
    recursive(tree, str(basepath)) # 再帰的に下っていき、spaファイルにたどり着いたら変換

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
