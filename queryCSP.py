import os
import glob
import time
import shutil
import requests
import tarfile
import numpy as np

def download(save_loc: str):
    """
    :param save_loc: Location to save tarball and 'DR3/' folder with all CSP-I Third Release data.
    :return: None
    """
    url = "https://csp.obs.carnegiescience.edu/data/CSP_Photometry_DR3.tgz"
    tar_name = "allCSP-I_tarball.tgz"
    path = save_loc+tar_name

    # Check if already downloded
    if os.path.exists(save_loc+tar_name):
        print("[+++] File previous downloaded! Skipping...")
    else:
        # Use request to get file
        response = requests.get(url)

        # Save file
        if response.status_code == 200:
            with open(path, "wb") as f:
                f.write(response.content)
            print(f"[+++] File '{tar_name}' downloaded successfully to {path}.")
        else:
            raise(f"[!!!] Failed to download the file. Status code: {response.status_code}")

    # Unpack tarball
    extract_path = save_loc
    if os.path.exists(extract_path+'DR3/'):
        print("[+++] File previous extracted! Skipping...")
    else:
        print(f"[+++] Extracting tarball to {extract_path}...")
        with tarfile.open(path, "r:gz") as tar:
            tar.extractall(extract_path)
    return
def seperate_data(tar_list: str, data_loc: str, save_loc: str):
    """
    :param tar_list: List of target names to seperate
    :param data_loc:
    :param save_loc:
    :return:
    """
    targets = np.genfromtxt(tar_list, delimiter=',', skip_header=1, dtype=str)
    hdr = list(targets[0])
    all_names = list(targets[:, hdr.index('Name')])[1:]

    # Sort through data
    print(f"[+++] Relocating 1991bg-like SNe from '{data_loc+'DR3/'}' to '{save_loc}'...")
    for path in glob.glob(data_loc+'DR3/*.txt'):
        name = 'SN ' + path.split('/')[-1].split('_')[0][2:]
        if name in all_names:
            shutil.copy(path, save_loc)
    return


if __name__ == "__main__":
    start = time.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(time.time() - start, 4), 'seconds\n|---------------------------|')
