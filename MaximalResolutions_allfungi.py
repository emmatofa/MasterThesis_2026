# Check maximal resolutions available

# Import packages
import cooler
import numpy as np

# Function for finding fraction of bins with more than 1000 contacts at different resolutions

def fraction_bins_with_contacts(clr, min_contacts=1000):
    pixels = clr.pixels()[:] # Contacts 
    bins = clr.bins()[:] # List of bins at given resolution
    coverage = np.bincount( 
        np.concatenate([pixels["bin1_id"], pixels["bin2_id"]]),
        weights=np.concatenate([pixels["count"], pixels["count"]]),
        minlength=len(bins)
    )
    # Computes the fraction of bins that have the minimum number of contacts (1000)
    fraction = (coverage >= min_contacts).sum() / len(bins)
    return fraction

# *A. rouxii*
## Path to mcool file
ARouxii_mcool_file = "/Users/emma/Documents/NMBU/Master/Master/Koding_python/ggAmyRoux1.mcool"

## list all available resolutions 
res_paths_ARouxii = cooler.fileops.list_coolers(ARouxii_mcool_file)

threshold_fraction = 0.8 ## at least 80% of bins should have at least 1000 contacts

## Loops through all resolutions in the mcool file 
for res_path in res_paths_ARouxii: 
    print(f"\nChecking resolution: {res_path}")
    ## Open the mcool file at the different resolutions
    clr = cooler.Cooler(f"{ARouxii_mcool_file}::{res_path}")
    ##Finds the fraction of bins that have the minimum number and prints the resolutions and the fraction using two decimals 
    fraction = fraction_bins_with_contacts(clr, min_contacts=1000)
    print(f"Resolution: {clr.binsize} bp, fraction ≥1000 contacts: {fraction:.2f}")

# *M. flavus*
## Path to mcool file
MFlavus_mcool_file = "/Users/emma/Documents/NMBU/Master/Master/gzMucFlav1/Subgenome2/gzMucFlav1_sub2.mcool"

## list all resolution paths
res_paths_MFlavus = cooler.fileops.list_coolers(MFlavus_mcool_file)

threshold_fraction = 0.8

for res_path in res_paths_MFlavus:
    print(f"\nChecking resolution: {res_path}")
    ## open the mcool group using the full path from list_coolers
    clr = cooler.Cooler(f"{MFlavus_mcool_file}::{res_path}")
    
    fraction = fraction_bins_with_contacts(clr, min_contacts=1000)
    print(f"Resolution: {clr.binsize} bp, fraction ≥1000 contacts: {fraction:.2f}")

# *M. racemosus*
# Path to mcool file
MRacemosus_mcool_file = "/Users/emma/Documents/NMBU/Master/Master/gzMucFlav1/nf-core:HiC pipeline files/gzMucRace1.mcool"

## list all resolution paths
res_paths_MRacemosus = cooler.fileops.list_coolers(MRacemosus_mcool_file)

threshold_fraction = 0.8

for res_path in res_paths_MRacemosus:
    print(f"\nChecking resolution: {res_path}")
    ## open the mcool group using the full path from list_coolers
    clr = cooler.Cooler(f"{MRacemosus_mcool_file}::{res_path}")
    
    fraction = fraction_bins_with_contacts(clr, min_contacts=1000)
    print(f"Resolution: {clr.binsize} bp, fraction ≥1000 contacts: {fraction:.2f}")
    
    if fraction >= threshold_fraction and map_resolution is None:
        map_resolution = clr.binsize
        print(f"Map resolution = {map_resolution} bp")