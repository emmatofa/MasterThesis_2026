# %% [markdown]
# # TADs: Insulation and boundaries in *A. rouxii* at 10 kb

# %% [markdown]
# ## Importing libraries

# %%
# Import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# %%
# Import python libraries for working with cooler files and tools for analysis
import cooler
import cooltools.lib.plotting
from cooltools import insulation # For TAD calling
import cooltools

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.4'):
    raise AssertionError("tutorials rely on cooltools version 0.5.4 or higher,"+
                         "please check your cooltools version and update to the latest")


# %%
# Import libraries for plotting 
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe
import matplotlib as mpl

# %% [markdown]
# ## Set format for file

# %%
plt.rcParams['font.family'] = 'Times New Roman'

# %% [markdown]
# ## Importing Hi-C matrix

# %%
cooler.fileops.list_coolers('/Users/emma/Documents/NMBU/Master/Master/Koding_python/ggAmyRoux1.mcool')

# %%
# Import the cooler file
resolution = 10000 # Use 10 kb resolution for the TAD calling
res_10kb = cooler.Cooler(f'/Users/emma/Documents/NMBU/Master/Master/Koding_python/ggAmyRoux1.mcool::/resolutions/{resolution}')

# %% [markdown]
# ## Filtering

# %%
# Filtering out scaffolds
min_size = 100_000
# Keep only chromosomes/scaffolds that are large enough
chroms_to_keep = [c for c in res_10kb.chromnames if res_10kb.chromsizes[c] >= min_size and "Scaffold" not in c]

print("Chromosomes kept:", chroms_to_keep) # Print the chromosomes that we plot. 

# %%
# Get the bins for the chromosomes we keep
bins_res_10kb= res_10kb.bins()[:]
keep_mask_res_10kb= bins_res_10kb['chrom'].isin(chroms_to_keep)
# Fetch full genome-wide matrix (balanced)
full_matrix_res_10kb = res_10kb.matrix(balance=True)[:]

# Subset matrix to only keep large chromosomes
matrix_filtered_res_10kb = full_matrix_res_10kb[keep_mask_res_10kb.values, :][:, keep_mask_res_10kb.values]

# %%
# Create view dataframe that excludes the filtered out scaffolds when calculating insulation
view_df = bioframe.make_viewframe(
    [(c, 0, res_10kb.chromsizes[c]) for c in chroms_to_keep]
)

# %% [markdown]
# ## TAD calling using `Insulation` from Cooltools

# %% [markdown]
# TAD calling was performed following the insulation score tutorial found at cooltools: https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html

# %%
resolution = 10000
windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution] # The window sizes are 30 kb, 50 kb, 100 kb and 250 kb. 
insulation_table = insulation(
    res_10kb,
    window_bp=windows,
    view_df=view_df,   # This excludes scaffolds and only uses chromosomes in the filtered data
    ignore_diags=2,
    verbose=True
) # Using insulation on the mcool file at window sizes set above. Verbose is used to get progress reports during running.

# %%
# Inspect the dataframe with insulation score, 
# number of valid pixels, if the bin is valid, the boundary strength and whether a locus is called a boundary after thresholding, 
# for each of the window sizes which are 30 kb, 50 kb, 100 kb and 250 kb. 
# Get the information for the largest window, here 250 000. 
first_window_summary =insulation_table.columns[[ str(windows[0]) in i for i in insulation_table.columns]]


insulation_table[['chrom','start','end','region','is_bad_bin']+list(first_window_summary)].iloc[1000:1005]

# %% [markdown]
# ## Functions used for plotting (found in the cooltools tutorial)

# %%
# Functions to help with plotting (found in the cooltools tutorial)
# Function for creating a 45 degrees rotated heatmap using pcolormesh found Matplotlib
# Creates a grid and reshapes it
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs): 
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)] 
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im

# %%
from matplotlib.ticker import EngFormatter # Used to formate axis labels
bp_formatter = EngFormatter('b') # Here used to formate to get genomic coordinates in megabases
def format_ticks(ax, x=True, y=True, rotate=True): # To format axis ticks 
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

# %% [markdown]
# ## Insulation track for one chromosome at all window sizes

# %%
start = 0 # Define the start position, start with 0 since coordinates for each chromosome start with 0.
end = start+ 90*windows[0] 
region = ('subgenom1_SUPER_7', start, end) # Plots the chromosome from start to set end

norm = LogNorm(vmax=0.1, vmin=0.001) # Log normalization of the contact map
data = res_10kb.matrix(balance=True).fetch(region) # Uses the mcool file and balances it and then extracts the region set above

f, ax = plt.subplots(figsize=(18, 6)) # Creates the plot and set the size of the plot
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall') # Uses the pcolormesh_45deg function to rotate and apply log normalization. Uses fall as colour in Hi-C matrix. 

ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0]) # limits the height of the triangle (contact map)
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False) # Hides x-axis for the contact map since it will share with the second plot

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)

insul_region = bioframe.select(insulation_table, region) # Extract the information for the region from the insulation_table

ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax) # Create the bottom axis, the size is 50% of main plot and shares x-axis with contact map
ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0,1,5))))

ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), # Creates the insulation plot. 
            insul_region['log2_insulation_score_'+str(windows[0])], # Plots the log_2_insulation score
            label=f'Window {windows[0]} bp') # Label showing window size in basepairs

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4); 

format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])

# Add the other windows 
for res in windows[1:]:
    ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), insul_region[f'log2_insulation_score_{res}'], label=f'Window {res} bp') # Add label showing window size for all lines
ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4)

ax.set_title(f'Insulation score and contact map (10 k) for chromosome 3 in A. rouxii for different window sizes')


# %% [markdown]
# ## Boundary calling

# %%
# Valleys in insulation score correspond to highly insulating regions and are called as TAD boundaries.
# The potential boundaries are assigned a boundary strength which are used to find regions which insulate strongly and are added in is_boundary column.
start = 0
end = start+ 90*windows[0]
region = ('subgenom1_SUPER_7', start, end)
norm = LogNorm(vmax=0.1, vmin=0.001)
data = res_10kb.matrix(balance=True).fetch(region)

f, ax = plt.subplots(figsize=(20, 10))
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)

insul_region = bioframe.select(insulation_table, region) # Get information from the insulation_table for the given region

ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax) # Create the bottom axis, the size is 50% of main plot and shares x-axis with contact map

ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), # Create the insulation score plot
            insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

# Identify boundaries 
boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])] # Removes rows with missing values (NA) in boundary strength (for window size 30 k)

# Divide into weak and strong boundaries 
weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']] # Weak boundaries when is_boundary_3000 = False
strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']] # Strong boundaries when is_boundary_3000 = True

ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1), # Add scatter points for weak boundaries 
            weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1), # Add scatter points for strong boundaries 
            strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4, fontsize = 12); # Legend below insulation plot


format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])


ax.set_title('Contact map and insulation score across chromosome 3 (10 kb and 30 k windows)', fontsize = 16)

# %% [markdown]
# ## Plotting the whole chromosome 

# %%
# Print the chromosome sizes 
print(res_10kb.chromsizes)

# %%
chrom = 'subgenom1_SUPER_7'

start = 0
end = res_10kb.chromsizes[chrom]
region = (chrom, start, end)
norm = LogNorm(vmax=0.1, vmin=0.001)
data = res_10kb.matrix(balance=True).fetch(region)

f, ax = plt.subplots(figsize=(40, 15))


im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)
ax.tick_params(axis='y', labelsize=16)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)

insul_region = bioframe.select(insulation_table, region) # Get information from the insulation_table for the given region

ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax) # Create the bottom axis, the size is 50% of main plot and shares x-axis with contact map

ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), # Create the insulation score plot
            insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

# Identify boundaries 
boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])] # Removes rows with missing values (NA) in boundary strength (for window size 30 k)

# Divide into weak and strong boundaries 
weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']] # Weak boundaries when is_boundary_3000 = False
strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']] # Strong boundaries when is_boundary_3000 = True

ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1), # Add scatter points for weak boundaries 
            weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1), # Add scatter points for strong boundaries 
            strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4, fontsize = 20); # Legend below insulation plot
ins_ax.set_ylabel('Insulation Score\n(Log2)', fontsize = 20)
ins_ax.set_xlabel('Genomic position', fontsize = 20)

ins_ax.tick_params(axis='y', labelsize=16)
ins_ax.tick_params(axis='x', labelsize=16)
format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])


ax.set_title('Contact map and insulation score across chromosome 7 (10 kb and 30 k windows)', fontsize=30)

# %% [markdown]
# ## Adding coverage

# %%
# Import file with Hi-C coverage at 10 kb windows
coverage_file = "/Users/emma/Documents/NMBU/Master/Master/Results_all_fungi/ggAmyRoux1_coverage_10kb.txt"  
coverage_table = pd.read_csv(
    coverage_file,
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "coverage"]
)

# %%
# Sort so only chromosomes without scaffold is keept
coverage_table = coverage_table[coverage_table['chrom'].isin(chroms_to_keep)]

coverage_table = coverage_table.sort_values(['chrom', 'start']) # Sort the table by size

# %%
print(coverage_table.head())

# %% [markdown]
# ### Plot insulationscore together with coverage to look for correlation between low coverage and boundaries

# %%
chrom = 'subgenom1_SUPER_7' 
plt.rcParams['font.size'] = 16 # Set font size in plott 

# Plot the whole chromosome
start = 0
end = res_10kb.chromsizes[chrom]
region = (chrom, start, end)
norm = LogNorm(vmax=0.1, vmin=0.001)

data = res_10kb.matrix(balance=True).fetch(region) # Get the region from the Hi-C matrix


# Create the plot
f, ax = plt.subplots(figsize=(30, 15)) # Plot size
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall') # Uses the pcolormesh_45deg function to plot the hic map
ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)


# Hi-C colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)


# Insulation plot

insul_region = bioframe.select(insulation_table, region) # Get information from the insulation_table for the given region

ins_ax = divider.append_axes("bottom", size="50%", pad=0.01, sharex=ax) # Create the bottom axis, the size is 50% of main plot and shares x-axis with contact map

ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), # Create the insulation score plot
            insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

# Identify boundaries 
boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])] # Removes rows with missing values (NA) in boundary strength (for window size 30 k)

# Divide into weak and strong boundaries 
weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']] # Weak boundaries when is_boundary_3000 = False
strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']] # Strong boundaries when is_boundary_3000 = True

ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1), # Add scatter points for weak boundaries 
            weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1), # Add scatter points for strong boundaries 
            strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4); # Legend below insulation plot

ins_ax.set_ylabel('Insulation score \n (log2)', labelpad=20)
format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])


# Coverage plot 

cov_ax = divider.append_axes("bottom", size="50%", pad=0.05, sharex=ax) # Creates the plot below insulation score, shares x-axis with the other plots

# Gets coverage for given chromosome region
cov_region = coverage_table[
    (coverage_table['chrom'] == region[0]) &
    (coverage_table['end'] > region[1]) &
    (coverage_table['start'] < region[2])
]

x_cov = cov_region[['start','end']].mean(axis=1) # get the midpoint of the genomic interval, meaning the point is in the middle of the bin
y_cov = cov_region['coverage'] # Selects the coverage for that region 

# Draw lines connecting points
cov_ax.plot(x_cov, y_cov, color='green', lw=1.5, label='Coverage')
# Adds dots
cov_ax.scatter(x_cov, y_cov, color='green', s=20)

# Axis settings
cov_ax.set_ylabel("Coverage \n per 10 kb ", labelpad=3)
cov_ax.set_xlabel("Genomic position (Mb)")
format_ticks(cov_ax, y=True, rotate=False)
cov_ax.set_ylim(0, max(y_cov)*1.2)


# Add title and show plot
ax.set_title('Hi-C contact map, insulation score (30 kb), and coverage across chromosome 7 in A. rouxii', fontsize=30)
plt.tight_layout()
plt.show()


# %%
chrom = 'subgenom1_SUPER_7' 

plt.rcParams['font.size'] = 18 # Set font size in plott 

# Plot the whole chromosome
start = 0
end = 3000000
region = (chrom, start, end)
norm = LogNorm(vmax=0.1, vmin=0.001)

data = res_10kb.matrix(balance=True).fetch(region) # Get the region from the Hi-C matrix


# Create the plot
f, ax = plt.subplots(figsize=(30, 15)) # Plot size
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall') # Uses the pcolormesh_45deg function to plot the hic map
ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)


# Hi-C colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)


# Insulation plot

insul_region = bioframe.select(insulation_table, region) # Get information from the insulation_table for the given region

ins_ax = divider.append_axes("bottom", size="50%", pad=0.01, sharex=ax) # Create the bottom axis, the size is 50% of main plot and shares x-axis with contact map

ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), # Create the insulation score plot
            insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

# Identify boundaries 
boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{windows[0]}'])] # Removes rows with missing values (NA) in boundary strength (for window size 30 k)

# Divide into weak and strong boundaries 
weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']] # Weak boundaries when is_boundary_3000 = False
strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']] # Strong boundaries when is_boundary_3000 = True

ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1), # Add scatter points for weak boundaries 
            weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1), # Add scatter points for strong boundaries 
            strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4); # Legend below insulation plot

ins_ax.set_ylabel('Insulation score \n (log2)', labelpad=20)
format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])
ins_ax.xaxis.set_visible(False)

# Coverage plot 

cov_ax = divider.append_axes("bottom", size="50%", pad=0.5, sharex=ax) # Creates the plot below insulation score, shares x-axis with the other plots

# Gets coverage for given chromosome region
cov_region = coverage_table[
    (coverage_table['chrom'] == region[0]) &
    (coverage_table['end'] > region[1]) &
    (coverage_table['start'] < region[2])
]

x_cov = cov_region[['start','end']].mean(axis=1) # get the midpoint of the genomic interval, meaning the point is in the middle of the bin
y_cov = cov_region['coverage'] # Selects the coverage for that region 

# Draw lines connecting points
cov_ax.plot(x_cov, y_cov, color='green', lw=1.5, label='Coverage')
# Adds dots
cov_ax.scatter(x_cov, y_cov, color='green', s=20)

# Axis settings
cov_ax.set_ylabel("Coverage \n per 10 kb", labelpad=3)
cov_ax.set_xlabel("Genomic position (Mb)")
format_ticks(cov_ax, y=True, rotate=False)
cov_ax.set_ylim(0, max(y_cov)*1.2)


# Add title and show plot
ax.set_title('Hi-C contact map, insulation score (30 kb), and coverage in a region of chromosome 7 in M. flavus', fontsize=30)
plt.tight_layout()
plt.show()


# %% [markdown]
# ## Plot TADs as an overlay on the Hi-C matrix

# %%
# Function for extracting TADs from insulation table using adjacent boundaries. Found at cooltools
def extract_TADs(insulation_table, is_boundary_col, max_TAD_length=3_000_000):
    tads = bioframe.merge(insulation_table[insulation_table[is_boundary_col] == False]) # Use is_boundary in insulation table
    return tads[(tads["end"] - tads["start"]) <= max_TAD_length].reset_index(drop=True)[['chrom','start','end']]


# %%
TADs_table = extract_TADs(insulation_table, 'is_boundary_30000') # Extract TADs using is_boundary in 30 kb windows
TADs_table.head() # Shows the extracted TADs

# %% [markdown]
# ## Find sizes of TADs

# %%
# Calculate the size of each TAD in kb
TADs_table['TAD_size_kb'] = (TADs_table['end'] - TADs_table['start']) / 1000

# The total number of called TADs
total_TADs = len(TADs_table)
print("Total number of called TADs:", total_TADs)

# Smallest TAD
min_tad = TADs_table.loc[TADs_table['TAD_size_kb'].idxmin()]
print("Smallest TAD:")
print(min_tad[['chrom','start','end','TAD_size_kb']])

# Largest TAD
max_tad = TADs_table.loc[TADs_table['TAD_size_kb'].idxmax()]
print("\nLargest TAD:")
print(max_tad[['chrom','start','end','TAD_size_kb']])

# Mean size of TADs
mean_TAD_size = TADs_table['TAD_size_kb'].mean()
print("Mean TAD size (kb):", mean_TAD_size)


# %% [markdown]
# # Plotting TADs

# %%
# Visualizing the first 10 inter-boundary intervals as a grey overlay vs the Hi-C data
region = ('subgenom1_SUPER_7', TADs_table.iloc[0].start, TADs_table.iloc[10].end) # start of fist TAD and end of the 10th TAD in chromosome 3
norm = LogNorm(vmax=0.1, vmin=0.001) # Log scale for contact map
data = res_10kb.matrix(balance=True).fetch(region) # fetches the balanced contact matrix from the cooler file

f, ax = plt.subplots(figsize=(18, 15)) # creates the plot
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall') # plots the contact matrix rotated 45 degrees

# Set axes 
ax.set_aspect(0.5) 
ax.set_ylim(0, 13*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(True)

# Add titile
ax.set_title(f"First 10 TADs in Chromosome 3 ({resolution} kb, 15 000 window size)", fontsize=16)

# Adds the colorbar in the heatmap 
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)

# makes an overlay of the first 10 TADs with grey boxes
idx = 10
max_pos = TADs_table[:idx]['end'].max()/resolution
contact_matrix = np.zeros((int(max_pos), int(max_pos)))
contact_matrix[:] = np.nan


for _, row in TADs_table[:idx].iterrows():
    contact_matrix[int(row['start']/resolution):int(row['end']/resolution), int(row['start']/resolution):int(row['end']/resolution)] = 1
    contact_matrix[int(row['start']/resolution + 1):int(row['end']/resolution - 1), int(row['start']/resolution + 1):int(row['end']/resolution - 1)] = np.nan

# Shows the lines over the contact map in blue 
im = pcolormesh_45deg(ax, contact_matrix, start=0, resolution=resolution, cmap='Blues', vmax=1, vmin=-1, alpha=0.3)

plt.show()

# %% [markdown]
# ## Add TAD lines on top of coverage and insulation plot

# %%

f, ax = plt.subplots(figsize=(30, 15))
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)

# Add colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)


# Overlay first 10 TADs in grey
idx = 10  # first 10 TADs
max_pos = TADs_table[:idx]['end'].max() // resolution
contact_matrix = np.full((int(max_pos), int(max_pos)), np.nan)

for _, row in TADs_table[:idx].iterrows():
    start_idx = int(row['start']/resolution)
    end_idx = int(row['end']/resolution)
    # Outer border = 1, inner area transparent
    contact_matrix[start_idx:end_idx, start_idx:end_idx] = 1
    contact_matrix[start_idx+1:end_idx-1, start_idx+1:end_idx-1] = np.nan

# Plot overlay directly on the existing Hi-C heatmap
im_tads = pcolormesh_45deg(ax, contact_matrix, start=0, resolution=resolution,
                           cmap='Blues', vmax=1, vmin=-1, alpha=0.3)


# Add insulation and coverage
## Insulation plot
insul_region = bioframe.select(insulation_table, region)
ins_ax = divider.append_axes("bottom", size="50%", pad=0.01, sharex=ax)
ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
            insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

boundaries = insul_region.dropna(subset=[f'boundary_strength_{windows[0]}'])
weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']]
strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']]

ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
               weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
               strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4)
ins_ax.set_ylabel('Insulation score \n (log2)', labelpad=20)
format_ticks(ins_ax, y=False, rotate=False)
ins_ax.xaxis.set_visible(False)

### Add grid and vertical lines at strong boundaries in insulation plot
ins_ax.grid(True, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=0)
ins_ax.vlines(x=strong_boundaries[['start', 'end']].mean(axis=1),
              ymin=ins_ax.get_ylim()[0], ymax=ins_ax.get_ylim()[1],
              color='black', linestyle='--', linewidth=1, alpha=0.7, zorder=1)


## Coverage plot
cov_ax = divider.append_axes("bottom", size="50%", pad=0.5, sharex=ax)
cov_region = coverage_table[(coverage_table['chrom'] == region[0]) &
                            (coverage_table['end'] > region[1]) &
                            (coverage_table['start'] < region[2])]

x_cov = cov_region[['start','end']].mean(axis=1)
y_cov = cov_region['coverage']

cov_ax.plot(x_cov, y_cov, color='mediumpurple', lw=1.5, label='Coverage')
cov_ax.scatter(x_cov, y_cov, color='mediumpurple', s=20)
cov_ax.set_ylabel("Coverage \n per 10 kb", labelpad=3)
cov_ax.set_xlabel("Genomic position")
format_ticks(cov_ax, y=True, rotate=False)
cov_ax.set_ylim(0, max(y_cov)*1.2)

### Add grid and vertical lines at strong boundaries in coverage plot
cov_ax.grid(True, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=0)
cov_ax.vlines(x=strong_boundaries[['start', 'end']].mean(axis=1),
              ymin=cov_ax.get_ylim()[0], ymax=cov_ax.get_ylim()[1],
              color='black', linestyle='--', linewidth=1, alpha=0.7, zorder=1)


# Title and layout
ax.set_title('Hi-C contact map with first 10 TADs overlay, insulation score (log2), and coverage in chromosome 7 of A. rouxii', fontsize=30)
plt.tight_layout()
plt.show()


# %% [markdown]
# ## Plot TADs with gene density (calculated using 10 kb windows)

# %%
# Import gene density file for gene density per 10 kb window
gene_density= pd.read_csv('/Users/emma/Documents/NMBU/Master/Master/Data/Gene annotation/gene_density_ggAmyRoux1_10000.bed', 
                        sep='\t', header=None, names=['chrom','start','end', 'gene_count'])

# %%
# Plot TADs with gene density
region = ('subgenom1_SUPER_1', TADs_table.iloc[0].start, TADs_table.iloc[9].end)

f, ax = plt.subplots(figsize=(30, 15))
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)

# Add colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)

# Overlay first 10 TADs in blue
idx = 10  # first 10 TADs
max_pos = TADs_table[:idx]['end'].max() // resolution
contact_matrix = np.full((int(max_pos), int(max_pos)), np.nan)

for _, row in TADs_table[:idx].iterrows():
    start_idx = int(row['start']/resolution)
    end_idx = int(row['end']/resolution)
    # Outer border = 1, inner area transparent
    contact_matrix[start_idx:end_idx, start_idx:end_idx] = 1
    contact_matrix[start_idx+1:end_idx-1, start_idx+1:end_idx-1] = np.nan

# Plot overlay directly on the existing Hi-C heatmap
im_tads = pcolormesh_45deg(ax, contact_matrix, start=0, resolution=resolution,
                           cmap='Blues', vmax=1, vmin=-1, alpha=0.3)


# Insulation plot
insul_region = bioframe.select(insulation_table, region)
ins_ax = divider.append_axes("bottom", size="50%", pad=0.01, sharex=ax)
ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
            insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

boundaries = insul_region.dropna(subset=[f'boundary_strength_{windows[0]}'])
weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']]
strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']]

ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
               weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
               strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4)
ins_ax.set_ylabel('Insulation score\n (log2)', labelpad=12)
format_ticks(ins_ax, y=False, rotate=False)
ins_ax.xaxis.set_visible(False)

ins_ax.grid(True, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=0)
# Vertical lines at strong boundaries
ins_ax.vlines(x=strong_boundaries[['start', 'end']].mean(axis=1),
              ymin=ins_ax.get_ylim()[0], ymax=ins_ax.get_ylim()[1],
              color='black', linestyle='--', linewidth=1, alpha=0.7, zorder=1)

# Gene_density plot
gene_ax = divider.append_axes("bottom", size="50%", pad=0.5, sharex=ax)
gene_region = gene_density[(gene_density['chrom'] == region[0]) &
                            (gene_density['end'] > region[1]) &
                            (gene_density['start'] < region[2])]

x_gene = gene_region[['start','end']].mean(axis=1)
y_gene = gene_region['gene_count']

gene_ax.plot(x_gene, y_gene, color='mediumpurple', lw=1.5, label='Gene')
gene_ax.scatter(x_gene, y_gene, color='mediumpurple', s=20)
gene_ax.set_ylabel("Gene Density \n (gene count)", labelpad=25)
gene_ax.set_xlabel("Genomic position (Mb)")
gene_ax.set_ylim(0, max(y_gene)*1.2)

## Add boundary lines 
gene_ax.vlines(
    strong_boundaries[['start','end']].mean(axis=1),
    ymin=0,
    ymax=max(y_gene),
    color="black",
    linestyle="--",
    alpha=0.6
)

gene_ax.set_xlim(region[1], region[2])
gene_ax.grid(True, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=0)

# Final title and layout
ax.set_title('Chromosome 1 in A. rouxii', fontsize=30)
plt.tight_layout()
plt.show()


# %% [markdown]
# # For results

# %%
# Set fontsizes for result plots (to be the same across species)
mpl.rcParams.update({
    "font.size": 22,
    "axes.titlesize": 30,
    "axes.labelsize": 24,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
})

# %%
region = ('subgenom1_SUPER_1', TADs_table.iloc[0].start, TADs_table.iloc[9].end)

f, ax = plt.subplots(figsize=(30, 15))
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
ax.set_aspect(0.5)
ax.set_ylim(0, 10*windows[0])
format_ticks(ax, rotate=False)
ax.xaxis.set_visible(False)

# Add colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="1%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)

# Overlay first 10 TADs in blue
idx = 10  # first 10 TADs
max_pos = TADs_table[:idx]['end'].max() // resolution
contact_matrix = np.full((int(max_pos), int(max_pos)), np.nan)

for _, row in TADs_table[:idx].iterrows():
    start_idx = int(row['start']/resolution)
    end_idx = int(row['end']/resolution)
    # Outer border = 1, inner area transparent
    contact_matrix[start_idx:end_idx, start_idx:end_idx] = 1
    contact_matrix[start_idx+1:end_idx-1, start_idx+1:end_idx-1] = np.nan

# Plot overlay directly on the existing Hi-C heatmap
im_tads = pcolormesh_45deg(ax, contact_matrix, start=0, resolution=resolution,
                           cmap='Blues', vmax=1, vmin=-1, alpha=0.3)


# Insulation plot
insul_region = bioframe.select(insulation_table, region)
ins_ax = divider.append_axes("bottom", size="50%", pad=0.01, sharex=ax)
ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
            insul_region[f'log2_insulation_score_{windows[0]}'], label=f'Window {windows[0]} bp')

boundaries = insul_region.dropna(subset=[f'boundary_strength_{windows[0]}'])
weak_boundaries = boundaries[~boundaries[f'is_boundary_{windows[0]}']]
strong_boundaries = boundaries[boundaries[f'is_boundary_{windows[0]}']]

ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
               weak_boundaries[f'log2_insulation_score_{windows[0]}'], label='Weak boundaries')
ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
               strong_boundaries[f'log2_insulation_score_{windows[0]}'], label='Strong boundaries')

ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4)
ins_ax.set_ylabel('Insulation score\n (log2)', labelpad=12)
format_ticks(ins_ax, y=False, rotate=False)
ins_ax.xaxis.set_visible(False)

ins_ax.grid(True, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=0)
# Vertical lines at strong boundaries
ins_ax.vlines(x=strong_boundaries[['start', 'end']].mean(axis=1),
              ymin=ins_ax.get_ylim()[0], ymax=ins_ax.get_ylim()[1],
              color='black', linestyle='--', linewidth=1, alpha=0.7, zorder=1)

# Gene_density plot
gene_ax = divider.append_axes("bottom", size="50%", pad=0.5, sharex=ax)
gene_region = gene_density[(gene_density['chrom'] == region[0]) &
                            (gene_density['end'] > region[1]) &
                            (gene_density['start'] < region[2])]

x_gene = gene_region[['start','end']].mean(axis=1)
y_gene = gene_region['gene_count']

gene_ax.plot(x_gene, y_gene, color='mediumpurple', lw=1.5, label='Gene')
gene_ax.scatter(x_gene, y_gene, color='mediumpurple', s=20)
gene_ax.set_ylabel("Gene Density \n (10 kb windows)", labelpad=25)
gene_ax.set_xlabel("Genomic position (Mb)")
gene_ax.set_ylim(0, max(y_gene)*1.2)

## Add boundary lines 
gene_ax.vlines(
    strong_boundaries[['start','end']].mean(axis=1),
    ymin=0,
    ymax=max(y_gene),
    color="black",
    linestyle="--",
    alpha=0.6
)

gene_ax.set_xlim(region[1], region[2])
gene_ax.grid(True, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=0)

# Final title and layout
ax.set_title('Chromosome 1 in A. rouxii', fontsize=30)
plt.tight_layout()
plt.show()



