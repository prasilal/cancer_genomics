import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the output of samtools depth
def read_depth_file(depth_file):
    df = pd.read_csv(depth_file, sep='\t', header=None, names=['chrom', 'pos', 'depth'])
    return df

# Calculate mean read depth for a window of 1000 bases
def calculate_mean_depth(df):
    window_size = 1000
    df['window'] = df['pos'] // window_size
    mean_depth = df.groupby(['chrom', 'window'])['depth'].mean().reset_index()
    mean_depth['pos'] = mean_depth['window'] * window_size + window_size // 2
    return mean_depth

# Calculate log2 ratio between tumor and wild-type
def calculate_log2_ratio(df_tumor, df_wildtype):
    merged_df = pd.merge(df_tumor, df_wildtype, on=['chrom', 'pos'], suffixes=('_tumor', '_wildtype'))
    merged_df['log2_ratio'] = np.log2(merged_df['depth_tumor'] / merged_df['depth_wildtype'])
    return merged_df

# Plot log2 ratio and save as image
def save_log2_ratio_plot(df):
    plt.figure(figsize=(10, 6))
    plt.plot(df['pos'], df['log2_ratio'], color='blue', linestyle='-', marker='.', markersize=2)
    plt.xlabel('Position in genome')
    plt.ylabel('Log2 Ratio (Tumor/Wild-type)')
    plt.title('Read-Depth Plot')
    plt.grid(True)
    plt.savefig('read_depth_plot.png', bbox_inches='tight')  # Save plot as image file
    plt.close()


tumor_depth_file = 'tu_read_depth.txt'
wildtype_depth_file = 'wt_read_depth.txt'

df_tumor = read_depth_file(tumor_depth_file)
df_wildtype = read_depth_file(wildtype_depth_file)

df_mean_depth_tumor = calculate_mean_depth(df_tumor)
df_mean_depth_wildtype = calculate_mean_depth(df_wildtype)

df_log2_ratio = calculate_log2_ratio(df_mean_depth_tumor, df_mean_depth_wildtype)

save_log2_ratio_plot(df_log2_ratio)
