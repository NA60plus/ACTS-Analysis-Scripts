import pandas as pd
import os
import os
import shutil
import glob
XSizeTot = 135.996
YSizeTot = 136.948 # readout side
DeadXRight = 1.5  # right end cap  (readout)
DeadXLeft = 4.5   # left end cap
DeadYBottom = 0.525 # bottom dead zone
DeadYTop = 0.525 # top dead zone
DeadTopBotHalves = 0.12 # dead space between top and bottom halves
DeadXTile = 0.02 # dead space between tiles in X (apart from DeadXDataBackbone after each 3 tiles)
DeadXDataBackbone = 0.06 # dead space between segments
NXTiles = 36
NXSegments = 12 # group of 3 tiles
NYSensors = 7
DYSens = (YSizeTot - DeadYBottom - DeadYTop) / NYSensors
DXSegment = (XSizeTot - DeadXRight - DeadXLeft) / NXSegments
DXTile = (DXSegment - DeadXDataBackbone) / 3


offsX = [2.9, -3.1,  -2.9 ,  3.1]
offsY = [3.1,  2.9,  -3.1 , -2.9]

def isInAcc(x, y, z):
  if z > 400:
     return True
  if x < offsX[0] and y > offsY[0]:
    geoId = 0
  elif x < offsX[1] and y < offsY[1]:
    geoId = 1
  elif x > offsX[2] and y < offsY[2]:
    geoId = 2
  elif x > offsX[3] and y > offsY[3]:
    geoId = 3
  else:
    return False
  x -= offsX[geoId]
  y -= offsY[geoId]

  if (x < 0.):
    x += XSizeTot # relate to chip local left/bottom corner
  else:
    x = XSizeTot - x # relate to chip local left/bottom corner

  if y<0:
    y += YSizeTot
  else:
    y = YSizeTot - y
  
  if x < DeadXLeft or x > XSizeTot - DeadXRight or y < DeadYBottom or y > YSizeTot - DeadYTop:
    return False
  
  x -= DeadXLeft
  y -= DeadYBottom

  sens = int(y/DYSens)
  y -= sens*DYSens
  if y < DeadYBottom or y > DYSens - DeadYTop :
    return False
  segm = int(x/DXSegment)
  x -= segm*DXSegment
  if x > DXSegment - DeadXDataBackbone:
    return False
  tile = int(x/DXTile)
  x -= tile*DXTile
  if x<DeadXTile:
    return False

  return True

# Define bitmask and shift details
kLayerMask = 0x0000FFF000000000  # Example: Adjust based on actual bit positions
layer_shift = 10

def shift_layer(geometry_id):
    """
    Extracts the layer bits, shifts them by 10, and reconstructs the new geometry_id.
    
    :param geometry_id: The original encoded geometry identifier.
    :return: The modified geometry identifier.
    """
    layer_value = (geometry_id & kLayerMask) >> 36  # Extract layer
    new_layer_value = (layer_value + layer_shift) & 0xFFF  # Ensure it fits in 36 bits
    new_geometry_id = (geometry_id & ~kLayerMask) | (new_layer_value << 36)  # Replace layer
    return new_geometry_id




def update_geometry_ids(file_path, output_file):
    """
    Reads a CSV file, shifts the layer in 'geometry_id' by 10, and saves the modified file.

    :param file_path: Path to the input CSV file.
    :param output_file: Path to save the modified CSV file.
    """
    df = pd.read_csv(file_path)
    # Apply layer shift to geometry_id
    df["geometry_id"] = df["geometry_id"].apply(shift_layer)

    # Save the modified CSV
    df.to_csv(output_file, index=False)

    return df

import pandas as pd
import numpy as np

def filter_dataframe(df, keep_probability, output_file):
    """
    Randomly filters rows from a DataFrame based on a given probability.

    Args:
        df (pd.DataFrame): The input DataFrame.
        keep_probability (float): Value between 0 and 1. Each row has this chance to be kept.

    Returns:
        pd.DataFrame: A new DataFrame with only the kept rows.
    """
    if not (0 <= keep_probability <= 1):
        raise ValueError("keep_probability must be between 0 and 1")

    # Generate a random number for each row
    random_values = np.random.rand(len(df))

    # Keep rows where the random number is less than or equal to the probability
    mask = random_values <= keep_probability
    df = df[mask].reset_index(drop=True)
    # Save the modified CSV
    df.to_csv(output_file, index=False)

    return df

def filter_csv(df, output_file):
    
    # Filter rows where the 'vx' column is not 0
    df_filtered = df[df['particle_id'] != 0]
    
    # Write the filtered DataFrame to a new CSV file
    df_filtered.to_csv(output_file, index=False)

    return df_filtered

def apply_dead_zones(df, output_file):
    
    # Filter rows where the 'vx' column is not 0
    df_filtered = df[df.apply(lambda row: isInAcc(row["tx"], row["ty"], row["tz"]), axis=1)]
    
    # Write the filtered DataFrame to a new CSV file
    df_filtered.to_csv(output_file, index=False)

def filter_by_particle_id(input_file, reference_file, output_file):
    # Read the CSV files into DataFrames
    df = pd.read_csv(input_file)
    ref_df = pd.read_csv(reference_file)
    
    # Get the 'particle_id' values from the reference file
    reference_ids = ref_df['particle_id'].unique()
    
    # Filter rows where 'particle_id' is present in the reference file
    df_filtered_by_id = df[df['particle_id'].isin(reference_ids)]
    
    # Write the filtered DataFrame to a new CSV file
    df_filtered_by_id.to_csv(output_file, index=False)

def filter_by_binary_prefix(input_file, output_file, n, bit_width=64):
    """
    Reads a CSV file and filters rows where the binary representation (with fixed width)
    of the 'particle_id' starts with n zeros.
    
    Parameters:
      input_file: Path to the input CSV.
      output_file: Path for the filtered output CSV.
      n: Number of leading zeros required.
      bit_width: The bit width for the binary representation (default is 64).
    """
    df = pd.read_csv(input_file)
    
    def has_leading_zeros(val):
        # Convert particle_id to integer (in case it's read as float)
        particle_int = int(val)
        # Format as a binary string with fixed bit width (padded with zeros)
        binary_str = format(particle_int, '0{}b'.format(bit_width))
        return binary_str.startswith('0' * n)
    
    # Apply the filtering function on the 'particle_id' column
    df_filtered = df[df['particle_id'].apply(has_leading_zeros)]
    df_filtered.to_csv(output_file, index=False)

def copy_and_rename_files(src_dir, dest_dir):
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    
    for filename in os.listdir(src_dir):
        if filename.endswith("primaries.csv"):
            new_filename = filename.replace("primaries", "particles")
            src_path = os.path.join(src_dir, filename)
            dest_path = os.path.join(dest_dir, new_filename)
            
            shutil.copy(src_path, dest_path)
            print(f"Copied and renamed: {filename} -> {new_filename}")


# Directory containing the CSV files
input_directory_list = ["simulatedEvents/rubenxprino_40GeV_omega_k0s_lambda"
                        ]

eff = 0.99
for input_directory in input_directory_list: # Change this to your input directory path
  output_directory = input_directory+"_filtered"
  output_directory_eff = input_directory+"_filtered_eff"+str(eff)
  output_directory_dead = input_directory+"_deadzones"
  # Create the output directory if it doesn't exist
  os.makedirs(output_directory, exist_ok=True)
  os.makedirs(output_directory_eff, exist_ok=True)
  os.makedirs(output_directory_dead, exist_ok=True)
  # Get all CSV files in the directory ending with "_hits.csv"
  csv_files_hits = glob.glob(os.path.join(input_directory, "*-hits.csv"))
  csv_files_particles = glob.glob(os.path.join(input_directory, "*-particles.csv"))
  # Process each file
  for hits,particles in zip(csv_files_hits, csv_files_particles):

      # Read the CSV file into a DataFrame
      df_hits = pd.read_csv(hits)
      df_particles = pd.read_csv(particles)

      output_file = output_directory+hits.replace(input_directory,"")
      df_hits_filtered = filter_csv(df_hits, output_file)
      #updates the geo id and saves it in output_directory
      df_shifted = update_geometry_ids(output_file, output_file)

      output_file = output_directory_eff+hits.replace(input_directory,"")
      _ = filter_dataframe(df_shifted, eff, output_file)

      output_file = output_directory_dead+hits.replace(input_directory,"")
      apply_dead_zones(df_hits_filtered, output_file)

      #copy the generated particles in the other directories
      output_file = output_directory+particles.replace(input_directory,"")
      shutil.copy(particles, output_file)

      output_file = output_directory_eff+particles.replace(input_directory,"")
      shutil.copy(particles, output_file)

      output_file = output_directory_dead+particles.replace(input_directory,"")
      shutil.copy(particles, output_file)



  output_directory = input_directory+"_primaries"
  os.makedirs(output_directory, exist_ok=True)

  copy_and_rename_files(input_directory, output_directory)


