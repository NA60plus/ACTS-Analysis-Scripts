import os
import shutil
import argparse

def copy_files_in_directory(src_dir, n):
    """
    Copies each file in the given directory `n` times.
    
    Parameters:
    src_dir (str): Source directory path.
    n (int): Number of copies to create for each file.
    """
    try:
        # List all files in the source directory
        files = [f for f in os.listdir(src_dir) if os.path.isfile(os.path.join(src_dir, f))]
        
        number_of_files = len(files)
        # Loop through each file in the directory
        for file in files:
            src_file_path = os.path.join(src_dir, file)
            num_str = file[5:14]  # "000000002" in "event000000002-hits.csv"
            num = int(num_str)
            # Create `n` copies of each file
            for i in range(1, n + 1):
                new_num_str = f"{(num+number_of_files*i):09d}"  # Zero-padded to 9 digits
                # Construct the new filename
                new_filename = f"event{new_num_str}-hits.csv"
                # Define destination file path with an appended copy number
                dest_file_path = os.path.join(src_dir, new_filename)
                
                # Copy the file
                shutil.copy2(src_file_path, dest_file_path)
                print(f"Created copy {i} for file {file}")
    
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", dest="dir", type=str, help="Source directory containing the files to rename.", default="/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/event_generation/simulatedEvents/bkghits_120GeV_vt/")
    parser.add_argument("-n", dest="ncopies", type=int, help="Value to add to the number in each filename.", default=1)
    
    # Parse arguments
    args = parser.parse_args()

    copy_files_in_directory(args.dir, args.ncopies)
