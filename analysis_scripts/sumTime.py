import csv

# Path to the .tsv file
ts_file_path = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_omega_full_roberta/timing.tsv"

# Initialize sums
time_total_sum = 0.0
time_perevent_sum = 0.0

# Read the .tsv file and sum the values
with open(ts_file_path, newline='') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    next(reader)  # Skip header
    for row in reader:
        time_total_sum += float(row[1])
        time_perevent_sum += float(row[2])

# Print results
print(f"Total time_total_s: {time_total_sum:.6f} s")
print(f"Total time_perevent_s: {time_perevent_sum:.6f} s")