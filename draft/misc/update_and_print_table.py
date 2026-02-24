import sys
from pathlib import Path
import polars as pl
import os

# Add src to python path
current_dir = Path.cwd()
src_path = current_dir / "src"
sys.path.append(str(src_path))

from bench.utils.functions import get_seq_stats_pl

sim_dir = current_dir / "results/simulated/ns_500_nc_5000_HIGH_INSERTION_RATE/simulated_data"
contigs_file = sim_dir / "simulated_contigs.fa"
spacers_file = sim_dir / "simulated_spacers.fa"

if not contigs_file.exists():
    print(f"Error: {contigs_file} not found")
    sys.exit(1)

if not spacers_file.exists():
    print(f"Error: {spacers_file} not found")
    sys.exit(1)

# specific column order for CSV
cols = ["file", "n_seqs", "median_length", "avg_length", "sum_length", "min_length", "max_length", "gc_frac_mean", "type", "name"]

# existing CSV
csv_path = current_dir / "draft/supp/figures/tables/all_dataset_sequence_stats.csv"

try:
    existing_df = pl.read_csv(csv_path, has_header=True)
    # Strip whitespace from column names if necessary (Polars usually handles this but good to be safe)
    existing_df.columns = [c.strip() for c in existing_df.columns]
except Exception as e:
    print(f"Error reading CSV: {e}")
    sys.exit(1)

# Check if we already have these specific entries to avoid duplicates?
# Use 'name' column. Expecting 'ns_500_nc_5000_HIGH_INSERTION_RATE'
run_name = "ns_500_nc_5000_HIGH_INSERTION_RATE"

# Compute stats
print("Computing stats for spacers...")
spacers_stats = get_seq_stats_pl(str(spacers_file))
spacers_stats = spacers_stats.with_columns(
    pl.lit("spacers").alias("type"),
    pl.lit(run_name).alias("name")
)

print("Computing stats for contigs...")
contigs_stats = get_seq_stats_pl(str(contigs_file))
contigs_stats = contigs_stats.with_columns(
    pl.lit("contigs").alias("type"),
    pl.lit(run_name).alias("name")
)

# Align columns
spacers_stats = spacers_stats.select(cols)
contigs_stats = contigs_stats.select(cols)

# Combine
new_rows = pl.concat([spacers_stats, contigs_stats])

# Cast to match commonly inferred types or just safely cast to String for CSV storage if mixed
# But meaningful types are better.
# Let's inspect schema of existing_df
print("Existing schema:", existing_df.schema)
print("New rows schema:", new_rows.schema)

# Force cast new_rows to match existing_df for critical columns if possible
# Or cast both to strict types
# Numeric columns: n_seqs, median_length, avg_length, sum_length, min_length, max_length, gc_frac_mean
num_cols = ["n_seqs", "median_length", "avg_length", "sum_length", "min_length", "max_length", "gc_frac_mean"]

for col in num_cols:
    existing_df = existing_df.with_columns(pl.col(col).str.strip_chars().cast(pl.Float64, strict=False))
    new_rows = new_rows.with_columns(pl.col(col).cast(pl.Float64, strict=False))

# String columns
str_cols = ["file", "type", "name"]
for col in str_cols:
    existing_df = existing_df.with_columns(pl.col(col).str.strip_chars().cast(pl.Utf8))
    new_rows = new_rows.with_columns(pl.col(col).str.strip_chars().cast(pl.Utf8))

# Filter out if already exists
existing_df = existing_df.filter(pl.col("name") != run_name)

# Concatenate
final_df = pl.concat([existing_df, new_rows])

# Write back
final_df.write_csv(csv_path)
print("Updated CSV successfully.")

# Now print the detailed markdown table for user
def format_num(n):
    try:
        f = float(n)
        if f.is_integer():
            return "{:,}".format(int(f))
        return "{:,.2f}".format(f)
    except:
        return str(n)

def format_sci(n):
    try:
        f = float(n)
        return "{:.1f}%".format(f * 100)
    except:
        return str(n)

print("\n| Dataset / Subset | N Seqs | Total Size (bp) | Length Range (bp) | Mean Length | GC % |")
print("|:---|---:|---:|:---|---:|---:|")

# We can reuse the same DF, iterating over rows
for row in final_df.iter_rows(named=True):
    name = row['name'].replace('_', ' ')
    
    # Custom display name logic if needed
    if row['type'] == 'spacers' and row['name'] != 'all_spacers':
         # Avoid duplicate printing if we want to group? 
         # But the requested table format in previous turn just listed everything.
         pass

    n_seqs = format_num(row['n_seqs'])
    sum_len = format_num(row['sum_length'])
    l_range = f"{format_num(row['min_length'])}-{format_num(row['max_length'])}"
    mean_len = format_num(row['avg_length'])
    gc = format_sci(row['gc_frac_mean'])
    
    # Add type to name for clarity?
    display_name = name
    if row['type'] == 'spacers' and 'spacers' not in name:
         display_name += " spacers"
    if row['type'] == 'contigs' and 'contigs' not in name: # "contigs" is usually not in the name like "fraction_0.1"
         pass # Actually look at previous table: "fraction 0.005", "ns 10000 nc 1000".
         # The 'type' column helps modify the display name if needed.
    
    # Check previous output:
    # | ns 100000 nc 10000 | 100,000 | 3,201,222 | 25-40 | 32.01 | 49.0% | -> This was spacers
    # | ns 100000 nc 10000 | 10,000 | ... -> This was contigs
    # It would be better to distinguish them.
    
    if row['type'] == 'spacers':
        full_name = f"{display_name} (spacers)"
    else:
        full_name = f"{display_name} (contigs)"
        
    if row['name'] == 'all_spacers' and row['type'] == 'spacers':
         full_name = "Real Spacers (iPHoP filtered)"
    elif 'fraction' in row['name']:
         full_name = f"Real Data {display_name}"

    print(f"| {full_name} | {n_seqs} | {sum_len} | {l_range} | {mean_len} | {gc} |")
