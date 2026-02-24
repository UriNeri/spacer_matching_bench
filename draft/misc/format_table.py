import csv
import sys

from pathlib import Path


#get file from 1st positional argument

input_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('draft/supp/figures/tables/all_dataset_sequence_stats.csv')

def format_num(n):
    try:
        f = float(n)
        if f.is_integer():
            return "{:,.0f}".format(f)
        return "{:,.2f}".format(f)
    except:
        return n


def format_sci(n):
    try:
        f = float(n)
        return "{:.1f}%".format(f * 100)
    except:
        return n

with open(input_file, 'r') as f:
    reader = csv.DictReader(f, skipinitialspace=True)
    reader.fieldnames = [name.strip() for name in reader.fieldnames]
    rows = list(reader)
has_type = 'type' in reader.fieldnames

if input_file.name == 'all_dataset_sequence_stats.csv' or has_type:
    if has_type:
        print("| Dataset / Subset | Type | N Seqs | Total Size (bp) | Length Range (bp) | Mean Length | GC % |")
        print("|:---|:---|---:|---:|:---|---:|---:|")
    else:
        print("| Dataset / Subset | N Seqs | Total Size (bp) | Length Range (bp) | Mean Length | GC % |")
        print("|:---|---:|---:|:---|---:|---:|")
else:
    char_counts = {col: sum(len(row[col]) for row in rows) for col in reader.fieldnames}
    header_str = " |: ".join(reader.fieldnames)
    width_str = " | ".join(f"{{.{char_counts[col]}}}" for col in reader.fieldnames)
    print(f"| {header_str} |")
    print(f"| {width_str} |")


for row in rows:
    name = row['name'].replace('_', ' ')
    n_seqs = format_num(row['n_seqs'])
    sum_len = format_num(row['sum_length'])
    l_range = f"{format_num(row['min_length'])}-{format_num(row['max_length'])}"
    mean_len = format_num(row['avg_length'])
    gc = format_sci(row['gc_frac_mean'])
    if has_type:
        type_label = row['type'].strip().capitalize()
        print(f"| {name} | {type_label} | {n_seqs} | {sum_len} | {l_range} | {mean_len} | {gc} |")
    else:
        print(f"| {name} | {n_seqs} | {sum_len} | {l_range} | {mean_len} | {gc} |")
