tar --use-compress-program="zstd -19 -T10" -cf real_data_v2.tar.zst \
  --exclude='*/tools_results.tsv' \
  results/real_data/subsamples/fraction_*/raw_outputs \
  results/real_data/subsamples/fraction_*/subsampled_data \
  results/real_data/subsamples/fraction_*/multi_metric_distances.parquet \
  results/real_data/subsamples/fraction_*/performance_results_hamming_mm5.tsv \
  results/real_data/subsamples/fraction_*/logs \
  results/real_data/subsamples/fraction_*/slurm_logs \
  results/real_data/subsamples/fraction_*/bash_scripts \
  results/real_data/subsamples/fraction_*/job_scripts


tar --use-compress-program="zstd -19 -T8" -cf simulated_data_v2.tar.zst \
  --exclude='*/tools_results.tsv' \
  results/simulated/ns_*/raw_outputs \
  results/simulated/ns_*/simulated_data \
  results/simulated/ns_*/multi_metric_distances.parquet \
  results/simulated/ns_*/performance_results_hamming_mm5.tsv \
  results/simulated/ns_*/logs \
  results/simulated/ns_*/slurm_logs \
  results/simulated/ns_*/bash_scripts \
  results/simulated/ns_*/job_scripts