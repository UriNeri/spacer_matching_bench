"""
Results comparison and validation command (compares the tool outputs against the ground truth).
"""

import os
import logging
from pathlib import Path
from typing import Dict, Optional
import polars as pl
# import polars_bio as pb
from rich.console import Console

from bench.utils.functions import (
    read_fasta,
    read_results,
    test_alignment,
    populate_pldf_withseqs_needletail,
    prettify_alignment_edit,
    prettify_alignment_gap_affine,
    prettify_alignment_hamming,
    get_seq_from_fastx,
    vstack_easy,
    match_intervals_with_tolerance,
    calculate_hamming_distance,
    calculate_edit_distance,
    calculate_gap_affine_edit,

)
from bench.utils.tool_commands import load_tool_configs

# Logger is configured by cli.py with RichHandler
logger = logging.getLogger(__name__)
console = Console(soft_wrap=False)

# Polars display config (keep wide to avoid truncation in debug prints)
pl.Config.set_tbl_rows(123123)
pl.Config.set_tbl_cols(123123)
pl.Config.set_fmt_str_lengths(2100)
pl.Config.set_tbl_width_chars(2100)


def merge_overlapping_intervals(df: pl.DataFrame, tolerance: int = 0) -> pl.DataFrame:
    """Merge overlapping intervals per spacer/contig/strand with an optional tolerance."""
    required_cols = ["spacer_id", "contig_id", "start", "end", "strand"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns for interval merge: {missing}")

    if df.height == 0:
        return df.select(required_cols)

    base = df.select(required_cols).sort(required_cols)

    def _merge_list(intervals: list[dict]) -> list[dict]:
        merged: list[dict] = []
        current_start: Optional[int] = None
        current_end: Optional[int] = None
        for item in intervals:
            s = int(item["start"])
            e = int(item["end"])
            if current_start is None:
                current_start, current_end = s, e
                continue
            if s <= current_end + tolerance:
                current_end = max(current_end, e)
            else:
                merged.append({"start": current_start, "end": current_end})
                current_start, current_end = s, e
        if current_start is not None:
            merged.append({"start": current_start, "end": current_end})
        return merged

    merged = (
        base.group_by(["spacer_id", "contig_id", "strand"], maintain_order=True)
        .agg(pl.struct(["start", "end"]).alias("intervals"))
        .with_columns(
            pl.col("intervals").map_elements(
                _merge_list,
                return_dtype=pl.List(pl.Struct({"start": pl.Int64, "end": pl.Int64})),
            )
        )
        .explode("intervals")
        .unnest("intervals")
        .select(required_cols)
    )

    return merged


def classify_unique_alignments_across_tools(
    ground_truth: pl.DataFrame,
    all_tool_results: pl.DataFrame,
    contigs_file: str,
    spacers_file: str,
    max_distance: int,
    distance_metric: str,
    gap_open_penalty: int,
    gap_extend_penalty: int,
    coordinate_tolerance: int,
    multi_metric: bool = False,
) -> pl.DataFrame:
    """Classify unique alignments across all tools with distance recalculation and GT matching.
    
    Args:
        multi_metric: If True, don't filter by distance here (let performance calc do it per metric)
    """

    if all_tool_results.height == 0:
        return pl.DataFrame(
            {
                "alignment_idx": [],
                "region_idx": [],
                "spacer_id": [],
                "contig_id": [],
                "start": [],
                "end": [],
                "strand": [],
                "classification": [],
                "start_gt": [],
                "end_gt": [],
                "mismatches_gt": [],
                "planned_mismatches": [],
                "distance_hamming": [],
                "distance_edit": [],
                "distance_gap_affine": [],
                "recalculated_distance": [],
            },
            schema={
                "alignment_idx": pl.Int64,
                "region_idx": pl.Int64,
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "start": pl.Int64,
                "end": pl.Int64,
                "strand": pl.Boolean,
                "classification": pl.Utf8,
                "start_gt": pl.Int64,
                "end_gt": pl.Int64,
                "mismatches_gt": pl.Int64,
                "planned_mismatches": pl.Int64,
                "distance_hamming": pl.Int64,
                "distance_edit": pl.Int64,
                "distance_gap_affine": pl.Int64,
                "recalculated_distance": pl.Int64,
            },
        )

    unique_alignments = (
        all_tool_results.select(["spacer_id", "contig_id", "start", "end", "strand"])
        .unique()
        .with_row_index("alignment_idx")
    )

    merged_regions = merge_overlapping_intervals(
        unique_alignments, tolerance=coordinate_tolerance
    ).with_row_index("region_idx")

    unique_alignments = (
        unique_alignments.join(
            merged_regions,
            on=["spacer_id", "contig_id", "strand"],
            how="left",
            suffix="_region",
        )
        .filter(
            (pl.col("start") >= pl.col("start_region") - coordinate_tolerance)
            & (pl.col("start") <= pl.col("end_region") + coordinate_tolerance)
            & (pl.col("end") >= pl.col("start_region") - coordinate_tolerance)
            & (pl.col("end") <= pl.col("end_region") + coordinate_tolerance)
        )
        .select(["spacer_id", "contig_id", "start", "end", "strand", "alignment_idx", "region_idx"])
    )

    matches = match_intervals_with_tolerance(
        ground_truth=ground_truth,
        tool_results=unique_alignments,
        tolerance=coordinate_tolerance,
    )

    alignments = unique_alignments.join(
        matches,
        on=["spacer_id", "contig_id", "start", "end", "strand"],
        how="left",
        suffix="_match",
    )

    alignments = alignments.join(
        ground_truth.rename({"start": "start_gt", "end": "end_gt"}).select(
            ["spacer_id", "contig_id", "start_gt", "end_gt", "strand", "mismatches"]
        ),
        on=["spacer_id", "contig_id", "start_gt", "end_gt", "strand"],
        how="left",
    ).rename({"mismatches": "mismatches_gt"})

    alignments = alignments.with_columns(
        pl.col("classification").fill_null("needs_verification"),
        pl.col("start_gt").cast(pl.Int64),
        pl.col("end_gt").cast(pl.Int64),
    )

    alignments = populate_pldf_withseqs_needletail(
        alignments,
        seqfile=spacers_file,
        idcol="spacer_id",
        seqcol="spacer_seq",
        trim_to_region=False,
        reverse_by_strand_col=False,
    )

    alignments = populate_pldf_withseqs_needletail(
        alignments,
        seqfile=contigs_file,
        idcol="contig_id",
        seqcol="contig_seq",
        trim_to_region=True,
        reverse_by_strand_col=True,
    )

    alignments = alignments.with_columns(
        pl.struct([
            pl.col("spacer_seq"),
            pl.col("contig_seq"),
            pl.col("strand"),
        ]).map_elements(
            lambda x: {
                "distance_hamming": calculate_hamming_distance(
                    spacer_seq=x["spacer_seq"],
                    contig_seq=x["contig_seq"],
                    strand=x["strand"],
                ),
                "distance_edit": calculate_edit_distance(
                    spacer_seq=x["spacer_seq"],
                    contig_seq=x["contig_seq"],
                    strand=x["strand"],
                ),
                "distance_gap_affine": calculate_gap_affine_edit(
                    spacer_seq=x["spacer_seq"],
                    contig_seq=x["contig_seq"],
                    strand=x["strand"],
                    gap_open=gap_open_penalty,
                    gap_extend=gap_extend_penalty,
                ),
            },
            return_dtype=pl.Struct(
                {
                    "distance_hamming": pl.Int64,
                    "distance_edit": pl.Int64,
                    "distance_gap_affine": pl.Int64,
                }
            ),
        ).alias("distances")
    ).unnest("distances")

    distance_expr = {
        "edit": pl.col("distance_edit"),
        "gap_affine": pl.col("distance_gap_affine"),
    }.get(distance_metric, pl.col("distance_hamming"))

    # In multi-metric mode, mark as invalid only if ALL three metrics exceed threshold
    if multi_metric:
        alignments = alignments.with_columns(
            distance_expr.alias("recalculated_distance"),
            # For alignments not in ground truth, check if ALL metrics exceed threshold
            pl.when(pl.col("classification") == "needs_verification")
            .then(
                pl.when(
                    (pl.col("distance_hamming") > max_distance) &
                    (pl.col("distance_edit") > max_distance) &
                    (pl.col("distance_gap_affine") > max_distance)
                )
                .then(pl.lit("invalid_alignment"))
                .otherwise(pl.lit("positive_not_in_plan"))
            )
            .otherwise(pl.col("classification"))
            .alias("classification"),
            pl.when(pl.col("mismatches_gt").is_not_null())
            .then(pl.col("mismatches_gt"))
            .otherwise(pl.col("distance_hamming"))
            .alias("planned_mismatches"),
        )
    else:
        # Single-metric mode: filter by distance
        alignments = alignments.with_columns(
            distance_expr.alias("recalculated_distance"),
            pl.when(pl.col("classification") == "positive_in_plan")
            .then(
                pl.when(distance_expr <= max_distance)
                .then(pl.lit("positive_in_plan"))
                .otherwise(pl.lit("invalid_alignment"))
            )
            .otherwise(
                pl.when(distance_expr <= max_distance)
                .then(pl.lit("positive_not_in_plan"))
                .otherwise(pl.lit("invalid_alignment"))
            )
            .alias("classification"),
            pl.when(pl.col("mismatches_gt").is_not_null())
            .then(pl.col("mismatches_gt"))
            .otherwise(pl.col("distance_hamming"))
            .alias("planned_mismatches"),
        )

    # Drop sequences before returning to keep outputs compact
    return alignments.select(
        [
            "alignment_idx",
            "region_idx",
            "spacer_id",
            "contig_id",
            "start",
            "end",
            "strand",
            "classification",
            "start_gt",
            "end_gt",
            "mismatches_gt",
            "planned_mismatches",
            "distance_hamming",
            "distance_edit",
            "distance_gap_affine",
            "recalculated_distance",
        ]
    )


def unified_multi_metric_pipeline(
    tools_results: pl.DataFrame,
    ground_truth: pl.DataFrame,
    contigs_file: str,
    spacers_file: str,
    max_distance: int,
    distance_metric: str,
    coordinate_tolerance: int,
    gap_open_penalty: int,
    gap_extend_penalty: int,
    multi_metric: bool = False,
):
    """End-to-end unified multi-metric processing and performance calculation.
    
    Args:
        multi_metric: If True, calculate performance for both hamming and edit distance.
                     In multi-metric mode, alignments are marked invalid only if ALL
                     three metrics (hamming, edit, gap_affine) exceed max_distance.
    
    Returns:
        Tuple of (tools_results_out, tools_results_with_gt_out, region_results, performance_results)
        - region_results includes a 'tools' column (List[str]) with tool names that found each alignment
        - tools column excludes 'ground_truth' but includes 'all_tools_combined'
    """

    # Add ground truth as a "tool" for unified processing
    if ground_truth.height > 0:
        # Select only the common columns to avoid schema mismatch
        common_cols = ["spacer_id", "contig_id", "start", "end", "strand"]
        gt_as_tool = ground_truth.select(common_cols).with_columns(
            pl.lit("ground_truth").alias("tool")
        )
        tools_results_with_gt = vstack_easy(
            tools_results.select([*common_cols, "tool"]),
            gt_as_tool
        )
        # Rejoin with the full tools_results to get all columns back
        tools_results_with_gt = tools_results_with_gt.join(
            tools_results,
            on=["spacer_id", "contig_id", "start", "end", "strand", "tool"],
            how="left",
        )
    else:
        tools_results_with_gt = tools_results

    # Add "all_tools_combined" synthetic tool for tool-independent analysis
    # This combines all unique alignments from all real tools (excludes ground_truth)
    logger.debug("Creating 'all_tools_combined' synthetic tool...")
    real_tools_only = tools_results_with_gt.filter(pl.col("tool") != "ground_truth")
    
    if real_tools_only.height > 0:
        # Get unique alignments by coordinates (deduplicate across all tools)
        all_tools_combined = (
            real_tools_only
            .select([col for col in real_tools_only.columns if col != "tool"])
            .unique(subset=["spacer_id", "contig_id", "start", "end", "strand"])
            .with_columns(pl.lit("all_tools_combined").alias("tool"))
        )
        
        # Add to tools_results_with_gt for unified processing
        tools_results_with_gt = vstack_easy(tools_results_with_gt, all_tools_combined)
        logger.debug(f"Added 'all_tools_combined' with {all_tools_combined.height} unique alignments")
    else:
        logger.debug("No real tool results found - skipping 'all_tools_combined'")

    alignment_classifications = classify_unique_alignments_across_tools(
        ground_truth=ground_truth,
        all_tool_results=tools_results_with_gt,
        contigs_file=contigs_file,
        spacers_file=spacers_file,
        max_distance=max_distance,
        distance_metric=distance_metric,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty,
        coordinate_tolerance=coordinate_tolerance,
        multi_metric=multi_metric,
    )

    # Join classifications to tools_results (excluding ground_truth tool, but keeping all_tools_combined)
    # Filter out ground_truth from tools_results_with_gt to get real tools + all_tools_combined
    tools_results_excluding_gt = tools_results_with_gt.filter(pl.col("tool") != "ground_truth")
    tools_results_out = tools_results_excluding_gt.join(
        alignment_classifications,
        on=["spacer_id", "contig_id", "start", "end", "strand"],
        how="left",
    )

    # Also create version with ground_truth for display purposes
    tools_results_with_gt_out = tools_results_with_gt.join(
        alignment_classifications,
        on=["spacer_id", "contig_id", "start", "end", "strand"],
        how="left",
    )

    performance_results = calculate_all_tool_performance(
        tools_results=tools_results_out,
        alignment_classifications=alignment_classifications,
        ground_truth=ground_truth,
        max_distance=max_distance,
        distance_metric=distance_metric,
    )
    
    # If multi-metric mode, calculate performance for both hamming and edit
    if multi_metric:
        performance_hamming = calculate_all_tool_performance(
            tools_results=tools_results_out,
            alignment_classifications=alignment_classifications,
            ground_truth=ground_truth,
            max_distance=max_distance,
            distance_metric="hamming",
        )
        performance_edit = calculate_all_tool_performance(
            tools_results=tools_results_out,
            alignment_classifications=alignment_classifications,
            ground_truth=ground_truth,
            max_distance=max_distance,
            distance_metric="edit",
        )
        
        # Rename columns to add metric suffixes (except 'tool' and 'ground_truth_planned')
        rename_hamming = {col: f"{col}_hamming" for col in performance_hamming.columns if col not in ["tool", "ground_truth_planned"]}
        rename_edit = {col: f"{col}_edit" for col in performance_edit.columns if col not in ["tool", "ground_truth_planned"]}
        
        performance_hamming = performance_hamming.rename(rename_hamming)
        performance_edit = performance_edit.rename(rename_edit)
        
        # Join on tool  (also join on ground_truth_planned since it's the same)
        performance_results = performance_hamming.join(
            performance_edit,
            on=["tool", "ground_truth_planned"],
            how="outer",
        )

    # Add tools column to region_results for notebook compatibility
    # Group tools by alignment coordinates (exclude ground_truth, include all_tools_combined)
    tool_assignments = (
        tools_results_with_gt
        .filter(pl.col("tool") != "ground_truth")  # Exclude ground_truth from tools list
        .group_by(["spacer_id", "contig_id", "start", "end", "strand"])
        .agg(pl.col("tool").unique().sort().alias("tools"))
    )
    
    # Join tools to alignment classifications
    region_results = alignment_classifications.join(
        tool_assignments,
        on=["spacer_id", "contig_id", "start", "end", "strand"],
        how="left",
    )
    
    # For alignments from ground_truth only (not found by any real tool),
    # set tools to empty list
    region_results = region_results.with_columns(
        pl.col("tools").fill_null([])
    )

    return tools_results_out, tools_results_with_gt_out, region_results, performance_results

def run_compare_results(
    input_dir,
    max_mismatches=5,
    output_file=None,
    threads=4,
    skip_tools="",
    only_tools=None,
    contigs=None,
    spacers=None,
    distance_metric="hamming",
    gap_open_penalty=5,
    gap_extend_penalty=5,
    logfile=None,
    skip_hyperfine=False,
    tools_results_out_file=None,
    coordinate_tolerance: int = 5,
    multimetric_output: Optional[str] = None,
    ground_truth_file: Optional[str] = None,
    verbose=False,
    multi_metric=False,
):
    """
    Unified multi-metric compare-results entrypoint. Works for simulated and real datasets.
    Always runs the unified pipeline: region grouping with tolerance, per-alignment distances,
    region-level classifications, and per-tool performance summary.
    
    The performance results will include an 'all_tools_combined' synthetic tool entry that
    represents the union of all real tool results (deduplicated by coordinates). This is
    useful for tool-independent analysis in notebooks.
    
    Example CLI usage:
        # Single metric (hamming distance)
        pixi run spacer_bencher compare-results -i results/simulated/run1 -mm 5 --distance hamming
        
        # Multi-metric mode (both hamming and edit distance)
        pixi run spacer_bencher compare-results -i results/simulated/run1 -mm 5 --multi-metric
        
        # Real data with custom files
        pixi run spacer_bencher compare-results -i results/real_data/subsamples/fraction_0.001 \\
            --contigs results/real_data/subsamples/fraction_0.001/subsampled_data/subsampled_contigs.fa \\
            --spacers ./imgvr4_data/spacers/iphop_filtered_spacers.fna \\
            -mm 5 --distance edit --skip-hyperfine --multi-metric
    
    # Example commented-out code block for notebooks:
    # from bench.commands.compare_results import run_compare_results
    # tools_results, region_results, performance = run_compare_results(
    #     input_dir="results/simulated/run1",
    #     max_mismatches=5,
    #     multi_metric=True,  # Calculate both hamming and edit metrics
    #     distance_metric="hamming",  # Default metric for single-metric mode
    # )
    # # Filter for all_tools_combined synthetic tool
    # combined_perf = performance.filter(pl.col("tool") == "all_tools_combined")
    """

    logger.debug(f"Processing tool results from {input_dir}")

    if not os.path.exists(input_dir):
        logger.error(f"Input directory {input_dir} does not exist")
        raise FileNotFoundError(f"Input directory {input_dir} does not exist")

    contigs_path = contigs if contigs else f"{input_dir}/simulated_data/simulated_contigs.fa"
    spacers_path = spacers if spacers else f"{input_dir}/simulated_data/simulated_spacers.fa"

    tools_results_output = tools_results_out_file or f"{input_dir}/tools_results.tsv"
    region_output_path = multimetric_output or f"{input_dir}/multi_metric_distances.parquet"
    perf_output = output_file or f"{input_dir}/performance_results_{distance_metric}_mm{max_mismatches}.tsv"

    logger.debug(f"Contigs: {contigs_path}")
    logger.debug(f"Spacers: {spacers_path}")

    logger.debug("Loading tool configurations...")
    tools = load_tool_configs(
        results_dir=input_dir,
        threads=threads,
        contigs_file=contigs,
        spacers_file=spacers,
    )
    logger.debug(f"Loaded {len(tools)} tool configurations")

    if skip_tools:
        skip_list = skip_tools.split(",")
        tools = {k: v for k, v in tools.items() if k not in skip_list}
        logger.debug(f"Remaining tools after skip: {len(tools)}")

    if only_tools:
        only_list = only_tools.split(",")
        tools = {k: v for k, v in tools.items() if k in only_list}
        logger.debug(f"Tools to process: {len(tools)}")

    logger.debug("Reading spacer sequences...")
    spacers_dict = read_fasta(spacers_path)
    spacer_lendf = pl.DataFrame(
        {
            "spacer_id": spacers_dict.keys(),
            "length": [len(seq) for seq in spacers_dict.values()],
        }
    )

    logger.debug("Reading tool alignment results...")
    tools_results = read_results(
        tools,
        max_mismatches=max_mismatches,
        spacer_lendf=spacer_lendf,
        ref_file=contigs_path,
    )
    logger.info(f"Read {tools_results.height} total alignment results")

    # Ground truth (optional). Empty for real datasets.
    if ground_truth_file is None:
        default_gt = f"{input_dir}/simulated_data/planned_ground_truth.tsv"
        alt_gt = f"{input_dir}/simulated_data/ground_truth.tsv"
        ground_truth_path = default_gt if os.path.exists(default_gt) else alt_gt
    else:
        ground_truth_path = ground_truth_file

    if os.path.exists(ground_truth_path):
        logger.debug(f"Reading ground truth from {ground_truth_path}")
        ground_truth = pl.read_csv(ground_truth_path, separator="\t")
        if ground_truth.height > 0:
            ground_truth = ground_truth.filter(pl.col("mismatches") <= max_mismatches)
            logger.info(
                f"{ground_truth.height} ground truth annotations within max mismatches of {max_mismatches}"
            )
        else:
            ground_truth = pl.DataFrame(
                {
                    "spacer_id": pl.Series(dtype=pl.Utf8),
                    "contig_id": pl.Series(dtype=pl.Utf8),
                    "start": pl.Series(dtype=pl.Int64),
                    "end": pl.Series(dtype=pl.Int64),
                    "strand": pl.Series(dtype=pl.Boolean),
                    "mismatches": pl.Series(dtype=pl.Int64),
                }
            )
            logger.info("Ground truth file empty; proceeding with empty GT")
    else:
        ground_truth = pl.DataFrame(
            {
                "spacer_id": pl.Series(dtype=pl.Utf8),
                "contig_id": pl.Series(dtype=pl.Utf8),
                "start": pl.Series(dtype=pl.Int64),
                "end": pl.Series(dtype=pl.Int64),
                "strand": pl.Series(dtype=pl.Boolean),
                "mismatches": pl.Series(dtype=pl.Int64),
            }
        )
        logger.info("No ground truth provided; proceeding with empty GT")

    tools_results_out, tools_results_with_gt_out, region_results, performance_results = unified_multi_metric_pipeline(
        tools_results=tools_results,
        ground_truth=ground_truth,
        contigs_file=contigs_path,
        spacers_file=spacers_path,
        max_distance=max_mismatches,
        distance_metric=distance_metric,
        coordinate_tolerance=coordinate_tolerance,
        gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty,
        multi_metric=multi_metric,
    )

    if logger.isEnabledFor(logging.DEBUG):
        try:
            display_example_alignments(
                alignment_classifications=region_results,
                tools_results=tools_results_with_gt_out,
                contigs_file=contigs_path,
                spacers_file=spacers_path,
                max_distance=max_mismatches,
                num_examples=3,
                distance_metric=distance_metric,
                gap_open_penalty=gap_open_penalty,
                gap_extend_penalty=gap_extend_penalty,
            )
        except Exception as e:
            logger.debug(f"Skipping example alignments display: {e}")

    tools_results_out.write_csv(tools_results_output, separator="\t")
    logger.info(f"Wrote tool results to {tools_results_output}")

    region_results.write_parquet(region_output_path)
    logger.info(f"Wrote region-level results to {region_output_path}")

    performance_results.write_csv(perf_output, separator="\t")
    logger.info(f"Wrote performance results to {perf_output}")

    # Optional stdout summary
    if performance_results.height > 0:
        if multi_metric:
            # In multi-metric mode, display all metrics for both hamming and edit
            summary_cols = [
                "tool",
                "ground_truth_planned",
                "ground_truth_augmented_hamming",
                "recall_planned_hamming",
                "recall_augmented_hamming",
                "all_true_positives_hamming",
                "planned_true_positives_hamming",
                "positives_not_in_plan_hamming",
                "invalid_alignments_hamming",
                "false_negatives_planned_hamming",
                "false_negatives_augmented_hamming",
                "ground_truth_augmented_edit",
                "recall_planned_edit",
                "recall_augmented_edit",
                "all_true_positives_edit",
                "planned_true_positives_edit",
                "positives_not_in_plan_edit",
                "invalid_alignments_edit",
                "false_negatives_planned_edit",
                "false_negatives_augmented_edit",
            ]
            sort_col = "recall_planned_hamming"
        else:
            summary_cols = [
                "tool",
                "ground_truth_planned",
                "ground_truth_augmented",
                "recall_planned",
                "recall_augmented",
                "all_true_positives",
                "planned_true_positives",
                "positives_not_in_plan",
                "invalid_alignments",
                "false_negatives_planned",
                "false_negatives_augmented",
            ]
            sort_col = "recall_planned"
        
        summary_table = performance_results.select([c for c in summary_cols if c in performance_results.columns]).sort(
            sort_col, descending=True
        )
        pl.Config.set_tbl_cols(-1)
        pl.Config.set_tbl_rows(15)
        # logger.info("")
        print(summary_table)
    else:
        logger.warning("No performance results to display")

    return tools_results_out, region_results, performance_results


def display_example_alignments(
    alignment_classifications: pl.DataFrame,
    tools_results: pl.DataFrame,
    contigs_file: str,
    spacers_file: str,
    max_distance: int = 5,
    num_examples: int = 3,
    distance_metric: str = "hamming",
    gap_open_penalty: int = 5,
    gap_extend_penalty: int = 5,
):
    """
    Display example alignments for each classification type.

    Args:
        alignment_classifications: DataFrame with classified alignments
        tools_results: DataFrame with tool results to show which tools found each alignment
        contigs_file: Path to contigs FASTA file
        spacers_file: Path to spacers FASTA file
        max_mismatches: Maximum allowed mismatches
        num_examples: Number of examples to show per classification
        distance_metric: Distance metric for the performance evaluation: 'hamming' (substitutions only) or 'edit' (substitutions + indels) note that debug prints tracebacks for all.
        gap_open_penalty: Gap open penalty for alignment
        gap_extend_penalty: Gap extend penalty for alignment
    Note:
        if the query region (spacer) is of a different length than the target region (area on the contig where it matched), the hamming distance is not meaningful.
    """
    logger.debug("[bold cyan]EXAMPLE ALIGNMENTS BY CLASSIFICATION[/bold cyan]")

    # Get examples for each classification type
    classifications_to_show = (
        alignment_classifications["classification"].unique().to_list()
    )

    all_run_tools = set(tools_results["tool"].unique().to_list()) if "tool" in tools_results.columns else set()

    # first display examples of each classification type
    for classification in classifications_to_show:
        examples = alignment_classifications.filter(
            pl.col("classification") == classification
        ).head(num_examples)

        total_count = alignment_classifications.filter(
            pl.col("classification") == classification
        ).height

        if examples.height == 0:
            classification_name = classification.upper().replace("_", " ")
            logger.debug(f"\n[dim]{classification_name} (0 found)[/dim]")
            continue

        logger.debug(
            f"\n[bold yellow]{classification.upper().replace('_', ' ')} ({total_count} total, showing {min(num_examples, examples.height)}):[/bold yellow]"
        )

        for idx, row in enumerate(examples.iter_rows(named=True), 1):
            spacer_id = row["spacer_id"]
            contig_id = row["contig_id"]
            start = row["start"]
            end = row["end"]
            strand = row["strand"]

            try:
                spacer_df = get_seq_from_fastx(
                    spacers_file,
                    [spacer_id],
                    return_df=True,
                    idcol="spacer_id",
                    seqcol="spacer_seq",
                )
                contig_df = get_seq_from_fastx(
                    contigs_file,
                    [contig_id],
                    return_df=True,
                    idcol="contig_id",
                    seqcol="contig_seq",
                )

                if spacer_df.height == 0 or contig_df.height == 0:
                    logger.warning(
                        f"Could not find sequences for {spacer_id} or {contig_id}"
                    )
                    continue

                spacer_seq = spacer_df["spacer_seq"][0]
                contig_seq = contig_df["contig_seq"][0]
            except Exception as e:
                logger.error(f"Error fetching sequences: {e}")
                continue

            # Find which tools reported this exact alignment (alignment_idx)
            alignment_idx = row.get("alignment_idx")
            region_idx = row.get("region_idx")
            
            tools_with_this_alignment = tools_results.filter(
                (pl.col("spacer_id") == spacer_id)
                & (pl.col("contig_id") == contig_id)
                & (pl.col("start") == start)
                & (pl.col("end") == end)
                & (pl.col("strand") == strand)
            )
            tool_names_this_alignment = (
                tools_with_this_alignment["tool"].unique().sort().to_list()
                if "tool" in tools_with_this_alignment.columns
                else []
            )
            
            # Find which tools reported this region (region_idx)
            tools_with_this_region = tools_results.filter(
                (pl.col("region_idx") == region_idx)
            ) if region_idx is not None and "region_idx" in tools_results.columns else pl.DataFrame()
            
            tool_names_this_region = (
                tools_with_this_region["tool"].unique().sort().to_list()
                if tools_with_this_region.height > 0 and "tool" in tools_with_this_region.columns
                else []
            )

            # When fetching from FASTA, contig_seq is full-length, so pass start/end
            edit_distance = test_alignment(
                spacer_seq,
                contig_seq,
                strand=strand,
                start=start,
                end=end,
                distance_metric="edit",
                gap_cost=gap_open_penalty,
                gap_extend=gap_extend_penalty,
            )
            # Hamming distance already includes length penalty in calculate_hamming_distance()
            hamming_distance = test_alignment(
                spacer_seq,
                contig_seq,
                strand=strand,
                start=start,
                end=end,
                distance_metric="hamming",
                gap_cost=gap_open_penalty,
                gap_extend=gap_extend_penalty,
            )

            gap_affine_distance: Optional[int] = None
            if distance_metric == "gap_affine":
                gap_affine_distance = test_alignment(
                    spacer_seq,
                    contig_seq,
                    strand=strand,
                    start=start,
                    end=end,
                    distance_metric="gap_affine",
                    gap_cost=gap_open_penalty,
                    gap_extend=gap_extend_penalty,
                )

            alignment_hamming = prettify_alignment_hamming(
                spacer_seq,
                contig_seq,
                strand=strand,
                start=start,
                end=end,
            )
            alignment_edit = prettify_alignment_edit(
                spacer_seq,
                contig_seq,
                strand=strand,
                start=start,
                end=end,
            )
            alignment_gap_affine = prettify_alignment_gap_affine(
                spacer_seq,
                contig_seq,
                strand=strand,
                start=start,
                end=end,
                gap_cost=gap_open_penalty,
                extend_cost=gap_extend_penalty,
            )

            strand_str = "(-)" if strand else "(+)"
            location_str = f"{contig_id}:{start}-{end} {strand_str}"
            
            logger.debug(f"[bold]Example {idx}:[/bold]")
            for metric_name, metric_value in [
                ("hamming", hamming_distance),
                ("edit", edit_distance),
                ("gap_affine", gap_affine_distance if gap_affine_distance is not None else edit_distance),
            ]:
                # Select alignment based on which metric we're displaying (not the global distance_metric)
                if metric_name == "gap_affine":
                    alignment_str = alignment_gap_affine
                elif metric_name == "edit":
                    alignment_str = alignment_edit
                elif metric_name == "hamming":
                    alignment_str = alignment_hamming
                lines = alignment_str.split("\n")
                # Add sequence IDs on the right side
                max_len = max(len(line) for line in lines)
                if distance_metric == metric_name:
                    logger.debug(f"[bold] USED {metric_name.capitalize()} :[/bold]")
                else:
                    logger.debug(f"[bold] (NOT USED) {metric_name.capitalize()}:[/bold]")

                logger.debug(f"  {lines[0]:<{max_len}}  [cyan]{spacer_id}[/cyan]")
                logger.debug(f"  {lines[1]:<{max_len}}  [dim](alignment)[/dim]")
                logger.debug(f"  {lines[2]:<{max_len}}  [cyan]{location_str}[/cyan]")

                # show distance
                logger.debug(
                    f"  {metric_name.capitalize()} distance: [{'bold green' if metric_value <= max_distance else 'bold red'}]{metric_value}[/{'bold green' if metric_value <= max_distance else 'bold red'}]"
                )

            logger.debug(f"  Max allowed mismatches: [bold]{max_distance}[/bold]")
            logger.debug(
                "  Metric used for validation: [bold]" +
                f"{distance_metric.capitalize()} distance:[/bold] " +
                f"{gap_affine_distance if distance_metric == 'gap_affine' else edit_distance if distance_metric == 'edit' else hamming_distance} "
            )

            # Display which tools found this alignment and region
            if tool_names_this_alignment:
                logger.debug(f"  Tools reporting this alignment_idx ({alignment_idx}): [magenta]{', '.join(tool_names_this_alignment)}[/magenta]")
            # if tool_names_this_region and tool_names_this_region != tool_names_this_alignment:
                logger.debug(f"  Tools reporting this region_idx ({region_idx}): [cyan]{', '.join(tool_names_this_region)}[/cyan]")
            
            # Highlight "planned" false negatives (in ground_truth but not in other tools)
            if "ground_truth" in tool_names_this_alignment:
                other_tools = [t for t in tool_names_this_alignment if t != "ground_truth"]
                if not other_tools:
                    logger.debug(f"  [bold red]FALSE NEGATIVE: Missed by all tools[/bold red]")
                else:
                    all_real_tools = set([t for t in tool_names_this_region if t != "ground_truth"])
                    missed_tools = all_real_tools - set(other_tools)
                    if missed_tools:
                        logger.debug(f"  [yellow]Missed by: {', '.join(sorted(missed_tools))}[/yellow]")


            # Highlight "non-planned" false negatives (not in ground_truth but  in other tools)
            if "ground_truth" not in tool_names_this_alignment:
                missed =  set(all_run_tools) - set(tool_names_this_alignment) 
                logger.debug(f"  [yellow]Missed by: {', '.join(sorted(missed))}[/yellow]")
    
    # Second, for each tool, display an example REGION that was completely missed (if any)
    logger.debug("\n[bold cyan]PER-TOOL COMPLETELY MISSED REGIONS:[/bold cyan]")
    
    for tool in sorted(all_run_tools):
        # Get all true positive regions (unique region_idx values)
        if "region_idx" not in alignment_classifications.columns:
            continue
            
        true_positive_regions = (
            alignment_classifications
            .filter(pl.col("classification").is_in(["positive_in_plan", "positive_not_in_plan"]))
            .select("region_idx")
            .unique()
        )
        
        if true_positive_regions.height == 0:
            continue
        
        # Get regions found by this specific tool (any alignment in the region counts)
        tool_regions = (
            tools_results
            .filter(pl.col("tool") == tool)
            .select("region_idx")
            .unique()
        )
        
        # Find regions completely missed by this tool (anti-join)
        missed_regions = true_positive_regions.join(
            tool_regions,
            on="region_idx",
            how="anti"
        )
        
        if missed_regions.height == 0:
            logger.debug(f"\n[bold green]{tool}: Found all {true_positive_regions.height} true positive regions[/bold green]")
            continue
        
        # Select one random missed region
        missed_region_idx = missed_regions.sample(1)["region_idx"][0]
        
        # Get one example alignment from this missed region
        region_alignments = alignment_classifications.filter(
            pl.col("region_idx") == missed_region_idx
        )
        example = region_alignments.sample(1).row(0, named=True)
        
        logger.debug(f"\n[bold red]{tool}: Completely missed {missed_regions.height}/{true_positive_regions.height} true positive regions[/bold red]")
        logger.debug(f"[bold red]Example region_idx {missed_region_idx} (alignment_idx {example.get('alignment_idx')}) missed by {tool}:[/bold red]")
        
        # Show which tools DID find this region
        tools_that_found_region = (
            tools_results
            .filter(pl.col("region_idx") == missed_region_idx)
            .select("tool")
            .unique()
            .sort("tool")
        )
        if tools_that_found_region.height > 0:
            found_by = tools_that_found_region["tool"].to_list()
            logger.debug(f"  [green]Region found by: {', '.join(found_by)}[/green]")
        
        # Display the example alignment from the missed region
        spacer_id = example["spacer_id"]
        contig_id = example["contig_id"]
        start = example["start"]
        end = example["end"]
        strand = example["strand"]
        classification = example["classification"]

        logger.debug(f"  Classification: [yellow]{classification}[/yellow]")
        logger.debug(f"  Spacer: [cyan]{spacer_id}[/cyan]")
        logger.debug(f"  Location: [cyan]{contig_id}:{start}-{end} {'(-)' if strand else '(+)'}[/cyan]")
        try:
            spacer_df = get_seq_from_fastx(
                spacers_file,
                [spacer_id],
                return_df=True,
                idcol="spacer_id",
                seqcol="spacer_seq",
            )
            contig_df = get_seq_from_fastx(
                contigs_file,
                [contig_id],
                return_df=True,
                idcol="contig_id",
                seqcol="contig_seq",
            )

            if spacer_df.height == 0 or contig_df.height == 0:
                logger.warning(
                    f"Could not find sequences for {spacer_id} or {contig_id}"
                )
                continue

            spacer_seq = spacer_df["spacer_seq"][0]
            contig_seq = contig_df["contig_seq"][0]
        except Exception as e:
            logger.error(f"Error fetching sequences: {e}")
            continue
        
        # Display alignments
        alignment_edit = prettify_alignment_edit(
            spacer_seq,
            contig_seq,
            strand=strand,
            start=start,
            end=end,
        )
        alignment_hamming = prettify_alignment_hamming(
            spacer_seq,
            contig_seq,
            strand=strand,
            start=start,
            end=end,
        )
        
        logger.debug(f"  Alignment (edit):\n{alignment_edit}")
        logger.debug(f"  Alignment (hamming):\n{alignment_hamming}")
        
        # Show which tools DID find this alignment
        tools_that_found = tools_results.filter(
            (pl.col("spacer_id") == spacer_id) &
            (pl.col("contig_id") == contig_id) &
            (pl.col("start") == start) &
            (pl.col("end") == end) &
            (pl.col("strand") == strand)
        )
        if tools_that_found.height > 0:
            found_by = tools_that_found["tool"].unique().sort().to_list()
            logger.debug(f"  [green]Found by: {', '.join(found_by)}[/green]")



def calculate_all_tool_performance(
    tools_results: pl.DataFrame,
    alignment_classifications: pl.DataFrame,
    ground_truth: pl.DataFrame,
    max_distance: int,
    distance_metric: str = "hamming",
) -> pl.DataFrame:
    """
    Calculate performance metrics for all tools at once using group_by aggregations.

    Always calculates both planned-only and augmented (including non-planned) metrics.

    Args:
        tools_results: DataFrame with all tool results (must have 'tool' column)
        alignment_classifications: DataFrame with classifications for unique alignments
                                  (columns: spacer_id, contig_id, start, end, strand, classification, recalculated_distance)
        ground_truth: Original ground truth DataFrame (needed to identify false negatives)
        max_distance: Maximum distance threshold
        distance_metric: Which distance metric to use for filtering ("hamming", "edit", "gap_affine")

    Returns:
        DataFrame with performance metrics for all tools (both planned-only and augmented)
    """
    # Select which distance column to use based on distance_metric
    distance_col = {
        "edit": "distance_edit",
        "gap_affine": "distance_gap_affine",
    }.get(distance_metric, "distance_hamming")
    
    # FIRST: Reclassify alignment_classifications based on the selected distance metric
    # This determines the global augmented GT for this metric
    reclassified_alignments = alignment_classifications.with_columns(
        pl.when(pl.col("classification") == "positive_in_plan")
        .then(
            pl.when(pl.col(distance_col) <= max_distance)
            .then(pl.lit("positive_in_plan"))
            .otherwise(pl.lit("invalid_alignment"))
        )
        .when(pl.col("classification") == "positive_not_in_plan")
        .then(
            pl.when(pl.col(distance_col) <= max_distance)
            .then(pl.lit("positive_not_in_plan"))
            .otherwise(pl.lit("invalid_alignment"))
        )
        .otherwise(pl.col("classification"))
        .alias(f"classification_{distance_metric}")
    )
    
    # Join tool results with reclassified alignments
    classified_results = tools_results.join(
        reclassified_alignments,
        on=["spacer_id", "contig_id", "start", "end", "strand"],
        how="left",
    ).with_columns(
        # Use the metric-specific classification
        pl.col(f"classification_{distance_metric}").alias("classification")
    )

    # Handle cases where some results don't have classifications (shouldn't happen, but defensive)
    # This can occur if alignment_classifications is incomplete
    if classified_results.filter(pl.col("classification").is_null()).height > 0:
        logger.warning(
            f"Found {classified_results.filter(pl.col('classification').is_null()).height} tool results without classifications - skipping them"
        )
        classified_results = classified_results.filter(
            pl.col("classification").is_not_null()
        )

    # If no classified results remain, return empty performance metrics
    if classified_results.height == 0:
        logger.warning(
            "No classified results found - returning empty performance metrics"
        )
        return pl.DataFrame(
            {
                "tool": [],
                "recall_planned": [],
                "recall_augmented": [],
                "all_true_positives": [],
                "planned_true_positives": [],
                "positives_not_in_plan": [],
                "invalid_alignments": [],
                "false_negatives_planned": [],
            }
        )

    # Ground truth counts (constant across all tools)
    planned_gt_count = ground_truth.height
    
    # For augmented GT, count verified non-planned alignments (overlap-merged) FOR THIS METRIC
    # Use the reclassified alignments, not the original
    not_in_plan_df = reclassified_alignments.filter(
        pl.col(f"classification_{distance_metric}") == "positive_not_in_plan"
    ).select(["spacer_id", "contig_id", "start", "end", "strand"])
    merged_not_in_plan = merge_overlapping_intervals(not_in_plan_df, tolerance=3)
    unique_positives_not_in_plan = merged_not_in_plan.height
    augmented_gt_count = planned_gt_count + unique_positives_not_in_plan

    # Aggregate by tool to count each classification type
    # Important: We need to count unique genomic locations, not duplicate tool reports
    # For positive_in_plan, use GROUND TRUTH coordinates (start_gt, end_gt) to avoid counting
    # the same GT location multiple times when a tool reports it with slightly different coordinates
    
    # For positives_not_in_plan, we need to merge overlapping intervals per tool
    # to avoid counting the same spurious match multiple times when reported with slight boundary differences
    tool_not_in_plan_counts = []
    for tool_name in classified_results["tool"].unique().to_list():
        tool_not_in_plan = (
            classified_results.filter(
                (pl.col("tool") == tool_name)
                & (pl.col("classification") == "positive_not_in_plan")
            )
            .select(["spacer_id", "contig_id", "start", "end", "strand"])
        )
        merged = merge_overlapping_intervals(tool_not_in_plan, tolerance=3)
        tool_not_in_plan_counts.append(
            {"tool": tool_name, "positives_not_in_plan": merged.height}
        )

    not_in_plan_df = pl.DataFrame(tool_not_in_plan_counts)

    performance = classified_results.group_by("tool").agg(
        [
            # Count unique GROUND TRUTH locations for positive_in_plan (use start_gt, end_gt)
            pl.struct(
                [
                    "spacer_id",
                    "contig_id",
                    "start_gt",
                    "end_gt",
                    "strand",
                    "classification",
                ]
            )
            .filter(pl.col("classification") == "positive_in_plan")
            .n_unique()
            .alias("planned_true_positives"),
            # For invalid alignments, use tool-reported coordinates (exact count is fine here)
            pl.struct(
                ["spacer_id", "contig_id", "start", "end", "strand", "classification"]
            )
            .filter(pl.col("classification") == "invalid_alignment")
            .n_unique()
            .alias("invalid_alignments"),
        ]
    )
    
    # Join the merged positives_not_in_plan counts
    performance = performance.join(not_in_plan_df, on="tool", how="left")

    # Add ground truth counts and calculate metrics
    # Note: augmented_gt is GLOBAL for this metric (not per-tool)
    performance = performance.with_columns(
        [
            pl.lit(planned_gt_count).alias("ground_truth_planned"),
            pl.lit(augmented_gt_count).alias("ground_truth_augmented"),
            (pl.col("planned_true_positives") + pl.col("positives_not_in_plan")).alias(
                "all_true_positives"
            ),
        ]
    )

    performance = performance.with_columns(
        [
            # False negatives: GT - TPs detected by this tool
            (pl.lit(planned_gt_count) - pl.col("planned_true_positives")).clip(lower_bound=0).alias(
                "false_negatives_planned"
            ),
            (pl.lit(augmented_gt_count) - pl.col("all_true_positives")).clip(lower_bound=0).alias(
                "false_negatives_augmented"
            ),
        ]
    )

    performance = performance.with_columns(
        [
            # Recall: TPs detected by tool / total GT
            (pl.col("planned_true_positives") / pl.lit(planned_gt_count))
            .fill_null(0.0)
            .alias("recall_planned"),
            (pl.col("all_true_positives") / pl.lit(augmented_gt_count))
            .fill_null(0.0)
            .alias("recall_augmented"),
        ]
    )
    return performance



def calculate_all_distances_for_alignment(
    spacer_seq: str,
    contig_seq: str,
    strand: bool,
    gap_open: int = 5,
    gap_extend: int = 5,
) -> dict:
    """
    Calculate all three distance metrics for a single alignment.
    
    Args:
        spacer_seq: Query sequence
        contig_seq: Target sequence (already trimmed to region)
        strand: Reverse complement flag
        gap_open: Gap opening penalty for gap-affine
        gap_extend: Gap extension penalty for gap-affine
        
    Returns:
        Dictionary with hamming, edit, and gap_affine distances
    """
    # Calculate hamming distance
    hamming_dist = calculate_hamming_distance(
        spacer_seq=spacer_seq,
        contig_seq=contig_seq,
        strand=strand,
        start=None,
        end=None,
        silent=True,
    )
    
    # Calculate edit distance (minimal)
    edit_dist = calculate_edit_distance(
        spacer_seq=spacer_seq,
        contig_seq=contig_seq,
        strand=strand,
        start=None,
        end=None,
    )
    
    # Calculate gap-affine distance
    gap_affine_dist = calculate_gap_affine_edit(
        spacer_seq=spacer_seq,
        contig_seq=contig_seq,
        strand=strand,
        start=None,
        end=None,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )
    
    return {
        'distance_hamming': hamming_dist,
        'distance_edit': edit_dist,
        'distance_gap_affine': gap_affine_dist,
    }

