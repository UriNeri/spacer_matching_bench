#!/usr/bin/env python3
"""
Convert Google Docs exported markdown (test.md) back to Quarto-compatible .qmd format.
This script:
1. Adds the YAML front matter from the original index.qmd
2. Converts superscript citation numbers (^1,2^) to Quarto [@key1; @key2] format
3. Converts heading syntax from "- - " prefix to proper ## / ### / #### headings
4. Converts Google Docs math ($$O\left(...\right)$$) to Quarto inline/block math
5. Removes Google Docs artifacts ([]{#anchor-N}, ::: {} blocks, reference list at end)
6. Converts Google Docs figure/table blocks to Quarto format
"""

import re
import sys

# --- Citation number to bib key mapping ---
# Built from the reference list at the bottom of test.md matched to references.bib keys
CITATION_MAP = {
    1: "Mojica_2005",
    2: "CRISPR_review",
    3: "CRISPR_classification",
    4: "Edwards2015_phage_host",
    5: "CRISPR_gene_editing_review",
    6: "soto_perez_crispr_2019",
    7: "Frith_2010",
    8: "Morgulis_2006",
    9: "Stern_2010",
    10: "Biswas2013",
    11: "Altschul1990_blast",
    12: "Shah2018",
    13: "Madden2018",
    14: "Needleman1970_global_alignment",
    15: "Smith1981_local_alignment",
    16: "Myers1999_bitparallel",
    17: "Langmead2009_bowtie",
    18: "Langmead2012_bowtie2",
    19: "Li2018_minimap2",
    20: "Sahlin2022_strobealign",
    21: "Steinegger2017_mmseqs2",
    22: "Jain_2016",
    23: "Deveau2008",
    24: "Semenova2011",
    25: "Fineran2014",
    26: "Schelling2023",
    27: "Lee_2012_mutation_rate_Ecoli",
    28: "Kucukyildirim_2021_high_indel_rate",
    29: "Hatfull_2011_bacteriophages_genomes",
    30: "Ha_2018",
    31: "Paez_Espino_2015",
    32: "Kupczok_2018_phage_genome_evolution",
    33: "Shmakov_2017",
    34: "camargo_img_vr4_2023",
    35: "Dion_2021",
    36: "Zhang_2021",
    37: "Roux2023_iphop",
    38: "camargo_genomad_2024",
    39: "zenodo_doi",
    40: "edgar_piler_cr_2007",
    41: "bland_crt_2007",
    42: "Daily2016_parasail",
    43: "Šošić_Šikić_2017_edlib",
    44: "Peter_hyperfine_2023",
    45: "SLURM_2002",
    46: "conda",
    47: "mamba",
    48: "Kosmopoulos_2023",
    49: "Maier_2018",
    50: "Zhang_2025",
    51: "Turgeman_Grott_2018",
    52: "Shmakov_2023",
    53: "Mitrofanov2025",
    54: "Vink2021",
}

YAML_HEADER = """---
title: "Computational Tool Choice Impacts CRISPR Spacer-Protospacer Detection"
author:
  - name: Uri Neri*^1^
    corresponding: true
  - name: Antonio Pedro Camargo^1^
  - name: Brian Bushnell^1^
  - name: Rick Beeloo^2^
  - name: Simon Roux^1^
format:
  pdf:
    documentclass: article
    geometry:
      - margin=1in
    fig-format: pdf
    embed-resources: true
    keep-tex: true
    fig-pos: 'H'
    number-sections: true
    cite-method: biblatex
    bibliography: references.bib
    pdf-engine: xelatex
  html:
    toc: true
    toc-depth: 3
    number-sections: true
  docx: default
execute:
  echo: false
  warning: false
jupyter: python3
---
"""


def convert_citations(text):
    """Convert superscript citations like ^1,2^ or ^23-26^ to [@key1; @key2] format."""

    def replace_citation(match):
        raw = match.group(1)
        nums = []
        # Handle ranges and comma-separated values: e.g. "23-26" or "1,2" or "29,32"
        parts = raw.split(",")
        for part in parts:
            part = part.strip()
            if "-" in part:
                start, end = part.split("-", 1)
                start, end = int(start.strip()), int(end.strip())
                nums.extend(range(start, end + 1))
            else:
                nums.append(int(part))

        keys = []
        for n in nums:
            if n in CITATION_MAP:
                keys.append(f"@{CITATION_MAP[n]}")
            else:
                keys.append(f"@UNKNOWN_REF_{n}")
        return "[" + "; ".join(keys) + "]"

    # Match ^N^ or ^N,M^ or ^N-M^ patterns (citation superscripts)
    # Be careful not to match author superscript affiliations like ^1^ in author lines
    # Only match when preceded by a word char, punctuation, or space (not in title/author context)
    text = re.sub(r"\^(\d+(?:[,-]\d+)*)\^", replace_citation, text)
    return text


def convert_headings(text):
    """Convert Google Docs heading syntax to Quarto markdown headings."""
    lines = text.split("\n")
    result = []
    i = 0
    while i < len(lines):
        line = lines[i]

        # "- - Results" or "- - Discussion" -> ## heading
        # "- - Abstract" -> ## Abstract {#sec-abstract}
        # "    - Tool Selection" -> ### heading (indented sub of ##)
        # "- - - heading" -> ### heading
        # "- - - - heading" -> #### heading

        # Pattern: "- - - - Real datasets" -> ####
        if re.match(r"^- - - - ", line):
            heading_text = line.replace("- - - - ", "", 1).strip()
            heading_id = make_section_id(heading_text)
            result.append(f"#### {heading_text} {{#{heading_id}}}")
            i += 1
            continue

        # Pattern: "- - - Heading" -> ###
        if re.match(r"^- - - ", line):
            heading_text = line.replace("- - - ", "", 1).strip()
            heading_id = make_section_id(heading_text)
            result.append(f"### {heading_text} {{#{heading_id}}}")
            i += 1
            continue

        # Pattern: "- - Heading" -> ##
        if re.match(r"^- - ", line):
            heading_text = line.replace("- - ", "", 1).strip()
            heading_id = make_section_id(heading_text)
            result.append(f"## {heading_text} {{#{heading_id}}}")
            i += 1
            continue

        # Pattern: "    - Tool Selection" (preceded by "- - Methods" style heading)
        if re.match(r"^    - ", line):
            heading_text = line.strip().lstrip("- ").strip()
            heading_id = make_section_id(heading_text)
            result.append(f"### {heading_text} {{#{heading_id}}}")
            i += 1
            continue

        # "      database" continuation of "- - - - " heading with wrapping
        # Not needed, multi-line headings are unlikely

        result.append(line)
        i += 1

    return "\n".join(result)


# Section ID mapping - match the original index.qmd IDs where possible
SECTION_IDS = {
    "abstract": "sec-abstract",
    "introduction": "sec-introduction",
    "methods": "sec-methods",
    "tool selection": "sec-tool-selection",
    "data generation and acquisition": "sec-data",
    "real datasets": "sec-real-data",
    "synthetic dataset generation": "sec-synthetic-data",
    "semi-synthetic dataset": "sec-semi-synthetic",
    "coordinate tolerance and unique region counting": "sec-coordinate-tolerance",
    "alignment verification and distance metric calculation": "sec-alignment-recalc",
    "performance definitions and calculation": "sec-performance-calc",
    "non-planned match rate estimation": "sec-non-planned-rate",
    "computational resource and runtime tracking": "sec-resource-tracking",
    "versioning and reproducibility": "sec-reproducibility",
    "extensibility": "sec-extensibility",
    "results": "sec-results",
    "selection of distance metric and threshold values": "sec-distance-metric-selection",
    "tool performance across distance thresholds": "sec-distance-performance",
    "computational resource requirements and scalability": "sec-resource-usage",
    "performance as a function of query (spacer) abundance in reference database": "sec-abundance-performance",
    "discussion": "sec-discussion",
    "distance metric choice and practical thresholds": "sec-distance-thresholds",
    "tool performance and algorithmic considerations": "sec-tool-performance-discussion",
    "biological interpretation of spacer-protospacer matches": "sec-biological-interpretation",
    "study limitations": "sec-limitations",
    "conclusion": "sec-conclusion",
    "code and data availability": "sec-code-data",
    "acknowledgements": "sec-acknowledgements",
}


def make_section_id(heading_text):
    key = heading_text.lower().strip()
    # Remove any trailing anchor refs
    key = re.sub(r"\[\]\{#[^}]+\}", "", key).strip()
    if key in SECTION_IDS:
        return SECTION_IDS[key]
    # Auto-generate from heading text
    slug = re.sub(r"[^a-z0-9]+", "-", key).strip("-")
    return f"sec-{slug}"


def remove_anchor_spans(text):
    """Remove Google Docs anchor spans like []{#anchor} and []{#anchor-N}."""
    return re.sub(r"\[\]\{#anchor(?:-\d+)?\}", "", text)


def remove_reference_list(text):
    """Remove the rendered reference list at the end (Google Docs export artifact)."""
    # The reference list starts with a series of ::: {} blocks containing numbered refs
    # Find the first ::: {} block that contains a numbered reference
    lines = text.split("\n")
    # Find the start of the reference list: look for first "::: {}" followed by "N\."
    ref_start = None
    for i, line in enumerate(lines):
        if line.strip() == "::: {}":
            # Check if next non-empty line starts with a number and period
            for j in range(i + 1, min(i + 3, len(lines))):
                if re.match(r"^\d+\\?\.", lines[j].strip()):
                    ref_start = i
                    break
            if ref_start is not None:
                break

    if ref_start is not None:
        # Trim everything from ref_start onward
        text = "\n".join(lines[:ref_start]).rstrip() + "\n"
    return text


def convert_google_doc_math(text):
    """Convert Google Docs math artifacts to clean Quarto math."""
    # $$O\left( \mathit{mn} \right)$$ -> $O(mn)$
    # These are inline but wrapped in $$...$$, convert to $...$
    # Also clean up \mathit{}, \left( \right), etc.

    def clean_math(match):
        content = match.group(1)
        # Remove \mathit{} wrapper
        content = re.sub(r"\\mathit\{([^}]+)\}", r"\1", content)
        # Remove \left and \right
        content = content.replace(r"\left(", "(").replace(r"\right)", ")")
        content = content.replace(r"\left", "").replace(r"\right", "")
        # Clean up excessive spaces
        content = re.sub(r"\s+", " ", content).strip()
        return f"${content}$"

    # Convert inline $$ ... $$ that appear mid-paragraph (not on their own line)
    # We look for $$ ... $$ patterns that don't span multiple lines
    text = re.sub(r"\$\$([^$\n]+?)\$\$", clean_math, text)
    return text


def convert_google_figure_blocks(text):
    """Convert Google Docs figure blocks to Quarto figure syntax."""
    # The Google Docs export wraps figures in table-like blocks:
    # +----------+
    # | ![](path){...} caption |
    # +----------+

    # Figure 1 - background normalization
    text = re.sub(
        r"\+[-]+\+\n\| !\[\]\(Pictures/[^)]+\)\{[^}]*\}\[\]\{#anchor-18\}\*Figure 1:([^*]+)\*\s*\|\n\+[-]+\+",
        r'![Non-planned alignment ("background") frequency: mean ± SD across 4 fully synthetic simulations (50,000--100,000 spacers × 5,000--20,000 contigs). (A) Absolute count of non-planned matches at each distance threshold (0--5). (B) Non-planned matches normalized per spacer per Gbp of target sequence. (C) Non-planned matches normalized per Gbp$^2$ of total search space (spacer-bp × contig-bp). Blue bars: hamming distance (substitutions only). Orange bars: edit distance (substitutions + indels). Error bars: standard deviation across simulations. At hamming ≤3, the per Gbp$^2$ rate is approximately 40,000, while at edit ≤3 the rate is approximately 270,000--290,000 (a 6--8 fold increase). Full per-simulation breakdowns are provided in Supplementary Note 3.](figures/search_space_background_normalization.svg){#fig-background-normalization}',
        text,
        flags=re.DOTALL,
    )

    # Figure 2 - tool performance heatmap
    text = re.sub(
        r"\+[-]+\+\n\| !\[\]\(Pictures/[^)]+\)\{[^}]*\}\[\]\{#anchor-21\}\*Figure 2:([^*]+)\*\s*\|\n\+[-]+\+",
        r'![Tool recall (detection fraction of unique valid aligned regions) across distance thresholds. Upper row shows IMG/VR4 (fractions) results and lower row shows synthetic datasets; columns compare hamming and edit distance analyses. An asterix "*" is noted after tools that did failed or timed out, and the total number of sucessful runs (subsampled fraction for the IMG/VR sets, and independent measurements for the simulated set) are noted in brackets after each tool\'name, and the value in each cell is the mean of those run specific recall values. Note - the values are at exact distance (unlike at a "max" distance), i.e. regions that aligned with distance n are not considered for distance n+1. ](figures/heatmap_real_and_sim_hamming_edit_grid.svg){#fig-tool-performance}',
        text,
        flags=re.DOTALL,
    )

    # Figure 3 - resource usage
    text = re.sub(
        r"\+[-]+\+\n\| !\[\]\(Pictures/[^)]+\)\{[^}]*\}\[\]\{#anchor-23\}\*Figure 3:([^*]+)\*\s*\|\n\+[-]+\+",
        r"![Total CPU time scaling with dataset size (log-log scale). (A) Real IMG/VR4 subsampled datasets (fractions 0.001--1.0, approximately 10 Mbp to 18.9 Gbp target sequence). (B) Simulated datasets with varying numbers of spacers and contigs (see Supplementary Table S4 for dataset details). Marker shape and color encode tool identity. BLASTn runtime for fraction_1 was extrapolated from partial completion (1.69M of 3.83M spacers processed within the 72h wall time limit; see Supplementary Table S8 note).](./figures/resource_usage/resource_usage_cpu_2panel.svg){#fig-resource-usage}",
        text,
        flags=re.DOTALL,
    )

    return text


def convert_google_table_blocks(text):
    """Convert Google Doc Table 1 block and Table 2 surrounding markup."""
    # Table 1 is in a boxed format in the Google Doc export - replace with
    # a reference to the original table from index.qmd
    # The box starts with +------+ and contains the table caption
    # We'll replace it with a placeholder note
    text = re.sub(
        r"\+[-]+\+\n\| \*Table 1: Evaluated Tools.*?\+[-]+\+",
        "<!-- Table 1: See index.qmd for the full pipe-table version of Table 1 (Evaluated Tools) -->",
        text,
        flags=re.DOTALL,
    )

    # Table 2 caption box
    text = re.sub(
        r"\+[-]+\+\n\| \*Table 2: Comparison of non-planned.*?\+[-]+\+",
        "<!-- Table 2 caption: See index.qmd for formatted table -->",
        text,
        flags=re.DOTALL,
    )

    return text


def convert_internal_refs(text):
    """Convert Google Docs internal references to Quarto cross-references."""
    # [*Figure 1*](#1hae98wie9xx) -> @fig-background-normalization
    # [*figure 2*](#xu0ta4x6wnu2) -> @fig-tool-performance
    # [Figure 3](#cm30dc28p9tj) -> @fig-resource-usage
    # [*figure 3*](#asqvfn36xcs3) -> @fig-resource-usage
    # [*Table 2*](#ec5v49j2nfmg) -> @tbl-semisynthetic-hamming3
    # [Table 2](#l1nml75n6754) -> @tbl-semisynthetic-hamming3
    # [Section 4.4](#1027qr61jhwv) -> @sec-abundance-performance

    # Generic figure references
    text = re.sub(
        r"\[\*?[Ff]igure 1\*?\]\(#[^)]+\)", "@fig-background-normalization", text
    )
    text = re.sub(r"\[\*?[Ff]igure 2\*?\]\(#[^)]+\)", "@fig-tool-performance", text)
    text = re.sub(r"\[\*?[Ff]igure 3\*?\]\(#[^)]+\)", "@fig-resource-usage", text)
    text = re.sub(r"\[\*?[Tt]able 2\*?\]\(#[^)]+\)", "@tbl-semisynthetic-hamming3", text)
    text = re.sub(r"\[Section 4\.4\]\(#[^)]+\)", "@sec-abundance-performance", text)

    # Fix ?@ref patterns that appear in the Google Doc (broken quarto refs)
    # These should become @ref
    text = re.sub(r"\*\*\?@([\w-]+)\*\*", r"@\1", text)

    return text


def remove_title_and_authors(text):
    """Remove the title line and author lines that will be in YAML header."""
    lines = text.split("\n")
    # Skip first few lines: title, blank, authors, affiliations
    # Title is line 0: "Computational Tool Choice..."
    # Authors are next few lines up to and including the affiliation/email block
    # Find where the actual content starts (after affiliations)
    start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("\\* Uri Neri"):
            start = i + 1
            break

    # Skip any blank lines after authors
    while start < len(lines) and lines[start].strip() == "":
        start += 1

    return "\n".join(lines[start:])


def clean_remaining_artifacts(text):
    """Clean up remaining Google Docs artifacts."""
    # Remove empty ::: {} blocks
    text = re.sub(r":::\s*\{\}\s*\n:::\s*\n?", "", text)
    text = re.sub(r":::\s*\{\}\s*$", "", text, flags=re.MULTILINE)

    # Clean up Table 2 simple table format (dataset characteristics)
    # The simple table with --- separators is fine for Quarto, keep it

    return text


def main():
    input_file = "/clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_bench/draft/main/test.md"
    output_file = "/clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_bench/draft/main/test_converted.qmd"

    with open(input_file, "r") as f:
        text = f.read()

    # Step 1: Remove the title and author block (will be in YAML)
    text = remove_title_and_authors(text)

    # Step 2: Remove the reference list at the end
    text = remove_reference_list(text)

    # Step 3: Convert Google Docs figure blocks to Quarto format
    text = convert_google_figure_blocks(text)

    # Step 4: Convert table blocks
    text = convert_google_table_blocks(text)

    # Step 5: Convert internal cross-references
    text = convert_internal_refs(text)

    # Step 6: Convert citation superscripts to Quarto format
    text = convert_citations(text)

    # Step 7: Convert heading syntax
    text = convert_headings(text)

    # Step 8: Remove anchor spans
    text = remove_anchor_spans(text)

    # Step 9: Convert math notation
    text = convert_google_doc_math(text)

    # Step 10: Clean remaining artifacts
    text = clean_remaining_artifacts(text)

    # Add YAML header
    output = YAML_HEADER.lstrip() + "\n" + text

    with open(output_file, "w") as f:
        f.write(output)

    print(f"Converted {input_file} -> {output_file}")
    print(f"Output size: {len(output)} bytes, {output.count(chr(10))} lines")


if __name__ == "__main__":
    main()
