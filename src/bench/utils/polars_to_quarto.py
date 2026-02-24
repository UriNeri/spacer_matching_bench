#!/usr/bin/env python3
"""Convert Polars text tables (box-drawing format) to Quarto markdown pipe tables.

Usage:
    # from stdin (paste or pipe)
    cat table.txt | python scripts/polars_to_quarto.py
    python scripts/polars_to_quarto.py < table.txt

    # from file
    python scripts/polars_to_quarto.py table.txt

    # with options
    python scripts/polars_to_quarto.py table.txt --align right --float-fmt "{:.2f}" --caption "My table"

    # copy result to clipboard (if xclip available)
    python scripts/polars_to_quarto.py table.txt | xclip -selection clipboard
"""
import argparse
import sys
from pathlib import Path

# allow running from project root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from bench.utils.functions import polars_text_to_quarto_md


def main():
    parser = argparse.ArgumentParser(
        description="Convert Polars text table to Quarto markdown pipe table."
    )
    parser.add_argument(
        "infile",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Input file with Polars text table (default: stdin)",
    )
    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file (default: stdout)",
    )
    parser.add_argument(
        "--align",
        default="left",
        help="Column alignment: left, right, center, or per-column string of l/r/c (default: left)",
    )
    parser.add_argument(
        "--float-fmt",
        default="{:.4f}",
        help="Format string for float values (default: '{:.4f}')",
    )
    parser.add_argument(
        "--caption",
        default=None,
        help="Optional table caption",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Optional Quarto cross-reference label (e.g. tbl-my-table)",
    )
    args = parser.parse_args()

    text = args.infile.read()
    if not text.strip():
        print("Error: empty input", file=sys.stderr)
        sys.exit(1)

    result = polars_text_to_quarto_md(
        text,
        align=args.align,
        float_fmt=args.float_fmt,
        caption=args.caption,
        label=args.label,
    )
    args.outfile.write(result)


if __name__ == "__main__":
    main()
