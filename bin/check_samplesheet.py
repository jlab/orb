#!/usr/bin/env python
"""Provide a command line tool to validate and transform tabular samplesheets."""

import argparse
import csv
import logging
import sys
from pathlib import Path

logger = logging.getLogger()

def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (file object): A handle to a file object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.
    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


class RowChecker:
    def __init__(self, sample_col="sample", assemblers_col="assemblers", **kwargs):
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._assemblers_col = assemblers_col
        self.modified = []

    def validate_and_transform(self, row):
        self._validate_sample(row)
        self._validate_assemblers(row)
        self.modified.append(row)

    def _validate_sample(self, row):
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_assemblers(self, row):
        assemblers = row[self._assemblers_col].split(';')
        if not assemblers:
            raise AssertionError("At least one assembler is required.")
        row[self._assemblers_col] = assemblers

def check_samplesheet(file_in, file_out):
    required_columns = {"sample", "assemblers"}
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        print(reader.fieldnames)
        print(reader.line_num)
        print(sniff_format(in_handle))
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        checker = RowChecker()
        print("test")
        for row in reader:
            print(row)
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
    header = list(reader.fieldnames)
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            row['assemblers'] = ';'.join(row['assemblers'])
            writer.writerow(row)

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)

# Rest of the code remains the same.

if __name__ == "__main__":
    sys.exit(main())
