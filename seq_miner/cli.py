from .extract_reads import extract_from_bam, extract_from_fastq, read_id_list
import argparse
import textwrap
import sys
import json
import csv
from datetime import datetime

VERSION = "1.0.0"  # update with your version

def write_summary(file_path, summary_data):
    if file_path.endswith(".json"):
        with open(file_path, "w") as f:
            json.dump(summary_data, f, indent=2)
    elif file_path.endswith(".csv"):
        with open(file_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summary_data.keys())
            writer.writeheader()
            writer.writerow(summary_data)
    else:
        print(f"[WARNING] Unsupported summary format: {file_path}")

def log_output(logfile, message):
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    with open(logfile, "a") as f:
        f.write(f"{timestamp} {message}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Extract and filter reads from BAM or FASTQ files based on read ID, quality score, and length.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
            Examples:
              # Extract reads by ID from BAM
              seq-miner -i input.bam -o passed.bam -f bam -r read_ids.txt --min-qscore 10 --min-length 200

              # Filter FASTQ reads (parallel)
              seq-miner -i input.fastq.gz -o filtered.fastq -f fastq --min-qscore 15 --min-length 1000 --threads 4 --summary result.json
        """)
    )

    parser.add_argument("-i", "--input", required=True, help="Input BAM or FASTQ file")
    parser.add_argument("-o", "--output", required=True, help="Output file for passed reads")
    parser.add_argument("-f", "--format", required=True, choices=["bam", "fastq"], help="Input file format")
    parser.add_argument("-r", "--read-ids", help="File with list of read IDs (one per line, optional)")
    parser.add_argument("--min-qscore", type=float, default=0.0, help="Minimum mean quality score per read")
    parser.add_argument("--min-length", type=int, default=0, help="Minimum read length")
    parser.add_argument("--threads", type=int, help="Number of CPU threads to use (only for FASTQ)")
    parser.add_argument("--summary", help="Output summary to CSV or JSON file")
    parser.add_argument("--log", help="Log file path")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    args = parser.parse_args()
    read_ids = read_id_list(args.read_ids) if args.read_ids else None

    def verbose(msg):
        if args.verbose:
            print(f"[INFO] {msg}")
        if args.log:
            log_output(args.log, msg)

    verbose(f"Format        : {args.format}")
    verbose(f"Input file    : {args.input}")
    verbose(f"Output file   : {args.output}")
    verbose(f"Read ID file  : {args.read_ids or 'None'}")
    verbose(f"Min Q-score   : {args.min_qscore}")
    verbose(f"Min length    : {args.min_length}")
    if args.threads:
        verbose(f"Threads       : {args.threads}")
    if args.log:
        verbose(f"Log enabled   : {args.log}")

    if args.format == "bam":
        passed, low_q, short = extract_from_bam(
            args.input, args.output, read_ids,
            min_qscore=args.min_qscore,
            min_length=args.min_length
        )
    else:
        passed, low_q, short = extract_from_fastq(
            args.input, args.output, read_ids,
            min_qscore=args.min_qscore,
            min_length=args.min_length,
            threads=args.threads
        )

    summary = {
        "input": args.input,
        "output": args.output,
        "format": args.format,
        "read_id_file": args.read_ids or "None",
        "min_qscore": args.min_qscore,
        "min_length": args.min_length,
        "threads": args.threads if args.format == "fastq" else "N/A",
        "passed_reads": len(passed),
        "low_quality_reads": len(low_q),
        "short_reads": len(short),
    }

    print("\nSummary:")
    print(f"Passed reads     : {summary['passed_reads']}")
    print(f"Low-quality reads : {summary['low_quality_reads']}")
    print(f"Short reads      : {summary['short_reads']}")

    if args.summary:
        write_summary(args.summary, summary)
        print(f"\n Summary written to: {args.summary}")

    if args.log:
        log_output(args.log, f"Job complete. {summary['passed_reads']} passed, {summary['low_quality_reads']} low Q, {summary['short_reads']} short.")
