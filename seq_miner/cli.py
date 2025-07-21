from .extract_reads import extract_from_bam, extract_from_fastq, read_id_list
import argparse
import textwrap

def main():
    parser = argparse.ArgumentParser(
        description="Extract and filter reads from BAM or FASTQ files based on read ID, quality score, and length.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
            Examples:
              # Extract reads by ID from BAM
              seq-miner -i input.bam -o passed.bam -f bam -r read_ids.txt --min-qscore 10 --min-length 200

              # Filter FASTQ reads (parallel)
              seq-miner -i input.fastq.gz -o filtered.fastq -f fastq --min-qscore 15 --min-length 1000 --threads 4
        """)
    )

    parser.add_argument("-i", "--input", required=True, help="Input BAM or FASTQ file")
    parser.add_argument("-o", "--output", required=True, help="Output file for passed reads")
    parser.add_argument("-f", "--format", required=True, choices=["bam", "fastq"], help="Input file format")
    parser.add_argument("-r", "--read-ids", help="File with list of read IDs (one per line, optional)")
    parser.add_argument("--min-qscore", type=float, default=0.0, help="Minimum mean quality score per read")
    parser.add_argument("--min-length", type=int, default=0, help="Minimum read length")
    parser.add_argument("--threads", type=int, help="Number of CPU threads to use (only for FASTQ)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()
    read_ids = read_id_list(args.read_ids) if args.read_ids else None

    if args.verbose:
        print(f"[INFO] Format        : {args.format}")
        print(f"[INFO] Input file    : {args.input}")
        print(f"[INFO] Output file   : {args.output}")
        print(f"[INFO] Read ID file  : {args.read_ids or 'None'}")
        print(f"[INFO] Min Q-score   : {args.min_qscore}")
        print(f"[INFO] Min length    : {args.min_length}")
        if args.threads:
            print(f"[INFO] Threads       : {args.threads}")

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

    print(f"\nSummary:")
    print(f"Passed reads     : {len(passed)}")
    print(f"Low-quality reads : {len(low_q)}")
    print(f"Short reads      : {len(short)}")
