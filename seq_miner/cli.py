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
        
              # Filter FASTQ reads
              seq-miner -i input.fastq -o passed.fastq -f fastq --min-qscore 15 --min-length 1000
        """)
    )

    parser.add_argument("-i", "--input", required=True, help="Input BAM or FASTQ file")
    parser.add_argument("-o", "--output", required=True, help="Output file for passed reads")
    parser.add_argument("-f", "--format", required=True, choices=["bam", "fastq"], help="Input file format")
    parser.add_argument("-r", "--read-ids", help="Path to file with read IDs (one per line, optional)")
    parser.add_argument("--min-qscore", type=float, default=0.0, help="Minimum mean quality score per read")
    parser.add_argument("--min-length", type=int, default=0, help="Minimum read length")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()
    read_ids = read_id_list(args.read_ids) if args.read_ids else None

    if args.verbose:
        print(f"[INFO] Format: {args.format}")
        print(f"[INFO] Input file: {args.input}")
        print(f"[INFO] Output file: {args.output}")
        print(f"[INFO] Read ID filtering: {'enabled' if read_ids else 'disabled'}")
        print(f"[INFO] Min Q-score: {args.min_qscore}")
        print(f"[INFO] Min read length: {args.min_length}")

    if args.format == "bam":
        passed, low_q, short = extract_from_bam(
            args.input, args.output, read_ids, args.min_qscore, args.min_length
        )
    else:
        passed, low_q, short = extract_from_fastq(
            args.input, args.output, read_ids, args.min_qscore, args.min_length
        )

    print(f"\nSummary:")
    print(f"Passed reads: {len(passed)}")
    print(f"Low-quality reads: {len(low_q)}")
    print(f"Short reads: {len(short)}")
