# seq_miner/cli.py

import argparse
import sys
from .extract_reads import extract_from_fastq, extract_from_bam, read_id_list
from .__version__ import __version__

def main():
    parser = argparse.ArgumentParser(
        prog='seq-miner',
        description='Extract sequence reads from FASTQ or BAM files by ID, Q-score, length, and barcode classification.',
        epilog='Example: seq-miner -i input.fastq -o output_dir -f fastq --min-qscore 12 --min-length 800 --verbose'
    )

    parser.add_argument('-i', '--input', required=True, help='Input FASTQ or BAM file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-f', '--format', required=True, choices=['fastq', 'bam'], help='Input format: fastq or bam')
    parser.add_argument('--id-file', help='Optional file containing read IDs to extract (one per line)')
    parser.add_argument('--min-qscore', type=float, default=0, help='Minimum mean quality score filter')
    parser.add_argument('--min-length', type=int, default=0, help='Minimum read length filter')
    parser.add_argument('--log', help='Optional log file path')
    parser.add_argument('--threads', type=int, default=None, help='Number of parallel processes to use')
    parser.add_argument('--verbose', action='store_true', help='Show detailed progress during processing')
    parser.add_argument('--version', action='version', version=f'seq-miner {__version__}')

    args = parser.parse_args()

    read_ids = read_id_list(args.id_file) if args.id_file else None

    print(f"[INFO] Format        : {args.format}")
    print(f"[INFO] Input file    : {args.input}")
    print(f"[INFO] Output file   : {args.output}")
    print(f"[INFO] Read ID file  : {args.id_file}")
    print(f"[INFO] Min Q-score   : {args.min_qscore}")
    print(f"[INFO] Min length    : {args.min_length}")
    if args.log:
        print(f"[INFO] Log enabled   : {args.log}")

    if args.format == "fastq":
        extract_from_fastq(
            input_fastq=args.input,
            output_dir=args.output,
            read_ids=read_ids,
            min_qscore=args.min_qscore,
            min_length=args.min_length,
            threads=args.threads,
            log_file=args.log,
            verbose=args.verbose,
        )
    elif args.format == "bam":
        extract_from_bam(
            input_bam=args.input,
            output_dir=args.output,
            read_ids=read_ids,
            min_qscore=args.min_qscore,
            min_length=args.min_length,
            log_file=args.log,
            verbose=args.verbose,
        )
    else:
        print(f"[ERROR] Unsupported format: {args.format}")
        sys.exit(1)
