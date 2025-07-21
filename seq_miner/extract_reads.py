import pysam
from Bio import SeqIO
from statistics import mean
import os

def read_id_list(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}

def extract_from_bam(input_bam, output_bam, read_ids=None, min_qscore=0, min_length=0):
    read_ids = set(read_ids) if read_ids else None
    low_q = []
    short = []
    passed = []

    with pysam.AlignmentFile(input_bam, "rb") as infile, \
         pysam.AlignmentFile(output_bam, "wb", template=infile) as out_passed:

        for read in infile:
            if read_ids and read.query_name not in read_ids:
                continue

            seq = read.query_sequence
            qual = read.query_qualities

            if qual:
                mean_q = mean(qual)
            else:
                mean_q = 0

            if len(seq) < min_length:
                short.append(read)
            elif mean_q < min_qscore:
                low_q.append(read)
            else:
                passed.append(read)
                out_passed.write(read)

    return passed, low_q, short

def extract_from_fastq(input_fastq, output_fastq, read_ids=None, min_qscore=0, min_length=0):
    read_ids = set(read_ids) if read_ids else None
    passed = []
    low_q = []
    short = []

    with open(output_fastq, "w") as out_passed:
        for record in SeqIO.parse(input_fastq, "fastq"):
            if read_ids and record.id not in read_ids:
                continue

            avg_q = mean(record.letter_annotations["phred_quality"])
            seq_len = len(record.seq)

            if seq_len < min_length:
                short.append(record)
            elif avg_q < min_qscore:
                low_q.append(record)
            else:
                passed.append(record)
                SeqIO.write(record, out_passed, "fastq")

    return passed, low_q, short

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract and filter reads from BAM or FASTQ.")
    parser.add_argument("--input", "-i", required=True, help="Input BAM or FASTQ file")
    parser.add_argument("--output", "-o", required=True, help="Output file for passed reads")
    parser.add_argument("--format", "-f", required=True, choices=["bam", "fastq"], help="Input file format")
    parser.add_argument("--read-ids", "-r", help="File with list of read IDs to extract (optional)")
    parser.add_argument("--min-qscore", type=float, default=0, help="Minimum mean Q score")
    parser.add_argument("--min-length", type=int, default=0, help="Minimum read length")

    args = parser.parse_args()

    read_ids = read_id_list(args.read_ids) if args.read_ids else None

    if args.format == "bam":
        passed, low_q, short = extract_from_bam(args.input, args.output, read_ids, args.min_qscore, args.min_length)
    else:
        passed, low_q, short = extract_from_fastq(args.input, args.output, read_ids, args.min_qscore, args.min_length)

    print(f"Passed reads: {len(passed)}")
    print(f"Low-quality reads: {len(low_q)}")
    print(f"Short reads: {len(short)}")
