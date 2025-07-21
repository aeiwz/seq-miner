from Bio import SeqIO
from statistics import mean
from multiprocessing import Pool, cpu_count
from itertools import islice
import gzip
import os


def read_id_list(path):
    """Read a list of read IDs from a file, one per line."""
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def smart_open(filename, mode="rt"):
    """Automatically open gzip or plain text files."""
    return gzip.open(filename, mode) if filename.endswith(".gz") else open(filename, mode)


def filter_record(record, min_qscore=0, min_length=0, read_ids=None):
    avg_q = mean(record.letter_annotations["phred_quality"])
    if len(record.seq) < min_length:
        return ("short", record)
    elif avg_q < min_qscore:
        return ("low_q", record)
    elif (read_ids is None) or (record.id in read_ids):
        return ("passed", record)
    return (None, None)


def process_batch(records, min_qscore, min_length, read_ids):
    results = {"passed": [], "low_q": [], "short": []}
    for record in records:
        status, r = filter_record(record, min_qscore, min_length, read_ids)
        if status:
            results[status].append(r)
    return results


def chunker(iterator, size):
    """Yield chunks of size from an iterator."""
    while True:
        chunk = list(islice(iterator, size))
        if not chunk:
            break
        yield chunk


def extract_from_fastq(input_fastq, output_fastq, read_ids=None, min_qscore=0, min_length=0, threads=None):
    """
    Extract and filter reads from a FASTQ file using multiprocessing.

    Parameters:
        input_fastq: str, input FASTQ or FASTQ.gz file path
        output_fastq: str, output path for passed reads
        read_ids: set or None
        min_qscore: float, minimum average Q-score
        min_length: int, minimum read length
        threads: int or None (uses all CPUs if None)
    """
    n_cpu = threads or cpu_count()
    batch_size = 10000
    passed = []
    low_q = []
    short = []

    with smart_open(input_fastq) as in_f, open(output_fastq, "w") as out_f:
        parser = SeqIO.parse(in_f, "fastq")
        with Pool(processes=n_cpu) as pool:
            jobs = [
                pool.apply_async(process_batch, (batch, min_qscore, min_length, read_ids))
                for batch in chunker(parser, batch_size)
            ]
            for job in jobs:
                result = job.get()
                passed.extend(result["passed"])
                low_q.extend(result["low_q"])
                short.extend(result["short"])
        SeqIO.write(passed, out_f, "fastq")

    return passed, low_q, short


def extract_from_bam(input_bam, output_bam, read_ids=None, min_qscore=0, min_length=0):
    """
    Extract and filter reads from a BAM file.

    Parameters:
        input_bam: str, path to input BAM
        output_bam: str, path to output BAM
        read_ids: set or None
        min_qscore: float
        min_length: int
    """
    import pysam

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
            mean_q = mean(qual) if qual else 0

            if len(seq) < min_length:
                short.append(read)
            elif mean_q < min_qscore:
                low_q.append(read)
            else:
                passed.append(read)
                out_passed.write(read)

    return passed, low_q, short
