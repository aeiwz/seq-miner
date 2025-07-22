import os
import gzip
import pysam
from Bio import SeqIO
from statistics import mean
from multiprocessing import Pool, cpu_count
from itertools import islice
from collections import defaultdict

def smart_open(filename, mode="rt"):
    return gzip.open(filename, mode) if filename.endswith(".gz") else open(filename, mode)

def read_id_list(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}

def get_barcode(record_id):
    parts = record_id.split("_barcode")
    if len(parts) > 1:
        return f"barcode{parts[-1].split()[0]}"
    return "unknown"

def chunker(iterator, size):
    while True:
        chunk = list(islice(iterator, size))
        if not chunk:
            break
        yield chunk

def filter_record(record, min_qscore=0, min_length=0, read_ids=None):
    try:
        avg_q = mean(record.letter_annotations["phred_quality"])
    except Exception:
        avg_q = 0
    if len(record.seq) < min_length:
        return ("short", record)
    elif avg_q < min_qscore:
        return ("low_q", record)
    elif (read_ids is None) or (record.id in read_ids):
        return ("passed", record)
    return (None, None)

def process_batch(records, min_qscore, min_length, read_ids):
    results = defaultdict(list)
    for record in records:
        status, r = filter_record(record, min_qscore, min_length, read_ids)
        if status:
            results[status].append(r)
    return results

def extract_from_fastq(input_fastq, output_dir, read_ids=None, min_qscore=0, min_length=0, threads=None, log=None, verbose=False):
    os.makedirs(output_dir, exist_ok=True)
    n_cpu = threads or cpu_count()
    batch_size = 10000
    counters = {"total": 0, "passed": 0, "short": 0, "low_q": 0}
    barcode_buckets = defaultdict(list)
    low_q_short = []

    with smart_open(input_fastq) as in_f:
        parser = SeqIO.parse(in_f, "fastq")
        with Pool(processes=n_cpu) as pool:
            jobs = [pool.apply_async(process_batch, (batch, min_qscore, min_length, read_ids))
                    for batch in chunker(parser, batch_size)]
            for job in jobs:
                result = job.get()
                for record in result["passed"]:
                    bc = get_barcode(record.id)
                    barcode_buckets[bc].append(record)
                    counters["passed"] += 1
                for record in result["low_q"] + result["short"]:
                    low_q_short.append(record)
                counters["low_q"] += len(result["low_q"])
                counters["short"] += len(result["short"])
                counters["total"] += sum(len(v) for v in result.values())
                if verbose:
                    print(f"[DEBUG] Processed: {counters['total']} | Passed: {counters['passed']} | Low Q: {counters['low_q']} | Short: {counters['short']}")

    for bc, records in barcode_buckets.items():
        out_path = os.path.join(output_dir, f"{bc}", f"{bc}.fastq")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as f:
            SeqIO.write(records, f, "fastq")

    if low_q_short:
        with open(os.path.join(output_dir, "low_Q_short.fastq"), "w") as f:
            SeqIO.write(low_q_short, f, "fastq")

    if log:
        with open(log, "w") as logf:
            for k, v in counters.items():
                logf.write(f"{k}: {v}\n")

    return barcode_buckets, low_q_short

def extract_from_bam(input_bam, output_dir, read_ids=None, min_qscore=0, min_length=0, log=None, verbose=False):
    os.makedirs(output_dir, exist_ok=True)
    infile = pysam.AlignmentFile(input_bam, "rb")
    counters = {"total": 0, "passed": 0, "short": 0, "low_q": 0}
    barcode_files = {}
    low_q_short = []

    for read in infile:
        if read_ids and read.query_name not in read_ids:
            continue

        counters["total"] += 1
        seq = read.query_sequence
        qual = read.query_qualities
        mean_q = mean(qual) if qual else 0

        if not seq or len(seq) < min_length:
            low_q_short.append(read)
            counters["short"] += 1
        elif mean_q < min_qscore:
            low_q_short.append(read)
            counters["low_q"] += 1
        else:
            barcode = get_barcode(read.query_name)
            if barcode not in barcode_files:
                barcode_dir = os.path.join(output_dir, barcode)
                os.makedirs(barcode_dir, exist_ok=True)
                barcode_files[barcode] = pysam.AlignmentFile(os.path.join(barcode_dir, f"{barcode}.bam"), "wb", template=infile)
            barcode_files[barcode].write(read)
            counters["passed"] += 1

        if verbose and counters["total"] % 1000 == 0:
            print(f"[DEBUG] Processed: {counters['total']} | Passed: {counters['passed']} | Low Q: {counters['low_q']} | Short: {counters['short']}")

    for f in barcode_files.values():
        f.close()

    if low_q_short:
        out_bam = os.path.join(output_dir, "low_Q_short.bam")
        with pysam.AlignmentFile(out_bam, "wb", template=infile) as f:
            for read in low_q_short:
                f.write(read)

    if log:
        with open(log, "w") as logf:
            for k, v in counters.items():
                logf.write(f"{k}: {v}\n")

    return counters
