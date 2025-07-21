import os
import gzip
import re
from statistics import mean
from itertools import islice
from multiprocessing import Pool, cpu_count

from Bio import SeqIO
import pysam


def read_id_list(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def smart_open(filename, mode="rt"):
    return gzip.open(filename, mode) if filename.endswith(".gz") else open(filename, mode)


def chunker(iterator, size):
    while True:
        chunk = list(islice(iterator, size))
        if not chunk:
            break
        yield chunk


def extract_barcode_from_id(read_id):
    match = re.search(r"(barcode\d+)", read_id)
    return match.group(1) if match else "unknown"


def extract_barcode_from_bam(read):
    try:
        return read.get_tag("CB")
    except KeyError:
        return extract_barcode_from_id(read.query_name) or "unknown"


def filter_record(record, min_qscore=0, min_length=0, read_ids=None):
    avg_q = mean(record.letter_annotations["phred_quality"])
    if len(record.seq) < min_length:
        return "short", record
    elif avg_q < min_qscore:
        return "low_q", record
    elif read_ids is None or record.id in read_ids:
        return "passed", record
    return None, None


def process_batch(records, min_qscore, min_length, read_ids):
    results = {"passed": [], "low_q": [], "short": []}
    for record in records:
        status, r = filter_record(record, min_qscore, min_length, read_ids)
        if status:
            results[status].append(r)
    return results


def extract_from_fastq(input_fastq, output_dir, read_ids=None, min_qscore=0, min_length=0, threads=None, verbose=False, gzip_out=False):
    os.makedirs(output_dir, exist_ok=True)
    n_cpu = threads or cpu_count()
    batch_size = 10000
    total = {"passed": 0, "low_q": 0, "short": 0}
    barcode_bins = {}
    low_q_short_file = os.path.join(output_dir, "low_Q_short.fastq" + (".gz" if gzip_out else ""))
    writer_kwargs = {"format": "fastq"}

    with smart_open(input_fastq) as in_f, smart_open(low_q_short_file, "wt") as lqf:
        parser = SeqIO.parse(in_f, "fastq")
        with Pool(processes=n_cpu) as pool:
            for i, batch in enumerate(chunker(parser, batch_size), 1):
                job = pool.apply_async(process_batch, (batch, min_qscore, min_length, read_ids))
                result = job.get()

                total["passed"] += len(result["passed"])
                total["low_q"] += len(result["low_q"])
                total["short"] += len(result["short"])

                for rec in result["passed"]:
                    barcode = extract_barcode_from_id(rec.id)
                    if barcode not in barcode_bins:
                        barcode_dir = os.path.join(output_dir, barcode)
                        os.makedirs(barcode_dir, exist_ok=True)
                        barcode_file = os.path.join(barcode_dir, f"{barcode}.fastq" + (".gz" if gzip_out else ""))
                        barcode_bins[barcode] = smart_open(barcode_file, "wt")

                    SeqIO.write(rec, barcode_bins[barcode], **writer_kwargs)

                SeqIO.write(result["low_q"] + result["short"], lqf, **writer_kwargs)

                if verbose:
                    print(f"[DEBUG] Batch {i} complete — Passed: {total['passed']}, LowQ: {total['low_q']}, Short: {total['short']}")

        for handle in barcode_bins.values():
            handle.close()

    return total["passed"], total["low_q"], total["short"]


def extract_from_bam(input_bam, output_dir, read_ids=None, min_qscore=0, min_length=0, verbose=False):
    os.makedirs(output_dir, exist_ok=True)
    read_ids = set(read_ids) if read_ids else None
    total = {"passed": 0, "low_q": 0, "short": 0}
    barcode_bins = {}
    low_q_short_file = os.path.join(output_dir, "low_Q_short.bam")

    with pysam.AlignmentFile(input_bam, "rb") as infile, \
         pysam.AlignmentFile(low_q_short_file, "wb", template=infile) as lqf:

        for i, read in enumerate(infile, 1):
            if read_ids and read.query_name not in read_ids:
                continue

            seq = read.query_sequence
            qual = read.query_qualities
            mean_q = mean(qual) if qual else 0

            if len(seq) < min_length:
                total["short"] += 1
                lqf.write(read)
            elif mean_q < min_qscore:
                total["low_q"] += 1
                lqf.write(read)
            else:
                total["passed"] += 1
                barcode = extract_barcode_from_bam(read)
                if barcode not in barcode_bins:
                    barcode_dir = os.path.join(output_dir, barcode)
                    os.makedirs(barcode_dir, exist_ok=True)
                    barcode_file = os.path.join(barcode_dir, f"{barcode}.bam")
                    barcode_bins[barcode] = pysam.AlignmentFile(barcode_file, "wb", template=infile)
                barcode_bins[barcode].write(read)

            if verbose and i % 10000 == 0:
                print(f"[DEBUG] Processed {i} reads: Passed={total['passed']}, LowQ={total['low_q']}, Short={total['short']}")

    for handle in barcode_bins.values():
        handle.close()

    return total["passed"], total["low_q"], total["short"]
