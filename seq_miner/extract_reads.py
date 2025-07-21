import os
import gzip
from Bio import SeqIO
from statistics import mean
from multiprocessing import Pool, cpu_count
from itertools import islice
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


def extract_barcode_from_name(name):
    for part in name.split():
        if part.startswith("barcode"):
            return part
    return None


def extract_barcode_from_read(read):
    try:
        return read.get_tag("CB")
    except Exception:
        return extract_barcode_from_name(read.query_name)


def extract_from_fastq(input_fastq, output_dir, read_ids=None, min_qscore=0, min_length=0, threads=None, verbose=False, gzip_out=False):
    n_cpu = threads or cpu_count()
    batch_size = 10000
    passed_by_barcode = {}
    rejected_reads = []
    count = {"total": 0, "passed": 0, "low_q": 0, "short": 0}

    def filter_record(record):
        count["total"] += 1
        avg_q = mean(record.letter_annotations["phred_quality"])
        if len(record.seq) < min_length:
            count["short"] += 1
            return "short", record
        elif avg_q < min_qscore:
            count["low_q"] += 1
            return "low_q", record
        elif (read_ids is None or record.id in read_ids):
            count["passed"] += 1
            return "passed", record
        return None, None

    with smart_open(input_fastq) as in_f:
        parser = SeqIO.parse(in_f, "fastq")
        for batch in chunker(parser, batch_size):
            for record in batch:
                status, r = filter_record(record)
                if status == "passed":
                    bc = extract_barcode_from_name(r.id) or "unknown"
                    passed_by_barcode.setdefault(bc, []).append(r)
                elif status in {"low_q", "short"}:
                    rejected_reads.append(r)

            if verbose:
                print(f"[INFO] Processed {count['total']} reads "
                      f"(Passed: {count['passed']}, Low-Q: {count['low_q']}, Short: {count['short']})")

    # Write passed by barcode
    for bc, recs in passed_by_barcode.items():
        bc_dir = os.path.join(output_dir, bc)
        os.makedirs(bc_dir, exist_ok=True)
        bc_file = os.path.join(bc_dir, f"{bc}.fastq.gz" if gzip_out else f"{bc}.fastq")
        with smart_open(bc_file, "wt") as f:
            SeqIO.write(recs, f, "fastq")

    # Write unclassified
    rej_file = os.path.join(output_dir, "unclassified.fastq.gz" if gzip_out else "unclassified.fastq")
    with smart_open(rej_file, "wt") as f:
        SeqIO.write(rejected_reads, f, "fastq")

    return count


def extract_from_bam(input_bam, output_dir, read_ids=None, min_qscore=0, min_length=0, verbose=False):
    passed_by_barcode = {}
    rejected_reads = []
    count = {"total": 0, "passed": 0, "low_q": 0, "short": 0}

    infile = pysam.AlignmentFile(input_bam, "rb")
    bam_writers = {}

    # Writer for rejected
    os.makedirs(output_dir, exist_ok=True)
    rejected_path = os.path.join(output_dir, "unclassified.bam")
    rejected_out = pysam.AlignmentFile(rejected_path, "wb", template=infile)

    for read in infile:
        count["total"] += 1
        if read_ids and read.query_name not in read_ids:
            continue
        if read.is_unmapped:
            continue
        seq = read.query_sequence
        qual = read.query_qualities
        mean_q = mean(qual) if qual else 0

        if len(seq) < min_length:
            count["short"] += 1
            rejected_out.write(read)
        elif mean_q < min_qscore:
            count["low_q"] += 1
            rejected_out.write(read)
        else:
            count["passed"] += 1
            bc = extract_barcode_from_read(read) or "unknown"
            bc_dir = os.path.join(output_dir, bc)
            os.makedirs(bc_dir, exist_ok=True)
            bc_file = os.path.join(bc_dir, f"{bc}.bam")

            if bc not in bam_writers:
                bam_writers[bc] = pysam.AlignmentFile(bc_file, "wb", template=infile)
            bam_writers[bc].write(read)

        if verbose and count["total"] % 5000 == 0:
            print(f"[INFO] Processed {count['total']} reads "
                  f"(Passed: {count['passed']}, Low-Q: {count['low_q']}, Short: {count['short']})")

    for writer in bam_writers.values():
        writer.close()
    rejected_out.close()
    infile.close()
    return count
