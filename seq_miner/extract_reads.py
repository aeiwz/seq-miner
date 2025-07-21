from Bio import SeqIO
from statistics import mean
from multiprocessing import Pool, cpu_count, Value, Lock
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

def smart_write(filename):
    return gzip.open(filename, "wt") if filename.endswith(".gz") else open(filename, "w")

def parse_barcode(read_id, key="barcode"):
    # Example: read_id = 'abc_barcode01_xyz', extract 'barcode01'
    for part in read_id.split("_"):
        if key in part:
            return part
    return "unknown"

def filter_record(record, min_qscore=0, min_length=0, read_ids=None):
    avg_q = mean(record.letter_annotations["phred_quality"])
    if len(record.seq) < min_length:
        return ("short", record)
    elif avg_q < min_qscore:
        return ("low_q", record)
    elif (read_ids is None) or (record.id in read_ids):
        return ("passed", record)
    return (None, None)

def process_batch(records, min_qscore, min_length, read_ids, barcode_key="barcode", counter=None, lock=None, verbose=False):
    results = {"passed": {}, "low_q": [], "short": []}

    for record in records:
        status, r = filter_record(record, min_qscore, min_length, read_ids)
        if status == "passed":
            bc = parse_barcode(record.id, barcode_key)
            results["passed"].setdefault(bc, []).append(r)
        elif status == "low_q":
            results["low_q"].append(r)
        elif status == "short":
            results["short"].append(r)

    if counter is not None:
        with lock:
            counter.value += len(records)
            if verbose and counter.value % 10000 == 0:
                print(f"[INFO] Processed {counter.value} reads...")

    return results

def chunker(iterator, size):
    """Yield chunks of size from an iterator."""
    while True:
        chunk = list(islice(iterator, size))
        if not chunk:
            break
        yield chunk

def extract_from_fastq(input_fastq, output_dir, read_ids=None, min_qscore=0, min_length=0,
                       threads=None, gzip_output=False, barcode_key="barcode", verbose=False):
    n_cpu = threads or cpu_count()
    batch_size = 10000
    passed_reads = {}
    low_q = []
    short = []

    os.makedirs(output_dir, exist_ok=True)
    unclassified_file = os.path.join(output_dir, "unclassified.fastq.gz" if gzip_output else "unclassified.fastq")

    progress_counter = Value("i", 0)
    counter_lock = Lock()

    if verbose:
        print(f"[INFO] Using {n_cpu} CPUs to process FASTQ in batches of {batch_size}")

    with smart_open(input_fastq) as in_f:
        parser = SeqIO.parse(in_f, "fastq")
        with Pool(processes=n_cpu) as pool:
            jobs = [
                pool.apply_async(
                    process_batch,
                    (batch, min_qscore, min_length, read_ids),
                    {"barcode_key": barcode_key, "counter": progress_counter, "lock": counter_lock, "verbose": verbose}
                )
                for batch in chunker(parser, batch_size)
            ]
            for job in jobs:
                result = job.get()
                for bc, recs in result["passed"].items():
                    passed_reads.setdefault(bc, []).extend(recs)
                low_q.extend(result["low_q"])
                short.extend(result["short"])

    for bc, recs in passed_reads.items():
        bc_dir = os.path.join(output_dir, bc)
        os.makedirs(bc_dir, exist_ok=True)
        bc_file = os.path.join(bc_dir, f"{bc}.fastq.gz" if gzip_output else f"{bc}.fastq")
        with smart_write(bc_file) as out_f:
            SeqIO.write(recs, out_f, "fastq")

    with smart_write(unclassified_file) as out_f:
        SeqIO.write(low_q + short, out_f, "fastq")

    if verbose:
        total = progress_counter.value
        n_passed = sum(len(v) for v in passed_reads.values())
        print(f"[INFO] Total reads processed  : {total}")
        print(f"[INFO] Passed reads           : {n_passed}")
        print(f"[INFO] Low Q-score reads      : {len(low_q)}")
        print(f"[INFO] Short length reads     : {len(short)}")

    return passed_reads, low_q, short
