# seq-miner

**seq-miner** is a fast, Python-based command-line tool to extract and filter sequence reads from BAM and FASTQ files by:

- Read ID (single or batch)
- Mean quality score
- Minimum read length

Built for researchers working in genomics, transcriptomics, and metagenomics.


## Installation

Install via `pip`:

```bash
pip install seq-miner
```


## Usage

```bash
seq-miner --input INPUT --output OUTPUT --format FORMAT [options]
```

### Example 1: Filter FASTQ reads by Q-score and length

```bash
seq-miner -i sample.fastq -o filtered.fastq -f fastq --min_qscore 15 --min_length 100
```

### Example 2: Extract specific read IDs from a BAM file

```bash
seq-miner -i reads.bam -o matched.bam -f bam -r read_ids.txt --min_qscore 10 --min_length 200
```


## Options

| Flag               | Description                                                  |
|--------------------|--------------------------------------------------------------|
| `-i`, `--input`     | Input BAM or FASTQ file                                      |
| `-o`, `--output`    | Output file to write filtered/passed reads                  |
| `-f`, `--format`    | File format: `bam` or `fastq`                                |
| `-r`, `--read_ids`  | File containing read IDs (one per line, optional)            |
| `--min_qscore`      | Minimum average quality score per read (default: `0`)        |
| `--min_length`      | Minimum length per read (default: `0`)                       |


## Input Examples

### FASTQ file (`.fastq`)

Supports gzipped or plain FASTQ format.

### BAM file (`.bam`)

Requires [pysam](https://github.com/pysam-developers/pysam) under the hood.

### Read ID file (optional)

```txt
read00001
read00044
read20398
```


## Output

- Filtered reads saved to the specified output file.
- CLI prints counts of:
  - Passed reads
  - Low-quality reads
  - Short reads


## Dependencies

- [Biopython](https://biopython.org/)
- [pysam](https://github.com/pysam-developers/pysam)

Installable automatically via `pip install seq-miner`.


## Publishing (for maintainers)

To publish:

```bash
python -m build
twine upload dist/*
```

Or use GitHub Actions (see `.github/workflows/pypi-release.yml`) for trusted publishing.


## License

MIT License © Theerayut  
See [LICENSE](LICENSE) for full text.


## Contact

For issues, please open an issue on [GitHub](https://github.com/aeiwz/seq-miner/issues).
