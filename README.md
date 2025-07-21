# seq-miner

**seq-miner** is a lightweight, fast, and parallelizable tool to extract and filter reads from **BAM** or **FASTQ** files based on:

- Specific read IDs
- Mean quality score threshold
- Minimum read length
- Multi-threading (FASTQ)
- JSON/CSV-ready summary (optional)
- GitHub release tagging and PyPI publish automation


## Installation

```bash
pip install seq-miner
```

Or clone from source:

```bash
git clone https://github.com/your-org/seq-miner.git
cd seq-miner
pip install .
```


## Usage

### Extract reads from BAM

```bash
seq-miner -i reads.bam -o filtered.bam -f bam -r read_ids.txt --min-qscore 10 --min-length 200
```

### Filter FASTQ reads in parallel

```bash
seq-miner -i reads.fastq -o filtered.fastq -f fastq --min-qscore 15 --min-length 1000 --threads 4
```

### Show version

```bash
seq-miner --version
```


## Command-line Options

| Option            | Description                                          |
|-------------------|------------------------------------------------------|
| `-i`, `--input`   | Input BAM or FASTQ file                              |
| `-o`, `--output`  | Output file for passed reads                         |
| `-f`, `--format`  | File format: `bam` or `fastq`                        |
| `-r`, `--read-ids`| Optional file with read IDs (one per line)          |
| `--min-qscore`    | Minimum mean Q-score (default: `0.0`)               |
| `--min-length`    | Minimum read length (default: `0`)                  |
| `--threads`       | Number of CPU threads (only used for FASTQ)         |
| `--verbose`       | Enable verbose logging                               |
| `--version`       | Print the current version                            |


## Output Summary

When finished, you'll see:

```
Summary:
Passed reads     : 12345
Low-quality reads : 54
Short reads      : 91
```

Optionally, you can pipe this to JSON or CSV (coming soon).


## Auto Version + Release

- Version is stored in [`seqminer/__version__.py`](seqminer/__version__.py)
- Tagged automatically with GitHub Actions on push to `main`
- Published to PyPI on GitHub release

## License

MIT License © Theerayut  
See [LICENSE](LICENSE) for full text.


## Contact

For issues, please open an issue on [GitHub](https://github.com/aeiwz/seq-miner/issues).
