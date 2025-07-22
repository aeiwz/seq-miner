from setuptools import setup, find_packages

setup(
    name="seq-miner",
    version="1.3.1",
    description="Extract and filter reads from BAM/FASTQ by ID, quality, and length.",
    author="Aeiwz",
    packages=find_packages(),
    install_requires=[
        "pysam",
        "biopython"
    ],
    entry_points={
        "console_scripts": [
            "seq-miner=seq_miner.cli:main"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
    ],
    python_requires='>=3.7',
)
