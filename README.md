# FastContext

![GitHub issues](https://img.shields.io/github/issues/regnveig/FastContext?style=flat-square)
![GitHub](https://img.shields.io/github/license/regnveig/FastContext?style=flat-square)

## Description

FastQ Nucleotide Context Statistics Calculation

## Dependencies

1. Python 3.8.10+
2. Python packages:
	- argparse 1.4.0+
	- BioPython 1.79+
	- pandas 1.2.5+
	- python-Levenshtein 0.12.2
	- tqdm 4.61.2+

Install dependencies with the one-liner (Linux):

```bash
python3 -m pip install argparse biopython pandas python-Levenshtein tqdm
```

## Usage

Optional arguments:

| Short tag | Long tag | Required | Default | Description |
|:---|:---|:---|
| `-1 INPUTFILE_R1` | `--r1 INPUTFILE_R1` | Yes | -- | FastQ input R1 file (may be gzipped or bzipped) |
| `-2 INPUTFILE_R2` | `--r2 INPUTFILE_R2` | Yes | -- | FastQ input R2 file (may be gzipped or bzipped) |
| `-p PATTERNS` | `--patterns PATTERNS` | Yes | -- | Patterns to look for, plain JSON format: `{"first": "GATC", "second": "CTCAGCGCTGAG"}`. Names must contain 2-16 small Latin symbols (a-z), sequences must contain more than one symbols ATGC. Order of patterns is order of search. |
| `-t OUTPUTTSV` | `--tsv OUTPUTTSV` | Yes | -- | Output TSV file. Contains only general statistics. |
| `-j OUTPUTJSON` | `--json OUTPUTJSON` | No | -- | Output JSON.GZ file (gzipped JSON, raw data). Useful if you need to see particular read structure |
| `-k KMERMAXSIZE` | `--kmer-size KMERMAXSIZE` | No | 0 | Max K-mer size | 
| `-u UNRECOGNIZEDSEQ` | `--unrecognized UNRECOGNIZEDSEQ` | No | "unknown" | Long unrecognized sequences replacement |
| `-m MAXREADS` | `--max-reads MAXREADS` | No | 1000000 | Max reads number to analyze (0 -- no limit)|
| `-f RATEFLOOR` | `--rate-floor RATEFLOOR` | No | 0.001 | Min rate to write pattern into stats TSV table |
| `-@ THREADSNUM` | `--threads THREADSNUM` | No | [cpu_count()](https://docs.python.org/3/library/multiprocessing.html#multiprocessing.cpu_count) | Threads number |
| `-d` | `--dont-check-read-names` | No | False | Don't check read names. Use this if you have unusual paired read names. Makes sense only in paired reads mode. |
| `-l` | `--levenshtein` | No | False | Calculate Levenshtein distances for each position in read. Notice that it highly increases the time of processing. |
| `-h` | `--help` | No | -- | Show help message and exit |
| `-v` | `--version` | No | -- | Show program's version number and exit |
