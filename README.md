# FastContext

![GitHub issues](https://img.shields.io/github/issues/regnveig/FastContext?style=flat-square)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/regnveig/FastContext/sandbox?style=flat-square)
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

<dl>
<dt>-1, --r1</dt>
<dd>

**Required.**<br>
*Format:* String<br>
*Description:*
FastQ input R1 file.
May be uncompressed, gzipped or bzipped.<br>
*Usage:* `-1 input.fastq.gz`

</dd>

<dt>-p, --patterns</dt>
<dd>

**Required.**

*Format:* Plain Javascript Object String (Key-Value).
Names must contain 2-16 small Latin and numeric symbols (a-z, 0-9), sequences must contain more than one symbols ATGC.

*Description:* Patterns to look for.
The order of patterns is the order of search.

*Usage:* `-p '{"first": "CTCAGCGCTGAG", "second": "AAAAAA", "third": "GATC"}'`

</dd>

<dt>-t, --tsv</dt>
<dd>

**Required.**

*Format:* String

*Description:*
Output TSV file.
Contains only general statistics (read structure counts and percentages).

*Usage:* `-t statistics.tsv`

</dd>

<dt>-2, --r2</dt>
<dd>

*Format:* String

*Description:*
FastQ input R2 file.
May be uncompressed, gzipped or bzipped.
If single-end mode, ignore this option.

*Usage:* `-2 input_R2.fastq.gz`

</dd>

<dt>-j, --json</dt>
<dd>

*Format:* String

*Description:*
Output JSON.GZ file (gzipped JSON).
Contains extended statistics on pattern sequences, each read or read pair: read structure, Levenshtein distances (see -l option).

*Usage:* `-j statistics.json.gz`

</dd>

<dt>-k, --kmer-size</dt>
<dd>

*Format:* Non-negative Integer

*Default:* `0`

*Description:*
Max size of unrecognized sequence to be written as K-mer of certain length.

*Usage:* `-k 9`

</dd>

<dt>-u, --unrecognized</dt>
<dd>

*Format:* 2-16 small Latin Chars

*Default:* `unknown`

*Description:*
Long unrecognized sequences replacement.

*Usage:* `-u genome`

</dd>

<dt>-m, --max-reads</dt>
<dd>

*Format:* Non-negative Integer

*Default:* `1000000`

*Description:*
Max reads number to analyze (0 -- no limit).
Notice that read number bigger than recommended may cause memory overflow.

*Usage:* `-m 1000`

</dd>

<dt>-f, --rate-floor</dt>
<dd>

*Format:* Float from 0 to 1

*Default:* `0.001`

*Description:*
Min rate to write read structure into statistics TSV table.

*Usage:* `-f 0.1`

</dd>

<dt>-@, --threads</dt>
<dd>

*Format:* Non-negative integer less than `2 * cpu_count()`

*Default:* `cpu_count()`

*Description:*
Threads number.

*Usage:* `-@ 10`

</dd>

<dt>-d, --dont-check-read-names</dt>
<dd>

*Description:*
Don't check read names.
Use this if you have unusual (non-Illumina) paired read names.
Makes sense only in paired-end mode.

*Usage:* -d

</dd>

<dt>-l, --levenshtein</dt>
<dd>

*Description:*
Calculate patterns Levenshtein distances for each position in read.
Results are written into extended statistics file (JSON.GZ).
Notice that it highly increases the time of processing.

*Usage:* `-l`

</dd>

<dt>-h, --help</dt>
<dd>

*Description:*
Show help message and exit.

*Usage:* -h

</dd>

<dt>`-v`, `--version`</dt>
<dd>

*Description:*
Show program's version number and exit.

*Usage:* -v

</dd>
</dl>

## Examples

### General statistics TSV file

Contains counts, percentage and read structures.
Length of K-mer or pattern strand (**F**orward or **R**everse) is displayed after the comma.

	# R1
	Count	Percentage	ReadStructure
	5197	48.807287753568744	{unknown}
	3297	30.963561232156273	{unknown}--{oligme:F}--{oligb:F}--{701:F}--{unknown}
	114	1.0706235912847484	{unknown}--{oligb:F}--{701:F}--{unknown}
	71	0.6667918858001503	{unknown}--{oligme:F}--{unknown}
	69	0.6480090157776108	{unknown}--{oligme:F}--{unknown}--{701:F}--{unknown}
	60	0.5634861006761833	{unknown}--{oligme:F}--{oligb:F}--{701:F}--{kmer:14bp}
	
	# R2
	Count	Percentage	ReadStructure
	7545	70.85837716003005	{unknown}
	616	5.785123966942149	{unknown}--{oligme:F}--{oliga:R}--{502:R}--{unknown}
	540	5.07137490608565	{unknown}--{oligme:F}--{unknown}
	441	4.141622839969948	{unknown}--{oligme:F}--{oliga:R}--{unknown}
	298	2.798647633358377	{unknown}--{oliga:R}--{unknown}
	263	2.469947407963937	{unknown}--{502:R}--{unknown}
	233	2.188204357625845	{unknown}--{oligme:F}--{kmer:14bp}--{502:R}--{unknown}
	163	1.5308039068369648	{unknown}--{oliga:R}--{502:R}--{unknown}
	56	0.5259203606311045	{unknown}--{502:F}--{unknown}

### Extended statistics JSON.GZ file

Contains extended statistics: run options, performance, pattern analysis, full summary without rate floor, each read analysis.
Example is shorten.

```json
{
	"FastQ": {
		"R1": "tests/standard_test_R1.fastq.gz",
		"R2": "tests/standard_test_R2.fastq.gz"
		},
	"ReadType": "Paired-end",
	"ReadsAnalyzed": 100,
	"MaxReads": 100,
	"Performance": {
		"Threads": 12,
		"Started": "2022-07-07T20:30:11.503120",
		"Finished": "2022-07-07T20:30:11.680769" },
	"Patterns": {
		"oligme": { "F": "CTGTCTCTTATACACATCT", "R": "AGATGTGTATAAGAGACAG", "Length": 19},
		"502": {"F": "CTCTCTAT", "R": "ATAGAGAG", "Length": 8 }
	},
	"KmerMaxSize": 15,
	"RateFloor": 0.001,
	"UnrecognizedSequence": "unknown",
	"PatternsLevenshtein": [
		{
			"Analysis": "reverse complement only",
			"FirstPattern": "oligme",
			"SecondPattern": "oligme",
			"FirstLength": 19,
			"SecondLength": 19,
			"LevenshteinAbsolute": 11,
			"LevenshteinSimilarity": 0.42105263157894735,
			"Type":"good",
			"Risk":"low"
		},
		{
			"Analysis": "full",
			"FirstPattern": "oligme",
			"SecondPattern": "502",
			"FirstLength": 19,
			"SecondLength": 8,
			"LevenshteinAbsolute": 2,
			"LevenshteinSimilarity": 0.75,
			"Type": "nested",
			"Risk":"medium"
		}
	],
	"Summary": {
		"R1": {
			"{unknown}--{oligme:F}--{oligb:F}--{701:F}--{unknown}": { "Count": 34, "Rate": 0.34 },
			"{unknown}--{oligme:F}--{unknown}--{701:F}--{unknown}": {"Count": 1, "Rate": 0.01 },
			"{unknown}--{701:F}--{unknown}": { "Count": 1, "Rate" :0.01 }
		},
		"R2": {
			"{unknown}": {"Count": 83, "Rate": 0.83},
			"{unknown}--{502:R}--{unknown}": { "Count": 4, "Rate": 0.04},
			"{unknown}--{oligme:F}": { "Count": 1, "Rate": 0.01}
		}
	},
	"RawDataset": [
		{
			"Name": "M02435:112:000000000-DFC9M:1:1101:16174:1450",
			"R1": {
				"Sequence": "GTTTAGCTGCAAGCAGGTTTTGTTTTGTGTTAGGGA...TCGTATGCCGTCTTCTGCTTGAAAAAAAAAA",
				"PhredQual": [18,29,32,32,33,37,33,37,37,33,37,37,"...",38,38,37,37,13,12,24,24,12,12],
				"Levenshtein":[
					{
						"name": "oligme",
						"strand": "F",
						"length": 19,
						"values": [11,11,12,13,14,14,13,14,15,13,14,13,13,"...",NaN,NaN,NaN,NaN,NaN]
					},
					{
						"name": "oligme",
						"strand": "R",
						"length": 19,
						"values":[12,13,13,12,12,14,16,15,14,14,13,13,13,"...",NaN,NaN,NaN,NaN,NaN]
					}
				],
				"ReadStructure": [
					{ "type": "unrecognized" },
					{ "type": "pattern", "name": "oligme", "strand": "F" },
					{ "type": "pattern", "name": "oligb", "strand": "F" },
					{ "type": "pattern", "name": "701", "strand": "F" },
					{ "type": "unrecognized" }
				],
				"TextReadStructure": "{unknown}--{oligme:F}--{oligb:F}--{701:F}--{unknown}"
			},
			"R2": "..."
		}
	]
}
```
