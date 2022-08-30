# FastContext

![Logo](https://github.com/regnveig/FastContext/blob/stable/logo.png)

![PyPI](https://img.shields.io/pypi/v/FastContext?style=flat-square)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/FastContext?style=flat-square)
![PyPI - Status](https://img.shields.io/pypi/status/FastContext?style=flat-square)
![PyPI - Downloads](https://img.shields.io/pypi/dm/FastContext?style=flat-square)
![GitHub issues](https://img.shields.io/github/issues/regnveig/FastContext?style=flat-square)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/regnveig/FastContext/sandbox?style=flat-square)
![GitHub](https://img.shields.io/github/license/regnveig/FastContext?style=flat-square)
![Keybase PGP](https://img.shields.io/keybase/pgp/regnveig?style=flat-square)

## Description

FastContext is a tool for identification of adapters and other sequence patterns in the next generation sequencing (NGS) data.
The algorithm parses FastQ files (in a single-end or paired-end mode), searches read or read pair for user-specified patterns, and then generates a human-readable representation of the search results, which we call "read structure".
Also FastContext gathers statistics on frequency of occurence for each read structure.

## Installation

```bash
python3 -m pip install FastContext
```

Check installation:

```bash
FastContext --help
```

## Usage

Optional arguments:

<dl>
<dt>-1, --r1</dt>
<dd>

*Required.* <br>
*Format:* String<br>
*Description:*
FastQ input R1 file.
May be uncompressed, gzipped or bzipped.<br>
*Usage:* `-1 input.fastq.gz`

</dd>

<dt>-p, --patterns</dt>
<dd>

*Required.* <br>
*Description:* Patterns to look for.
The order of patterns is the order of search.<br>
*Format:* Plain Javascript Object String (Key-Value).
Names must contain 2-24 Latin and numeric symbols, and -_-, sequences must contain more than one symbols ATGC.<br>
*Usage:* `-p '{"First": "CTCAGCGCTGAG", "Second": "AAAAAA", "Third": "GATC"}'`

</dd>

<dt>-s, --summary</dt>
<dd>

*Required.* <br>
*Description:*
Output HTML file.
Contains statistics summary in human-readable form.<br>
*Format:* String<br>
*Usage:* `-s statistics.htm`

</dd>

<dt>-2, --r2</dt>
<dd>

*Description:*
FastQ input R2 file.
May be uncompressed, gzipped or bzipped.
If single-end mode, ignore this option.<br>
*Format:* String<br>
*Usage:* `-2 input_R2.fastq.gz`

</dd>

<dt>-j, --json</dt>
<dd>

*Description:*
Output JSON.GZ file (gzipped JSON).
Contains extended statistics on pattern sequences, each read or read pair: read structure, Levenshtein distances (see -l option).<br>
*Format:* String<br>
*Usage:* `-j statistics.json.gz`

</dd>

<dt>-k, --kmer-size</dt>
<dd>

*Description:*
Max size of unrecognized sequence to be written as K-mer of certain length.<br>
*Format:* Non-negative Integer<br>
*Default:* `0` <br>
*Usage:* `-k 9`

</dd>

<dt>-u, --unrecognized</dt>
<dd>

*Description:*
Long unrecognized sequences replacement.<br>
*Format:* 2-24 Latin and numeric symbols, and -_-<br>
*Default:* `unknown` <br>
*Usage:* `-u genome`

</dd>

<dt>-m, --max-reads</dt>
<dd>

*Description:*
Max reads number to analyze (0 -- no limit).
Notice that read number bigger than recommended may cause memory overflow.<br>
*Format:* Non-negative Integer<br>
*Default:* `1000000` <br>
*Usage:* `-m 1000`

</dd>

<dt>-f, --rate-floor</dt>
<dd>

*Description:*
Min rate to write read structure into statistics TSV table.<br>
*Format:* Float from 0 to 1<br>
*Default:* `0.001` <br>
*Usage:* `-f 0.1`

</dd>

<dt>-@, --threads</dt>
<dd>

*Description:*
Threads number.<br>
*Format:* Non-negative integer less than `2 * cpu_count()` <br>
*Default:* `cpu_count()` <br>
*Usage:* `-@ 10`

</dd>

<dt>-d, --dont-check-read-names</dt>
<dd>

*Description:*
Don't check read names.
Use this if you have unusual (non-Illumina) paired read names.
Makes sense only in paired-end mode.<br>
*Usage:* `-d`

</dd>

<dt>-l, --levenshtein</dt>
<dd>

*Description:*
Calculate patterns Levenshtein distances for each position in read.
Results are written into extended statistics file (JSON.GZ).
Notice that it highly increases the time of processing.<br>
*Usage:* `-l`

</dd>

<dt>-h, --help</dt>
<dd>

*Description:*
Show help message and exit.<br>
*Usage:* `-h`

</dd>

<dt>-v, --version</dt>
<dd>

*Description:*
Show program's version number and exit.<br>
*Usage:* `-v`

</dd>
</dl>

## Examples

### Summary statistics table

Contains counts, percentage and read structures.
Length of K-mer or pattern strand (**F**orward or **R**everse) is displayed after the comma.

Example:

#### R1

| Count | Percentage | Read Structure                                         |
|------:|-----------:|:-------------------------------------------------------|
| 5,197 |     48.807 | {unknown}                                              |
| 3,297 |     30.963 | {unknown}--{oligme:F}--{oligb:F}--{701:F}--{unknown}   |
|   114 |      1.070 | {unknown}--{oligb:F}--{701:F}--{unknown}               |
|    71 |      0.666 | {unknown}--{oligme:F}--{unknown}                       |
|    69 |      0.648 | {unknown}--{oligme:F}--{unknown}--{701:F}--{unknown}   |
|    60 |      0.563 | {unknown}--{oligme:F}--{oligb:F}--{701:F}--{kmer:14bp} |

#### R2

| Count | Percentage | Read Structure                                         |
|------:|-----------:|:-------------------------------------------------------|
| 7,545 |     70.858 | {unknown}                                              |
|   616 |      5.785 | {unknown}--{oligme:F}--{oliga:R}--{502:R}--{unknown}   |
|   540 |      5.071 | {unknown}--{oligme:F}--{unknown}                       |
|   441 |      4.141 | {unknown}--{oligme:F}--{oliga:R}--{unknown}            |
|   298 |      2.798 | {unknown}--{oliga:R}--{unknown}                        |
|   263 |      2.469 | {unknown}--{502:R}--{unknown}                          |
|   233 |      2.188 | {unknown}--{oligme:F}--{kmer:14bp}--{502:R}--{unknown} |
|   163 |      1.530 | {unknown}--{oliga:R}--{502:R}--{unknown}               |
|    56 |      0.525 | {unknown}--{502:F}--{unknown}                          |

### Extended statistics JSON.GZ file

Contains extended statistics: run options, performance, pattern analysis, full summary without rate floor, each read analysis.
Example is shorten.

```json
{
	"FastQ": {
		"R1": "tests/standard_test_R1.fastq.gz",
		"R2": "tests/standard_test_R2.fastq.gz"
	},
	"RunData": {
		"Read Type": "Paired-end",
		"Max Reads": 100,
		"Rate Floor": 0.001
	},
	"Performance": {
		"Reads Analyzed": 100,
		"Threads": 4,
		"Started": "2022-07-13T18:15:48.277660",
		"Finished": "2022-07-13T18:15:48.964721"
	},
	"PatternsData": {
		"PatternsList": {
			"oligme": {
				"F": "CTGTCTCTTATACACATCT",
				"R": "AGATGTGTATAAGAGACAG",
				"Length": 19
			},
			"s502": {
				"F": "CTCTCTAT",
				"R": "ATAGAGAG",
				"Length": 8
			}
		},
		"PatternsAnalysis": [
			{
				"Analysis": "reverse complement only",
				"FirstPattern": "oligme",
				"SecondPattern": "oligme",
				"FirstLength": 19,
				"SecondLength": 19,
				"LevenshteinAbsolute": 11,
				"LevenshteinSimilarity": 0.42105263157894735,
				"Type": "good",
				"Risk": "low"
			},
			{
				"Analysis": "full",
				"FirstPattern": "oligme",
				"SecondPattern": "s502",
				"FirstLength": 19,
				"SecondLength": 8,
				"LevenshteinAbsolute": 2,
				"LevenshteinSimilarity": 0.75,
				"Type": "nested",
				"Risk": "medium"
			}
		],
		"Other": {
			"Unrecognized Sequence": "unknown",
			"K-mer Max Size": 15
		}
	},
	"Summary": {
		"R1": {
			"{unknown}--{oligme:F}--{oligb:F}--{s701:F}--{unknown}": {
				"Count": 34,
				"Percentage": 34.0,
				"ReadStructure": [
					{ "type": "unrecognized" },
					{ "type": "pattern", "name": "oligme", "strand": "F" },
					{ "type": "pattern", "name": "oligb", "strand": "F" },
					{ "type": "pattern", "name": "s701", "strand": "F" },
					{ "type": "unrecognized" }
				]
			},
			"{unknown}--{oligme:F}--{unknown}--{s701:F}--{unknown}": [ "..." ],
			"{unknown}--{s701:F}--{unknown}": [ "..." ]
		},
		"R2": [ "..." ]
	},
	"RawDataset": [
		{
			"Name": "M02435:112:000000000-DFC9M:1:1101:14970:1484",
			"R1": {
				"Sequence": "ACCTAGAAGAGCCAAAAGACTCT...AATCTCGTATGCCGTCT",
				"PhredQual": [29,32,32,33,33,37,37,37,37,"...",38,38,38,13],
				"Levenshtein": [
					{
						"name": "oligme",
						"strand": "F",
						"length": 19,
						"values": [14,14,12,13,12,12,12,"...",NaN,NaN,NaN]
					},
					{
						"name": "oligme",
						"strand": "R",
						"length": 19,
						"values": [12,11,10,9,9,9,10,10,"...",NaN,NaN,NaN]
					}
				],
				"ReadStructure": [
					{ "type": "unrecognized" },
					{ "type": "pattern", "name": "oligme", "strand": "F" },
					{ "type": "pattern", "name": "oligb", "strand": "F" },
					{ "type": "pattern", "name": "s701", "strand": "F" },
					{ "type": "unrecognized" }
				],
				"TextReadStructure": "{unknown}--{oligme:F}--{oligb:F}--{s701:F}--{unknown}"
			},
			"R2": "..." 
		}
	]
}
```
