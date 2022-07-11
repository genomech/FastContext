#!/bin/python3

__version__ = "0.9b"
__author__ = "Emil Viesná"

from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from Levenshtein import distance
import contextlib
import multiprocessing
import argparse
import bz2
import datetime
import functools
import gzip
import itertools
import json
import os
import re
import pandas
import string
import sys
import time
import tqdm
import warnings

# Make one-lined warnings and errors without traceback
# NOTE If debugging, comment the lines below
def OneLinerWarning(message, category, filename, lineno, file = None, line = None): return '%s: %s\n' % (category.__name__, message)
warnings.formatwarning = OneLinerWarning
sys.tracebacklimit = 0

# Default variables
GLOBAL_KMER_SIZE = 0
GLOBAL_UNRECOGNIZED = 'unknown'
GLOBAL_MAX = 1000000
GLOBAL_THREADS = multiprocessing.cpu_count()
GLOBAL_RATE_FLOOR = 0.001

## ------======| THREADING |======------

@contextlib.contextmanager
def Threading(Threads = multiprocessing.cpu_count()):
	pool = multiprocessing.Pool(Threads)
	yield pool
	pool.close()
	pool.join()
	del pool

## ------======| PATTERNS CHECK FUNC |======------

# ReverseComplement: Make reverse complement of DNA sequence
def ReverseComplement(Sequence): return Seq(Sequence).reverse_complement().__str__()

# Check if list is length-sorted
def LengthSorted(List): return all(len(List[index]) >= len(List[index + 1]) for index in range(len(List) - 1))

# Primary check of patterns JSON object
def PatternsCheck(Patterns):
	# Check if line is a valid JSON
	try:
		PatternsDict = json.loads(Patterns)
	except json.decoder.JSONDecodeError as Error:
		raise ValueError(f'Patterns string is not a valid JSON string: {Patterns}\n{Error}')
	# Check if JSON contains non-empty dictionary
	if type(PatternsDict) != dict: raise ValueError(f'Patterns string is not a JavaScript Object (dictionary): {Patterns}')
	if len(PatternsDict) < 1: raise ValueError(f'Patterns string does not contain any patterns: {Patterns}')
	# Check every key-value
	for Key, Value in PatternsDict.items():
		# Check if key and value are strings
		if (type(Key) != str): raise ValueError(f'Invalid type of pattern name ({type(Key)}), must be string: {Key}')
		if (type(Value) != str): raise ValueError(f'Invalid type of pattern sequence ({type(Value)}), must be string: {Value}')
		# Check if key contains only small Latin and numerics
		# It is important to avoid unwanted matches with ATGC
		if re.fullmatch('^[a-z0-9]{2,16}$', Key) is None: raise ValueError(f'Pattern name must contain 2-16 small Latin and numeric symbols: {Key}')
		# Check if value contains valid ATGC sequence
		if re.fullmatch('^[ATGC]+$', Value) is None: raise ValueError(f'Pattern sequence must contain more than one symbols ATGC: {Value}')
		# Warn if sequence is too long
		if len(Value) > 16: warnings.warn(f'Pattern sequence "{Key}" is very long ({len(Value)} bp). This may affect search. You may want to use shorter fragments of it.', RuntimeWarning)
	# Warn if patterns are not length-sorted
	if not LengthSorted(list(PatternsDict.values())): warnings.warn(f'Patterns sequences are not sorted from longest to shortest. Notice that order of patterns matters in search.', RuntimeWarning)

# Check Levenshtein distance between pattern sequences, including their reverse complements
def PatternsLevenshtein(Patterns):
	# Levenshtein distance to risk rank
	SimilarityType = { 0: 'full', 1: 'high', 2: 'medium' }
	Result = []
	# Check every possible pair of patterns
	for KeyFirst, KeySecond in itertools.combinations_with_replacement(Patterns.keys(), 2):
		# Create object with reverse complements
		BlockFR = {
			'F': {
				'First': str(Patterns[KeyFirst]),
				'Second': str(Patterns[KeySecond])
			},
			'R': {
				'First': ReverseComplement(Patterns[KeyFirst]),
				'Second': ReverseComplement(Patterns[KeySecond])
			}
		}
		# Set distance to infinity until calculated
		LevenshteinDistAbs = float('+inf')
		# Set longer and shorter pattern sequence
		Longer, Shorter = ('First', 'Second') if len(BlockFR['F']['First']) >= len(BlockFR['F']['Second']) else ('Second', 'First')
		
		# Calculate the minimum of Levenshtein distance
		for Second in ['F', 'R']:
			for FragmentStart in range(len(BlockFR['F'][Longer]) - len(BlockFR[Second][Shorter]) + 1):
				Dist = distance(BlockFR['F'][Longer][FragmentStart : FragmentStart + len(BlockFR[Second][Shorter])], BlockFR[Second][Shorter])
				# If less than minimum, set minimum
				# Check if we're not comparing a sequence with itself 
				if (Dist < LevenshteinDistAbs) and not ((KeyFirst == KeySecond) and (Second == 'F')): LevenshteinDistAbs = int(Dist)
		# Create report
		Data = {
			# Analysis: full or reverse complement only
			'Analysis': 'full' if KeyFirst != KeySecond else 'reverse complement only',
			'FirstPattern': str(KeyFirst),
			'SecondPattern': str(KeySecond),
			'FirstLength': len(Patterns[KeyFirst]),
			'SecondLength': len(Patterns[KeySecond]),
			'LevenshteinAbsolute': int(LevenshteinDistAbs),
			# Similarity
			'LevenshteinSimilarity': float((len(BlockFR['F'][Shorter]) - LevenshteinDistAbs) / len(BlockFR['F'][Shorter]))
			}
		# Set type of pattern (palindrome, good) or pattern pair (match, nested, good)
		if Data['LevenshteinAbsolute'] < 3:
			if KeyFirst == KeySecond: Data['Type'] = 'palindrome'
			else: 
				if Data['FirstLength'] == Data['SecondLength']: Data['Type'] = 'match'
				else: Data['Type'] = 'nested'
			Data['Risk'] = SimilarityType[Data['LevenshteinAbsolute']]
		else:
			Data['Type'] = 'good'
			Data['Risk'] = 'low'
		# Append new comparison item to stats
		Result.append(Data)
	# Return stats
	return Result

# Check risk of match, nesting and palindrome
def DistancesCheck(Distances):
	for item in Distances:
		if item['Risk'] == 'full':
			if item['Type'] == 'palindrome': warnings.warn(f'Pattern sequence "{item["FirstPattern"]}" is a palindrome. FR orientation is disabled.', RuntimeWarning)
			elif item['Type'] == 'nested': warnings.warn(f'Pattern sequences "{item["FirstPattern"]}" and "{item["SecondPattern"]}" are nested. Notice that order of patterns matters in search.', RuntimeWarning)
			elif item['Type'] == 'match': raise ValueError(f'Patterns "{item["FirstPattern"]}" and "{item["SecondPattern"]}" have the same sequences. Please remove one of them.')
		if item['Risk'] in ['high', 'medium']: 
			if item['Type'] == 'palindrome': warnings.warn(f'Pattern sequence "{item["FirstPattern"]}" has palindrome risk ({item["Risk"]}). This may impact statistics.', RuntimeWarning)
			elif item['Type'] == 'nested': warnings.warn(f'Pattern sequences "{item["FirstPattern"]}" and "{item["SecondPattern"]}" have nesting risk ({item["Risk"]}). This may impact statistics.', RuntimeWarning)
			elif item['Type'] == 'match': warnings.warn(f'Pattern sequences "{item["FirstPattern"]}" and "{item["SecondPattern"]}" have matching risk ({item["Risk"]}). This may impact statistics.', RuntimeWarning)

## ------======| I/O CHECK FUNC |======------

# Read compressed and uncompressed files
def ReadAnyway(FileName):
	# Check gzip and bz2 magic numbers
	GzipCheck = lambda FileName: open(FileName, 'rb').read(2).hex() == '1f8b'
	Bzip2Check = lambda FileName: open(FileName, 'rb').read(3).hex() == '425a68'
	# Choose module to handle file
	CompressedOpenFunc = lambda Type: { 'none': open, 'gzip': gzip.open, 'bz2': bz2.open }[Type]
	CheckFlags = GzipCheck(FileName = FileName), Bzip2Check(FileName = FileName)
	Type = { (1, 0): 'gzip', (0, 1): 'bz2', (0, 0): 'none' }[CheckFlags]
	# Return file stream object
	return CompressedOpenFunc(Type)(FileName, 'rt')

# Open FastQ files
def OpenFastQ(FileName):
	try:
		return SeqIO.parse(ReadAnyway(FileName), 'fastq')
	except Exception as Error:
		raise RuntimeError(f'Unable to open FastQ file: "{FileName}"\n{Error}')

# Check if output file is writable
def WritableCheck(FileName):
	with open(FileName, 'a') as File:
		if not File.writable(): raise OSError(f'File is not writable: {FileName}')

# Check K-mer size
def KmerSizeCheck(KmerMaxSize, Patterns):
	# Check positivity
	if KmerMaxSize < 0: raise ValueError(f'K-mer size must be non-negative: {KmerMaxSize}')
	# Warn if K-mer is longer than necessary
	LongestPattern = max([len(Value) for Value in Patterns.values()])
	if KmerMaxSize > LongestPattern: warnings.warn(f'K-mer size is bigger than the longest pattern sequence ({LongestPattern} bp). Notice that K-mer size bigger than necessary may impact statistics.', RuntimeWarning)

# Check unrecognized sequence alias
def UnrecognizedCheck(UnrecognizedSeq, Patterns):
	# Check if it contains only small Latin symbols
	if re.fullmatch('^[a-z]{2,16}$', UnrecognizedSeq) is None: raise ValueError(f'Unrecognized sequence name must contain 2-16 small Latin symbols: {UnrecognizedSeq}')
	# Check if there is no collision between aliases among patterns
	if UnrecognizedSeq in Patterns.keys(): raise ValueError(f'Unrecognized sequence name and one of pattern names are equal ({UnrecognizedSeq}). Rename one of them.')

# Check max read number to be processed
def MaxReadsCheck(MaxReads):
	# Check non-negativity
	if MaxReads < 0: raise ValueError(f'Max reads number must be non-negative: {MaxReads}')
	# Warn if max read number is greater than recommended
	if (MaxReads == 0) or (MaxReads > GLOBAL_MAX): warnings.warn(f'Max reads number is bigger than recommended. It may cause problems with device memory if FastQ files are huge.', RuntimeWarning)

# Check rate floor. It must be from 0 to 1
def RateFloorCheck(RateFloor):
	if (RateFloor >= 1) or (RateFloor < 0) or (RateFloor != RateFloor): raise ValueError(f'Rate floor must be float from 0 to 1: {RateFloor}')

# Check threads number
def ThreadsNumCheck(ThreadsNum):
	# Check positivity
	if ThreadsNum <= 0: raise ValueError(f'Threads number must be positive: {ThreadsNum}')
	# Warn if there is more threads than cores
	if GLOBAL_THREADS * 2 > ThreadsNum > GLOBAL_THREADS: warnings.warn(f'Threads number ({ThreadsNum}) is bigger than CPU count on your device ({GLOBAL_THREADS}). It may impact performance.', RuntimeWarning)
	# If threads number is much more than cores number, go on strike
	if ThreadsNum > GLOBAL_THREADS * 2: raise ValueError(f'Threads number is too big: {ThreadsNum}')

## ------======| PARSER |======------

def CreateParser():
	Default_parser = argparse.ArgumentParser(
			formatter_class = argparse.RawDescriptionHelpFormatter,
			description = f'FastContext: FastQ Context Analyzer, version {__version__}',
			epilog = f'Author: {__author__}'
			)
	Default_parser.add_argument ('-v', '--version', action = 'version', version = __version__)
	Default_parser.add_argument ('-1', '--r1', required = True, type = str, dest = 'InputFile_R1', help = f'FastQ input R1 file (may be gzipped or bzipped)')
	Default_parser.add_argument ('-2', '--r2', type = str, default = 'none', dest = 'InputFile_R2', help = f'FastQ input R2 file (may be gzipped or bzipped)')
	Default_parser.add_argument ('-p', '--patterns', required = True, type = str, dest = 'Patterns', help = f'Patterns to look for, plain JSON format: \'{{"first": "GATC", "second": "CTCAGCGCTGAG"}}\'. Names must contain 2-16 small Latin and numeric symbols (a-z, 0-9), sequences must contain more than one symbols ATGC. Order of patterns is order of search.')
	Default_parser.add_argument ('-t', '--tsv', required = True, type = str, dest = 'OutputTSV', help = f'Output TSV file. Contains only general statistics (read structure counts and percentages).')
	Default_parser.add_argument ('-j', '--json', type = str, default = 'none', dest = 'OutputJSON', help = f'Output JSON.GZ file (gzipped JSON). Contains extended statistics on every read')
	Default_parser.add_argument ('-k', '--kmer-size', default = GLOBAL_KMER_SIZE, type = int, dest = 'KmerMaxSize', help = f'Max size of unrecognized sequence to be written as K-mer of certain length. Default = {GLOBAL_KMER_SIZE}')
	Default_parser.add_argument ('-u', '--unrecognized', default = GLOBAL_UNRECOGNIZED, type = str, dest = 'UnrecognizedSeq', help = f'Long unrecognized sequences replacement. Default = "{GLOBAL_UNRECOGNIZED}"')
	Default_parser.add_argument ('-m', '--max-reads', default = GLOBAL_MAX, type = int, dest = 'MaxReads', help = f'Max reads number to analyze (0 - no limit). Default = {GLOBAL_MAX}')
	Default_parser.add_argument ('-f', '--rate-floor', default = GLOBAL_RATE_FLOOR, type = float, dest = 'RateFloor', help = f'Min rate to write read structure into stats TSV table. Default = {GLOBAL_RATE_FLOOR}')
	Default_parser.add_argument ('-@', '--threads', default = GLOBAL_THREADS, type = int, dest = 'ThreadsNum', help = f'Threads number. Default = {GLOBAL_THREADS}')
	Default_parser.add_argument ('-d', '--dont-check-read-names', action = 'store_true', dest = 'DontCheckReadNames', help = f'Don\'t check read names. Use this if you have unusual (non-Illumina) paired read names. Makes sense only in paired reads mode. Disabled by default.')
	Default_parser.add_argument ('-l', '--levenshtein', action = 'store_true', dest = 'CalculateLevenshtein', help = f'Calculate patterns Levenshtein distances for each position in read. Notice that it highly increases the time of processing.')
	return Default_parser

## ------======| ANALYSIS |======------

def AnalyzeSequence(Read, ToolKit):
	# Check every read in pair, if exists
	for Type in (('Read',) if ('Read' in Read) else ('R1', 'R2')):
		# Create Levenshtein data array
		if ToolKit['levenshtein']: Read[Type]['Levenshtein'] = []
		# Create two sequences: one for pattern, another for Levenshtein analysis
		Sequence = str(Read[Type]['Sequence']).upper()
		ProbabilitySequence = str(Read[Type]['Sequence']).upper()
		# For every pattern in toolkit
		for Pattern, PatternData in ToolKit['patterns'].items():
			# If pattern is a palindrome, disable strand
			Strands = ('F',) if (PatternData['F'] == PatternData['R']) else ('F', 'R')
			for Strand in Strands:
				# Create replacement string (JSON with all pattern data)
				Replacement = json.dumps({ 'type': 'pattern', 'name': str(Pattern), 'strand': None if PatternData['F'] == PatternData['R'] else str(Strand) }, separators = (',', ':')) 
				# Replace all matches
				Sequence = Sequence.replace(PatternData[Strand], Replacement)
				# Calculate Levenshtein distance for every position in read and every pattern
				if ToolKit['levenshtein']:
					LevenshteinVector = list()
					for FragmentStart in range(len(ProbabilitySequence)):
						Frag = ProbabilitySequence[FragmentStart:FragmentStart + len(PatternData[Strand])]
						LevenshteinVector.append(distance(Frag, PatternData[Strand]) if (len(Frag) == len(PatternData[Strand])) else float('NaN'))
					# Write Levenshtein distances data to stats
					Read[Type]['Levenshtein'].append({ 'name': str(Pattern), 'strand': None if PatternData['F'] == PatternData['R'] else str(Strand), 'length': int(PatternData['length']), 'values': LevenshteinVector.copy()})
		# Replace K-mers
		for Len in range(1, ToolKit['kmer_size']):
			# Create replacement string (JSON with all K-mer data)
			Replacement = json.dumps({ 'type': 'kmer', 'length': int(Len) }, separators = (',', ':')) 
			# Kostyl: add } and {, then replace K-mers, finally remove } and {
			Sequence = f'}}{Sequence}{{'
			Sequence = re.sub(f'\\}}[ATGCN]{{{str(Len)}}}\\{{', f'}}{Replacement}{{', Sequence)
			Sequence = Sequence[1:-1]
		# Replace unrecognized sequences
		Sequence = re.sub(f'[ATGCN]+', json.dumps({ 'type': 'unrecognized'}, separators = (',', ':')), Sequence)
		# Create a beautiful JSON from this Frankenstein's monster
		# Don't give a fuck, actually, it works and fuck it
		Sequence = Sequence.replace('}{', '},{')
		Sequence = f'[{Sequence}]'
		# Load JSON string
		Read[Type]['Pattern'] = json.loads(Sequence)
		# Create human-readable and machine-countable text pattern from dictionary
		# It will be used while calculating statistics
		TextPattern = list()
		for Sector in Read[Type]['Pattern']:
			if Sector['type'] == 'unrecognized': TextSector = str(GLOBAL_UNRECOGNIZED)
			elif Sector['type'] == 'pattern':
				if Sector['strand'] is None: TextSector = str(Sector['name'])
				else: TextSector = f'{Sector["name"]}:{Sector["strand"]}'
			elif Sector['type'] == 'kmer': TextSector = f'kmer:{Sector["length"]}bp'
			TextPattern.append(str(TextSector))
		TextPattern = '--'.join([f'{{{pat}}}' for pat in TextPattern])
		# Write text pattern to stats
		Read[Type]['TextPattern'] = str(TextPattern)
	# Return annotated read or read pair
	return Read

## ------======| CALCULATE STATISTICS |======------

# Count Summary
def SummaryCount(Statistics, Tag):
	Counts = Counter([item[Tag]['TextPattern'] for item in Statistics['RawDataset']])
	return { Key: { 'Count': Value, 'Rate': Value / Statistics['ReadsAnalyzed'] } for Key, Value in Counts.items() }

# Create TSV table from Summary
def Summary2Table(Summary, RateFloor):
	Table = pandas.DataFrame(Summary).transpose().reset_index()
	# Filter by rate floor and sort by rate
	Table = Table[Table['Rate'] >= RateFloor].sort_values(by=['Rate'], ascending=False)
	Table['Count'] = Table['Count'].apply(int)
	# Calculate percentage
	Table['Rate'] = Table['Rate'].apply(lambda x: x * 100)
	Table = Table[['Count', 'Rate', 'index']].rename(columns = {'index': 'Pattern'})
	# Return TSV as string
	return Table.to_csv(sep='\t', index=False)

## ------======| MAIN |======------

def Main(Namespace):
	
	# Open FastQ Files
	print(f'# Open FastQ files ...', end='\n')
	InputR1 = OpenFastQ(FileName = Namespace.InputFile_R1)
	# If single-end, set InputR2 to None
	InputR2 = None if (Namespace.InputFile_R2 == 'none') else OpenFastQ(FileName = Namespace.InputFile_R2)
	
	# Check Patterns
	print(f'# Check patterns ...', end='\n')
	PatternsCheck(Patterns = Namespace.Patterns)
	# OK, they looks good, now load
	Patterns = json.loads(Namespace.Patterns)
	# More deep analysis of patterns
	LevenshteinDistances = PatternsLevenshtein(Patterns = Patterns)
	DistancesCheck(Distances = LevenshteinDistances)
	
	# Check Output Files Writeability: TSV and Raw JSON.GZ
	print(f'# Check output files writeability ...', end='\n')
	WritableCheck(FileName = Namespace.OutputTSV)
	if Namespace.OutputJSON != 'none': WritableCheck(FileName = Namespace.OutputJSON)
	
	# Other Input Variables Check
	print(f'# Check input variables ...', end='\n')
	KmerSizeCheck(KmerMaxSize = Namespace.KmerMaxSize, Patterns = Patterns)
	UnrecognizedCheck(UnrecognizedSeq = Namespace.UnrecognizedSeq, Patterns = Patterns)
	MaxReadsCheck(MaxReads = Namespace.MaxReads)
	RateFloorCheck(RateFloor = Namespace.RateFloor)
	ThreadsNumCheck(ThreadsNum = Namespace.ThreadsNum)
	
	# Create Statistics Object
	print(f'# Create statistics object ...', end='\n')
	Statistics = {}
	Statistics['FastQ'] = os.path.realpath(Namespace.InputFile_R1) if InputR2 is None else { 'R1': os.path.realpath(Namespace.InputFile_R1), 'R2': os.path.realpath(Namespace.InputFile_R2) }
	Statistics['ReadType'] = 'Single-end' if InputR2 is None else 'Paired-end'
	Statistics['ReadsAnalyzed'] = 0
	Statistics['MaxReads'] = int(Namespace.MaxReads)
	Statistics['Performance'] = { 'Threads': int(Namespace.ThreadsNum), 'Started': datetime.datetime.now().isoformat() }
	Statistics['Patterns'] = { str(Key): { 'F': str(Value), 'R': ReverseComplement(Value), 'Length': len(Value) } for Key, Value in Patterns.items() }
	Statistics['KmerMaxSize'] = int(Namespace.KmerMaxSize)
	Statistics['RateFloor'] = float(Namespace.RateFloor)
	Statistics['UnrecognizedSequence'] = str(Namespace.UnrecognizedSeq)
	Statistics['PatternsLevenshtein'] = LevenshteinDistances
	Statistics['Summary'] = { 'Read': {} } if InputR2 is None else { 'R1': {}, 'R2': {} }
	Statistics['RawDataset'] = None
	
	# Create ToolKit Object
	# Only variables that are necessary for thread function
	ToolKit = {
		'kmer_size': int(Namespace.KmerMaxSize),
		'unrecognized': str(Namespace.UnrecognizedSeq),
		'patterns': { str(pattern): { 'F': str(seq), 'R': ReverseComplement(seq), 'length': len(seq) } for pattern, seq in Patterns.items() },
		'levenshtein': bool(Namespace.CalculateLevenshtein)
		}
	
	# Create Dataset
	Dataset = list()
	
	# Read FastQ Files. tqdm progressbar is enabled
	print(f'# Start FastQ parsing ...', end='\n')
	with tqdm.tqdm(desc='Reads parsed', unit='') as pbar:
		while 1:
			try:
				Read1 = next(InputR1)
				# If paired-end, read InputR2
				if InputR2 is not None: 
					Read2 = next(InputR2)
					# Check Read Names
					if (Read1.name != Read2.name) and not (Namespace.DontCheckReadNames): raise RuntimeError(f'Read names in pair are not the same: "{Read1.name}", "{Read2.name}". If you\'re sure your FastQ files are valid, run the script with -d option.')
					# Add read to dataset. Notice that Phred quals are added too
					Dataset.append({ 'Name': str(Read1.name), 'R1': { 'Sequence': Read1.seq.__str__(), 'PhredQual': Read1.letter_annotations['phred_quality'].copy() }, 'R2': { 'Sequence': Read2.seq.__str__(), 'PhredQual': Read2.letter_annotations['phred_quality'].copy() } })
				# If single-end, just add the read
				else: Dataset.append({ 'Name': str(Read1.name), 'Read': { 'Sequence': Read1.seq.__str__(), 'PhredQual': Read1.letter_annotations['phred_quality'].copy() } } )
				# Increment and check read counter
				Statistics['ReadsAnalyzed'] += 1
				pbar.update(1)
				if Statistics['ReadsAnalyzed'] == Namespace.MaxReads: break
			# Break if file has ended ...
			except StopIteration:
				break
			# ... or we have some troubles with that
			except ValueError as Error:
				raise ValueError(f'File is not a valid FastQ file ({Error})')
	
	# Empty files are also not allowed
	if not Dataset: raise ValueError(f'FastQ file is empty.')
	
	# Read analysis
	print(f'# Start read analysis on {Namespace.ThreadsNum} threads ...', end='\n')
	with Threading(Namespace.ThreadsNum) as Pool:
		Statistics['RawDataset'] = list(tqdm.tqdm(Pool.imap(functools.partial(AnalyzeSequence, ToolKit=ToolKit), Dataset), total=len(Dataset), desc=f'Reads analyzed'))
	
	# Calculate statistics
	Tags = ('Read',) if InputR2 is None else ('R1', 'R2')
	print(f'# Calculate statistics ...', end='\n')
	for Tag in Tags: Statistics['Summary'][Tag] = SummaryCount(Statistics, Tag)
	
	# Write stats to table
	print(f'# Create statistics TSV table ...', end='\n')
	Tables = {}
	for Tag in Tags: Tables[Tag] = Summary2Table(Statistics['Summary'][Tag], Namespace.RateFloor)
	TextTable = ''
	for Key, Value in Tables.items(): TextTable += f'# {Key}\n{Value}\n'
	with open(Namespace.OutputTSV, 'wt') as TSV: TSV.write(TextTable)
	
	# Timestamp
	Statistics['Performance']['Finished'] = datetime.datetime.now().isoformat()
	
	# Write JSON.GZ raw data if necessary
	if Namespace.OutputJSON != 'none':
		print(f'# Save extended statistics (JSON.GZ) data ...', end='\n')
		json.dump(Statistics, gzip.open(Namespace.OutputJSON, 'wt'), separators=(',', ':'))
	
	# Bye-bye, baby
	print(f'# FastQ analysis ended successfully.', end='\n')

if __name__ == '__main__':
	
	Parser = CreateParser()
	Namespace = Parser.parse_args(sys.argv[1:])
	
	Main(Namespace)
