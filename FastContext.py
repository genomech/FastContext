#!/bin/python3

__version__ = "0.9b"
__author__ = "Emil Viesná"

from Bio import SeqIO
from Bio.Seq import Seq
from Levenshtein import distance
from collections import Counter
import contextlib
import multiprocessing
import argparse
import bz2
import functools
import gzip
import itertools
import json
import logging
import os
import re
import pandas
import string
import sys
import time
import warnings
import tqdm

def ReverseComplement(Sequence): return Seq(Sequence).reverse_complement().__str__()

def OneLinerWarning(message, category, filename, lineno, file=None, line=None): return '%s: %s\n' % (category.__name__, message)
warnings.formatwarning = OneLinerWarning

sys.tracebacklimit = 0

GLOBAL_KMER_SIZE = 0
GLOBAL_UNRECOGNIZED = 'genome'
GLOBAL_MAX = 1000000
GLOBAL_THREADS = multiprocessing.cpu_count()
GLOBAL_RATE_FLOOR = 0.001

## ------======| I/O |======------

def CurrentDir(): return os.path.dirname(os.path.abspath(__file__))

def ReadAnyway(FileName):
	GzipCheck = lambda FileName: open(FileName, 'rb').read(2).hex() == '1f8b'
	Bzip2Check = lambda FileName: open(FileName, 'rb').read(3).hex() == '425a68'
	CompressedOpenFunc = lambda Type: { 'none': open, 'gzip': gzip.open, 'bz2': bz2.open }[Type]
	CheckFlags = GzipCheck(FileName = FileName), Bzip2Check(FileName = FileName)
	Type = { (1, 0): 'gzip', (0, 1): 'bz2', (0, 0): 'none' }[CheckFlags]
	return CompressedOpenFunc(Type)(FileName, 'rt')

## ------======| THREADING |======------

@contextlib.contextmanager
def Threading(Threads = multiprocessing.cpu_count()):
	StartTime = time.time()
	pool = multiprocessing.Pool(Threads)
	yield pool
	pool.close()
	pool.join()
	del pool

## ------======| CHECK FUNC |======------

def ParseFastQ(FileName):
	try:
		return SeqIO.parse(ReadAnyway(FileName), 'fastq')
	except Exception as Error:
		raise RuntimeError(f'Unable to open FastQ file: "{FileName}"\n{Error}')

def PatternsCheck(Patterns):
	LengthSorted = lambda x: all(len(x[i]) >= len(x[i + 1]) for i in range(len(x) - 1))
	try:
		PatternsDict = json.loads(Patterns)
	except json.decoder.JSONDecodeError as Error:
		raise ValueError(f'Patterns string is not a valid JSON string: {Patterns}\n{Error}')
	if type(PatternsDict) != dict: raise ValueError(f'Patterns string is not a JavaScript Object (dictionary): {Patterns}')
	if len(PatternsDict) < 1: raise ValueError(f'Patterns string does not contain any patterns: {Patterns}')
	InvertedPatterns = {}
	for Key, Value in PatternsDict.items():
		if (type(Key) != str): raise ValueError(f'Invalid type of pattern name ({type(Key)}), must be string: {Key}')
		if (type(Value) != str): raise ValueError(f'Invalid type of pattern sequence ({type(Value)}), must be string: {Value}')
		if re.fullmatch('^[a-z]{2,16}$', Key) is None: raise ValueError(f'Pattern name must contain 2-16 small Latin symbols: {Key}')
		if re.fullmatch('^[ATGC]+$', Value) is None: raise ValueError(f'Pattern sequence must contain more than one symbols ATGC: {Value}')
		if len(Value) > 16: warnings.warn(f'Pattern "{Key}" is very long ({len(Value)} bp). This may impact search. You may want to use shorter fragments of it.', RuntimeWarning)
	if not LengthSorted(list(PatternsDict.values())): warnings.warn(f'Patterns are not sorted from longest to shortest. Notice that order of patterns matters in search.', RuntimeWarning)

def PatternsLevenshtein(Patterns):
	SimilarityType = { 0: 'full', 1: 'high', 2: 'medium' }
	Result = []
	for KeyX, KeyY in itertools.combinations_with_replacement(Patterns.keys(), 2):
		BlockFR = {
			'FX': str(Patterns[KeyX]),
			'FY': str(Patterns[KeyY]),
			'RX': ReverseComplement(Patterns[KeyX]),
			'RY': ReverseComplement(Patterns[KeyY])
		}
		LevenshteinDistAbs = float('+inf')
		Longer, Shorter = ('X', 'Y') if len(BlockFR['FX']) >= len(BlockFR['FY']) else ('Y', 'X')
		for Second in ['F', 'R']:
			for FragmentStart in range(len(BlockFR['F' + Longer]) - len(BlockFR[Second + Shorter]) + 1):
				Dist = distance(BlockFR['F' + Longer][FragmentStart:FragmentStart + len(BlockFR[Second + Shorter])], BlockFR[Second + Shorter])
				if (Dist < LevenshteinDistAbs) and not ((KeyX == KeyY) and (Second == 'F')): LevenshteinDistAbs = int(Dist)
		Data = {
			'Analysis': 'full' if KeyX != KeyY else 'rc only',
			'Pattern #1': str(KeyX),
			'Pattern #2': str(KeyY),
			'Length #1': len(Patterns[KeyX]),
			'Length #2': len(Patterns[KeyY]),
			'Levenshtein Absolute': int(LevenshteinDistAbs),
			'Levenshtein Similarity': float((len(BlockFR['F' + Shorter]) - LevenshteinDistAbs) / len(BlockFR['F' + Shorter]))
			}
		if Data['Levenshtein Absolute'] < 3:
			if KeyX == KeyY: Data['Type'] = 'palindrome'
			else: 
				if Data['Length #1'] == Data['Length #2']: Data['Type'] = 'match'
				else: Data['Type'] = 'nested'
			Data['Risk'] = SimilarityType[Data['Levenshtein Absolute']]
		else:
			Data['Type'] = 'good'
			Data['Risk'] = 'low'
		Result.append(Data)
	return Result

def DistancesCheck(Distances):
	for item in Distances:
		if item['Risk'] == 'full':
			if item['Type'] == 'palindrome': warnings.warn(f'Pattern "{item["Pattern #1"]}" is a palindrome. FR orientation is disabled.', RuntimeWarning)
			elif item['Type'] == 'nested': warnings.warn(f'Patterns "{item["Pattern #1"]}" and "{item["Pattern #2"]}" are nested. Notice that order of patterns matters in search.', RuntimeWarning)
			elif item['Type'] == 'match': raise ValueError(f'Patterns "{item["Pattern #1"]}" and "{item["Pattern #2"]}" have the same sequences. Please remove one of them.')
		if item['Risk'] in ['high', 'medium']: 
			if item['Type'] == 'palindrome': warnings.warn(f'Pattern "{item["Pattern #1"]}" has palindrome risk ({item["Risk"]}). This may impact statistics.', RuntimeWarning)
			elif item['Type'] == 'nested': warnings.warn(f'Patterns "{item["Pattern #1"]}" and "{item["Pattern #2"]}" have nesting risk ({item["Risk"]}). This may impact statistics.', RuntimeWarning)
			elif item['Type'] == 'match': warnings.warn(f'Patterns "{item["Pattern #1"]}" and "{item["Pattern #2"]}" have matching risk ({item["Risk"]}). This may impact statistics.', RuntimeWarning)


def WritableCheck(FileName):
	with open(FileName, 'a') as File:
		if not File.writable(): raise OSError(f'File is not writable: {FileName}')

def KmerSizeCheck(KmerMaxSize, Patterns):
	if KmerMaxSize < 0: raise ValueError(f'K-mer size must be non-negative: {KmerMaxSize}')
	LongestPattern = max([len(Value) for Value in Patterns.values()])
	if KmerMaxSize > LongestPattern: warnings.warn(f'K-mer size is bigger than the longest pattern ({LongestPattern} bp). Notice that big K-mer size may impact statistics.', RuntimeWarning)

def UnrecognizedCheck(UnrecognizedSeq, Patterns):
	if re.fullmatch('^[a-z]{2,16}$', UnrecognizedSeq) is None: raise ValueError(f'Unrecognized sequence name must contain 2-16 small Latin symbols: {UnrecognizedSeq}')
	if UnrecognizedSeq in Patterns.keys(): raise ValueError(f'Unrecognized sequence name and one of pattern names are equal.')

def MaxReadsCheck(MaxReads):
	if MaxReads < 0: raise ValueError(f'Max reads number must be positive: {MaxReads}')
	if (MaxReads == 0) or (MaxReads > GLOBAL_MAX): warnings.warn(f'Max reads number is bigger than recommended. It may cause problems with device memory if FastQ files are huge.', RuntimeWarning)

def RateFloorCheck(RateFloor):
	if (RateFloor >= 1) or (RateFloor < 0) or (RateFloor != RateFloor): raise ValueError(f'Rate floor must be from 0 to 1: {RateFloor}')

def ThreadsNumCheck(ThreadsNum):
	if ThreadsNum <= 0: raise ValueError(f'Threads number must be positive: {ThreadsNum}')
	if GLOBAL_THREADS * 4 > ThreadsNum > GLOBAL_THREADS: warnings.warn(f'Threads number ({ThreadsNum}) is bigger than CPU count on your device ({GLOBAL_THREADS}). It may impact performance.', RuntimeWarning)
	if ThreadsNum > GLOBAL_THREADS * 2: raise ValueError(f'Threads number is too big: {ThreadsNum}')

## ------======| PARSER |======------

def CreateParser():
	Default_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=f"FastContext: FastQ Context Analyzer, version {__version__}", epilog=f"Author: {__author__}")
	Default_parser.add_argument ('-v', '--version', action='version', version=__version__)
	Default_parser.add_argument ('-1', '--r1', required=True, type=str, dest='InputFile_R1', help=f'FastQ input R1 file (may be gzipped or bzipped)')
	Default_parser.add_argument ('-2', '--r2', type=str, default='none', dest='InputFile_R2', help=f'FastQ input R2 file (may be gzipped or bzipped)')
	Default_parser.add_argument ('-p', '--patterns', required=True, type=str, dest='Patterns', help=f'Patterns to look for, plain JSON format: \'{{"first": "GATC", "second": "CTCAGCGCTGAG"}}\'. Names must contain 2-16 small Latin symbols (a-z), sequences must contain more than one symbols ATGC. Order of patterns is order of search.')
	Default_parser.add_argument ('-t', '--tsv', required=True, type=str, dest='OutputTSV', help=f'Output TSV file. Contains only general statistics.')
	Default_parser.add_argument ('-j', '--json', type=str, default='none', dest='OutputJSON', help=f'Output JSON.GZ file (gzipped JSON, raw data). Useful if you need to see particular read structure')
	Default_parser.add_argument ('-k', '--kmer-size', default=GLOBAL_KMER_SIZE, type=int, dest='KmerMaxSize', help=f'Max unrecognized K-mer size. Default = {GLOBAL_KMER_SIZE}')
	Default_parser.add_argument ('-u', '--unrecognized', default=GLOBAL_UNRECOGNIZED, type=str, dest='UnrecognizedSeq', help=f'Long unrecognized sequences replacement. Default = "{GLOBAL_UNRECOGNIZED}"')
	Default_parser.add_argument ('-m', '--max-reads', default=GLOBAL_MAX, type=int, dest='MaxReads', help=f'Max reads number to analyze (0 - no limit). Default = {GLOBAL_MAX}')
	Default_parser.add_argument ('-f', '--rate-floor', default=GLOBAL_RATE_FLOOR, type=float, dest='RateFloor', help=f'Min rate to write pattern into stats TSV table. Default = {GLOBAL_RATE_FLOOR}')
	Default_parser.add_argument ('-@', '--threads', default=GLOBAL_THREADS, type=int, dest='ThreadsNum', help=f'Threads number. Default = {GLOBAL_THREADS}')
	Default_parser.add_argument ('-d', '--dont-check-read-names', action='store_true', dest='DontCheckReadNames', help=f'Don\'t check read names. Use this if you have unusual paired read names. Makes sense only in paired reads mode. Disabled by default.')
	return Default_parser

## ------======| ANALYSIS |======------

def AnalyzeSequence(Read, ToolKit):
	for Type in (('Seq',) if ('Seq' in Read) else ('R1', 'R2')):
		Read['Probability ' + Type] = {}
		Sequence = str(Read[Type]).upper()
		ProbabilitySequence = str(Read[Type]).upper()
		for Pattern, PatternData in ToolKit['patterns'].items():
			if PatternData['F'] == PatternData['R']:
				Sequence = Sequence.replace(PatternData['F'], f'[{Pattern}:P]')
				Read['Probability ' + Type][f'{Pattern}:P'] = []
				for FragmentStart in range(len(ProbabilitySequence)):
					Frag = ProbabilitySequence[FragmentStart:FragmentStart + len(PatternData['F'])]
					Read['Probability ' + Type][f'{Pattern}:P'].append(distance(Frag, PatternData['F']) if (len(Frag) == len(PatternData['F'])) else float('NaN'))
			else:
				for Strand in ('F', 'R'):
					Sequence = Sequence.replace(PatternData[Strand], f"[{Pattern}:{Strand}]")
					Read['Probability ' + Type][f'{Pattern}:{Strand}'] = []
					for FragmentStart in range(len(ProbabilitySequence)):
						Frag = ProbabilitySequence[FragmentStart:FragmentStart + len(PatternData[Strand])]
						Read['Probability ' + Type][f'{Pattern}:{Strand}'].append(distance(Frag, PatternData[Strand]) if (len(Frag) == len(PatternData[Strand])) else float('NaN'))
		for Len in range(1, ToolKit['kmer_size']): Sequence = re.sub(f'(^|\\])[ATGCN]{{{str(Len)}}}($|\\[)', f'\\1[{str(Len)}bp]\\2', Sequence)
		Sequence = re.sub(f'[ATGCN]+', f'[{ToolKit["unrecognized"]}]', Sequence)
		Sequence = Sequence.replace('][', ']-[')
		Read['Pattern ' + Type] = str(Sequence)
	return Read

## ------======| STATISTICS |======------

def SummaryCount(Statistics, PatternTag):
	Counts = Counter([item[PatternTag] for item in Statistics['Raw Dataset']])
	return { Key: { 'Count': Value, 'Rate': Value / Statistics['Reads Analyzed'] } for Key, Value in Counts.items() }

def Summary2Table(Summary, RateFloor):
	Table = pandas.DataFrame(Summary).transpose().reset_index()
	Table = Table[Table['Rate'] >= RateFloor].sort_values(by=['Rate'], ascending=False)
	Table['Count'] = Table['Count'].apply(int)
	Table['Rate'] = Table['Rate'].apply(lambda x: x * 100)
	Table = Table[['Count', 'Rate', 'index']].rename(columns={'index': 'Pattern'})
	return Table.to_csv(sep='\t', index=False)

## ------======| MAIN |======------

def Main(Namespace):
	
	# Open FastQ Files
	print(f'# Open FastQ files...', end='\n')
	InputR1 = ParseFastQ(Namespace.InputFile_R1)
	InputR2 = None if (Namespace.InputFile_R2 == 'none') else ParseFastQ(Namespace.InputFile_R2)
	
	# Check & Load Patterns
	print(f'# Check patterns...', end='\n')
	PatternsCheck(Namespace.Patterns)
	Patterns = json.loads(Namespace.Patterns)
	LevenshteinDistances = PatternsLevenshtein(Patterns)
	DistancesCheck(LevenshteinDistances)
	
	# Check Output Files Writeability
	print(f'# Check output files writeability...', end='\n')
	WritableCheck(Namespace.OutputTSV)
	if Namespace.OutputJSON != 'none': WritableCheck(Namespace.OutputJSON)
	
	# Other Input Variables Check
	print(f'# Check input variables...', end='\n')
	KmerSizeCheck(Namespace.KmerMaxSize, Patterns)
	UnrecognizedCheck(Namespace.UnrecognizedSeq, Patterns)
	MaxReadsCheck(Namespace.MaxReads)
	RateFloorCheck(Namespace.RateFloor)
	ThreadsNumCheck(Namespace.ThreadsNum)
	
	# Create Statistics Object
	print(f'# Create statistics object...', end='\n')
	Statistics = {}
	Statistics['FastQ'] = os.path.realpath(Namespace.InputFile_R1) if InputR2 is None else { 'R1': os.path.realpath(Namespace.InputFile_R1), 'R2': os.path.realpath(Namespace.InputFile_R2) }
	Statistics['Read Type'] = 'Single-end' if InputR2 is None else 'Paired-end'
	Statistics['Reads Analyzed'] = 0
	Statistics['Patterns'] = { str(Key): { 'F': str(Value), 'R': ReverseComplement(Value), 'Length': len(Value) } for Key, Value in Patterns.items() }
	Statistics['K-mer Max Size'] = int(Namespace.KmerMaxSize)
	Statistics['Unrecognized Sequence'] = str(Namespace.UnrecognizedSeq)
	Statistics['Patterns Levenshtein'] = LevenshteinDistances
	Statistics['Summary'] = {} if InputR2 is None else { 'R1': {}, 'R2': {} }
	Statistics['Raw Dataset'] = None
	
	# Create ToolKit Object for Threading
	ToolKit = {
		'kmer_size': int(Namespace.KmerMaxSize),
		'unrecognized': str(Namespace.UnrecognizedSeq),
		'patterns': { str(pattern): { 'F': str(seq), 'R': ReverseComplement(seq) } for pattern, seq in Patterns.items() }
		}
	
	# Reads Dataset
	Dataset = []
	
	# Read FastQ Files
	print(f'# Start FastQ parsing...', end='\n')
	with tqdm.tqdm(desc='Reads parsed', unit='') as pbar:
		while 1:
			try:
				Read1 = next(InputR1)
				if InputR2 is not None: 
					Read2 = next(InputR2)
					# Check Read Names
					if (Read1.name != Read2.name) and not (Namespace.DontCheckReadNames): raise RuntimeError(f'Read names in pair are not the same: "{Read1.name}", "{Read2.name}". If you\'re sure your files are valid, run the script with -d option.')
					Dataset.append({ 'Name': str(Read1.name), 'R1': Read1.seq.__str__(), 'R2': Read2.seq.__str__() })
				else:
					Dataset.append({ 'Name': str(Read1.name), 'Seq': Read1.seq.__str__() })
				Statistics['Reads Analyzed'] += 1
				pbar.update(1)
				if Statistics['Reads Analyzed'] == Namespace.MaxReads: break
			except StopIteration:
				break
			except ValueError as Error:
				raise ValueError(f'File is not a valid FastQ file ({Error})')
	
	if not Dataset: raise ValueError(f'FastQ file is empty.')
	
	print(f'# Start read analysis on {Namespace.ThreadsNum} threads...', end='\n')
	with Threading(Namespace.ThreadsNum) as Pool:
		Statistics['Raw Dataset'] = list(tqdm.tqdm(Pool.imap(functools.partial(AnalyzeSequence, ToolKit=ToolKit), Dataset), total=len(Dataset), desc=f'Reads analyzed'))
	
	if InputR2 is None:
		print(f'# Calculate statistics... ', end='\n')
		Statistics['Summary'] = SummaryCount(Statistics, 'Pattern Seq')
		print(f'# Create statistics table... ', end='\n')
		Table = Summary2Table(Statistics['Summary'], Namespace.RateFloor)
		with open(Namespace.OutputTSV, 'wt') as TSV: TSV.write(Table)
	else:
		print(f'# Calculate statistics... ', end='\n')
		Statistics['Summary']['R1'] = SummaryCount(Statistics, 'Pattern R1')
		Statistics['Summary']['R2'] = SummaryCount(Statistics, 'Pattern R2')
		print(f'# Create statistics table... ', end='\n')
		TableR1 = Summary2Table(Statistics['Summary']['R1'], Namespace.RateFloor)
		TableR2 = Summary2Table(Statistics['Summary']['R2'], Namespace.RateFloor)
		with open(Namespace.OutputTSV, 'wt') as TSV: TSV.write('# R1\n' + TableR1 + '\n# R2\n' + TableR2)
	
	if Namespace.OutputJSON != 'none':
		print(f'# Save raw JSON data... ', end='\n')
		json.dump(Statistics, gzip.open(Namespace.OutputJSON, 'wt'), separators=(',', ':'))
		print(f'# FastQ analysis ended successfully.', end='\n')
	
if __name__ == '__main__':
	
	Parser = CreateParser()
	Namespace = Parser.parse_args(sys.argv[1:])
	
	Main(Namespace)
