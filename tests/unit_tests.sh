ScriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd; )

TestScript=$( realpath ${ScriptDir}/../FastContext.py )

echo '--- FILE NOT EXIST ---'
python3 ${TestScript} -1 ${ScriptDir}/file_not_exists.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- NOT A FASTQ ---'
python3 ${TestScript} -1 ${ScriptDir}/not_a_fastq.txt -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- EMPTY FILE ---'
python3 ${TestScript} -1 ${ScriptDir}/empty_file.txt -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- UNREADABLE FILE ---'
chmod 000 ${ScriptDir}/unreadable_file.fastq
python3 ${TestScript} -1 ${ScriptDir}/unreadable_file.fastq -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
chmod 755 ${ScriptDir}/unreadable_file.fastq
echo

echo '--- CORRUPTED NAMES ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_corruptednames_R1.fastq.gz -2 ${ScriptDir}/gzipped_corruptednames_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- UNUSUAL NAMES WITH -d ---'
python3 ${TestScript} -1 ${ScriptDir}/nocompression_unusualnames_R1.fq -2 ${ScriptDir}/nocompression_unusualnames_R2.fq -d -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- GZIPPED NORMAL SINGLE-END ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- BZIPPED NORMAL PAIRED-END ---'
python3 ${TestScript} -1 ${ScriptDir}/bzipped_normal_R1.fastq.bz2 -2 ${ScriptDir}/bzipped_normal_R2.fastq.bz2 -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- NOT JSON ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p 'not_json' -s /dev/null
echo

echo '--- EMPTY JSON ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{}' -s /dev/null
echo

echo '--- INVALID JSON FORMAT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{ "foo": [ "bar" ], "baz": null, "spam": "eggs" }' -s /dev/null
echo

echo '--- INVALID JSON TYPE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '[ "bar", "baz", "ATGC" ]' -s /dev/null
echo

echo '--- INVALID PATTERN NAME ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "Duck": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- INVALID PATTERN SEQUENCE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGGa"}' -s /dev/null
echo

echo '--- TOO LONG PATTERN SEQUENCE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACACCCGAGCCCACG", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- NOT SORTED PATTERN SEQUENCES ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"bar": "CCGAAAACACG", "foo": "CTGTCTCTTATACAC", "baz": "TCGTCGGG"}' -s /dev/null
echo

echo '--- PALINDROME ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTCAGCGCTGAG", "bar": "GATCGATCC", "baz": "GGTACC"}' -s /dev/null
echo

echo '--- FULL MATCH ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "oof": "GTGTATAAGAGACAG"}' -s /dev/null
echo

echo '--- MATCH DANGER ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "GTGTATATGAGACAG", "baz": "CTTTCTCTTATACAC"}' -s /dev/null
echo

echo '--- NESTING DANGER ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CTGTCTCT", "baz": "AACAC"}' -s /dev/null
echo

echo '--- UNWRITABLE OUTPUT ---'
chmod 444 ${ScriptDir}/unwriteable_output_file.htm
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s ${ScriptDir}/unwriteable_output_file.htm
chmod 755 ${ScriptDir}/unwriteable_output_file.htm
echo

echo '--- KMER SIZE NEGATIVE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -k '-100' -s /dev/null
echo
# LOL It takes -100 for an argument
# Not as expected, but ok

echo '--- KMER SIZE TOO BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -k 50 -s /dev/null
echo

echo '--- KMER SIZE FLOAT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -k 64.9 -s /dev/null
echo

echo '--- KMER SIZE STRING ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -k 'spam' -s /dev/null
echo

echo '--- UNRECOGNIZED ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -u 'baz' -s /dev/null
echo

echo '--- MAX READS NEGATIVE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -m '-100' -s /dev/null
echo
# It takes -100 for an argument

echo '--- MAX READS TOO BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -m 0 -s /dev/null
echo

echo '--- MAX READS FLOAT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -m 64.9 -s /dev/null
echo

echo '--- MAX READS STRING ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -m 'spam' -s /dev/null
echo

echo '--- RATE FLOOR NEGATIVE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f '-100.56' -s /dev/null
echo
# It takes -100 for an argument

echo '--- RATE FLOOR TOO BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f 2.67 -s /dev/null
echo

echo '--- RATE FLOOR E ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f 1.9e-14 -s /dev/null
echo

echo '--- RATE FLOOR INT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f 0 -s /dev/null
echo

echo '--- RATE FLOOR NAN ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f NaN -s /dev/null
echo

echo '--- RATE FLOOR INF ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f +inf -s /dev/null
echo

echo '--- RATE FLOOR STRING ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f 'spam' -s /dev/null
echo

echo '--- THREADS NUM NEGATIVE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ '-100' -s /dev/null
echo
# It takes -100 for an argument

echo '--- THREADS NUM A LITTLE BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ 13 -s /dev/null
echo

echo '--- THREADS NUM TOO BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ 800 -s /dev/null
echo

echo '--- THREADS NUM FLOAT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ 64.9e-10 -s /dev/null
echo

echo '--- THREADS NUM STRING ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ 'spam' -s /dev/null
echo

echo '--- EMPTY SUMMARY ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -s ${ScriptDir}/empty_summary_test.htm -f 0.999
echo

echo '--- NORMAL TEST SINGLE-END HTML EXAMPLE ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "s502": "CTCTCTAT", "s702": "CTAGTACG", "s701": "TCGCCTTA" }' -k 15 -s ${ScriptDir}/standard_test_stats_se.htm
echo

echo '--- NORMAL TEST PAIRED-END HTML EXAMPLE ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -2 ${ScriptDir}/standard_test_R2.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "s502": "CTCTCTAT", "s702": "CTAGTACG", "s701": "TCGCCTTA" }' -k 15 -s ${ScriptDir}/standard_test_stats_pe.htm
echo

echo '--- NORMAL TEST SINGLE-END JSON EXAMPLE ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "s502": "CTCTCTAT", "s702": "CTAGTACG", "s701": "TCGCCTTA" }' -k 15 -m 100 -s /dev/null -j ${ScriptDir}/standard_test_raw_se.json.gz
echo

echo '--- NORMAL TEST PAIRED-END JSON EXAMPLE ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -2 ${ScriptDir}/standard_test_R2.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "s502": "CTCTCTAT", "s702": "CTAGTACG", "s701": "TCGCCTTA" }' -k 15 -m 100 -s /dev/null -j ${ScriptDir}/standard_test_raw_pe.json.gz
echo

echo '--- NORMAL TEST SINGLE-END JSON EXAMPLE, LEVENSHTEIN ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "s502": "CTCTCTAT", "s702": "CTAGTACG", "s701": "TCGCCTTA" }' -k 15 -m 100 -s /dev/null -j ${ScriptDir}/standard_test_raw_levenshtein_se.json.gz -l
echo

echo '--- NORMAL TEST PAIRED-END JSON EXAMPLE, LEVENSHTEIN ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -2 ${ScriptDir}/standard_test_R2.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "s502": "CTCTCTAT", "s702": "CTAGTACG", "s701": "TCGCCTTA" }' -k 15 -m 100 -s /dev/null -j ${ScriptDir}/standard_test_raw_levenshtein_pe.json.gz -l
echo
