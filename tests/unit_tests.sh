ScriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd; )

TestScript=$( realpath ${ScriptDir}/../FastContext.py )

echo '--- FILE NOT EXIST ---'
python3 ${TestScript} -1 ${ScriptDir}/file_not_exists.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- NOT A FASTQ ---'
python3 ${TestScript} -1 ${ScriptDir}/not_a_fastq.txt -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- EMPTY FILE ---'
python3 ${TestScript} -1 ${ScriptDir}/empty_file.txt -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- UNREADABLE FILE ---'
chmod 000 ${ScriptDir}/unreadable_file.fastq
python3 ${TestScript} -1 ${ScriptDir}/unreadable_file.fastq -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
chmod 755 ${ScriptDir}/unreadable_file.fastq
echo

echo '--- CORRUPTED NAMES ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_corruptednames_R1.fastq.gz -2 ${ScriptDir}/gzipped_corruptednames_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- UNUSUAL NAMES WITH -d ---'
python3 ${TestScript} -1 ${ScriptDir}/nocompression_unusualnames_R1.fq -2 ${ScriptDir}/nocompression_unusualnames_R2.fq -d -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- GZIPPED NORMAL SINGLE-END ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- BZIPPED NORMAL PAIRED-END ---'
python3 ${TestScript} -1 ${ScriptDir}/bzipped_normal_R1.fastq.bz2 -2 ${ScriptDir}/bzipped_normal_R2.fastq.bz2 -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- NOT JSON ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p 'not_json' -t /dev/null
echo

echo '--- EMPTY JSON ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{}' -t /dev/null
echo

echo '--- INVALID JSON FORMAT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{ "foo": [ "bar" ], "baz": null, "spam": "eggs" }' -t /dev/null
echo

echo '--- INVALID JSON TYPE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '[ "bar", "baz", "ATGC" ]' -t /dev/null
echo

echo '--- INVALID PATTERN NAME ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "Duck": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- INVALID PATTERN SEQUENCE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGGa"}' -t /dev/null
echo

echo '--- TOO LONG PATTERN SEQUENCE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACACCCGAGCCCACG", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- NOT SORTED PATTERN SEQUENCES ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"bar": "CCGAAAACACG", "foo": "CTGTCTCTTATACAC", "baz": "TCGTCGGG"}' -t /dev/null
echo

echo '--- PALINDROME ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTCAGCGCTGAG", "bar": "GATCGATCC", "baz": "GGTACC"}' -t /dev/null
echo

echo '--- FULL MATCH ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "oof": "GTGTATAAGAGACAG"}' -t /dev/null
echo

echo '--- MATCH DANGER ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "GTGTATATGAGACAG", "baz": "CTTTCTCTTATACAC"}' -t /dev/null
echo

echo '--- NESTING DANGER ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -2 ${ScriptDir}/gzipped_normal_R2.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CTGTCTCT", "baz": "AACAC"}' -t /dev/null
echo

echo '--- UNWRITABLE OUTPUT ---'
chmod 444 ${ScriptDir}/unwriteable_output_file.tsv
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -t ${ScriptDir}/unwriteable_output_file.tsv
chmod 755 ${ScriptDir}/unwriteable_output_file.tsv
echo

echo '--- KMER SIZE NEGATIVE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -k '-100' -t /dev/null
echo
# LOL It takes -100 for an argument
# Not as expected, but ok

echo '--- KMER SIZE TOO BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -k 50 -t /dev/null
echo

echo '--- KMER SIZE FLOAT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -k 64.9 -t /dev/null
echo

echo '--- KMER SIZE STRING ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -k 'spam' -t /dev/null
echo

echo '--- UNRECOGNIZED ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -u 'baz' -t /dev/null
echo

echo '--- MAX READS NEGATIVE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -m '-100' -t /dev/null
echo
# It takes -100 for an argument

echo '--- MAX READS TOO BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -m 0 -t /dev/null
echo

echo '--- MAX READS FLOAT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -m 64.9 -t /dev/null
echo

echo '--- MAX READS STRING ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -m 'spam' -t /dev/null
echo

echo '--- RATE FLOOR NEGATIVE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f '-100.56' -t /dev/null
echo
# It takes -100 for an argument

echo '--- RATE FLOOR TOO BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f 2.67 -t /dev/null
echo

echo '--- RATE FLOOR E ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f 1.9e-14 -t /dev/null
echo

echo '--- RATE FLOOR INT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f 0 -t /dev/null
echo

echo '--- RATE FLOOR NAN ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f NaN -t /dev/null
echo

echo '--- RATE FLOOR INF ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f +inf -t /dev/null
echo

echo '--- RATE FLOOR STRING ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -f 'spam' -t /dev/null
echo

echo '--- THREADS NUM NEGATIVE ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ '-100' -t /dev/null
echo
# It takes -100 for an argument

echo '--- THREADS NUM A LITTLE BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ 13 -t /dev/null
echo

echo '--- THREADS NUM TOO BIG ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ 800 -t /dev/null
echo

echo '--- THREADS NUM FLOAT ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ 64.9e-10 -t /dev/null
echo

echo '--- THREADS NUM STRING ---'
python3 ${TestScript} -1 ${ScriptDir}/gzipped_normal_R1.fastq.gz -p '{"foo": "CTGTCTCTTATACAC", "bar": "CCGAAAACACG", "baz": "TCGTCGGG"}' -@ 'spam' -t /dev/null
echo

echo '--- NORMAL TEST SINGLE-END TSV EXAMPLE ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "fivezerotwo": "CTCTCTAT", "sevenzerotwo": "CTAGTACG", "sevenzeroone": "TCGCCTTA" }' -k 15 -t ${ScriptDir}/standard_test_stats_se.tsv
echo

echo '--- NORMAL TEST PAIRED-END TSV EXAMPLE ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -2 ${ScriptDir}/standard_test_R2.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "fivezerotwo": "CTCTCTAT", "sevenzerotwo": "CTAGTACG", "sevenzeroone": "TCGCCTTA" }' -k 15 -t ${ScriptDir}/standard_test_stats_pe.tsv
echo

echo '--- NORMAL TEST SINGLE-END JSON EXAMPLE ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "fivezerotwo": "CTCTCTAT", "sevenzerotwo": "CTAGTACG", "sevenzeroone": "TCGCCTTA" }' -k 15 -m 100 -t /dev/null -j ${ScriptDir}/standard_test_raw_se.json.gz
echo

echo '--- NORMAL TEST PAIRED-END JSON EXAMPLE ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -2 ${ScriptDir}/standard_test_R2.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "fivezerotwo": "CTCTCTAT", "sevenzerotwo": "CTAGTACG", "sevenzeroone": "TCGCCTTA" }' -k 15 -m 100 -t /dev/null -j ${ScriptDir}/standard_test_raw_pe.json.gz
echo

echo '--- NORMAL TEST SINGLE-END JSON EXAMPLE, LEVENSHTEIN ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "fivezerotwo": "CTCTCTAT", "sevenzerotwo": "CTAGTACG", "sevenzeroone": "TCGCCTTA" }' -k 15 -m 100 -t /dev/null -j ${ScriptDir}/standard_test_raw_levenshtein_se.json.gz -l
echo

echo '--- NORMAL TEST PAIRED-END JSON EXAMPLE, LEVENSHTEIN ---'
python3 ${TestScript} -1 ${ScriptDir}/standard_test_R1.fastq.gz -2 ${ScriptDir}/standard_test_R2.fastq.gz -p '{"oligme": "CTGTCTCTTATACACATCT", "oligb": "CCGAGCCCACGAGAC", "oliga": "TCGTCGGCAGCGTC", "fivezerotwo": "CTCTCTAT", "sevenzerotwo": "CTAGTACG", "sevenzeroone": "TCGCCTTA" }' -k 15 -m 100 -t /dev/null -j ${ScriptDir}/standard_test_raw_levenshtein_pe.json.gz -l
echo
