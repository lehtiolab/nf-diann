#!/usr/bin/env python3

# Precursor table output, from parquet to tsv, could possibly run directly
# inside DIA-NN process

import sys

from pyarrow import compute as pc
from pyarrow import parquet as pq
from pyarrow import csv as pcsv

precursors = pq.read_table(sys.argv[1])

if sys.argv[2] == 'trypsin':
    miscl_re = '[KR][^P]'
elif sys.argv[2] == 'trypsinp':
    miscl_re = '[KR][A-Z]'

precursors = precursors.append_column('missed_cleavages', pc.count_substring_regex(precursors['Stripped.Sequence'], miscl_re))

write_opt = pcsv.WriteOptions(delimiter='\t', quoting_style='none', quoting_header='none')
pcsv.write_csv(precursors, 'precursors.txt', write_opt)

