#!/usr/bin/env python3
'''
Make a samples.json file with sample names and file names.
'''

import json
from glob import glob

# Change this line to match your filenames.
fastqs = glob('/data/srlab/slowikow/src/snakefiles/data/fastq/*.fastq.gz')
FILES = {}

# Change this line to extract a sample name from each filename.
SAMPLES = [fastq.split('/')[-1].split('.')[0] for fastq in fastqs]

for sample in SAMPLES:
    # Change 'R1' and 'R2' to match the way your mate pairs are marked.
    mate1 = lambda fastq: sample in fastq and 'R1' in fastq
    mate2 = lambda fastq: sample in fastq and 'R2' in fastq
    FILES[sample] = {}
    FILES[sample]['R1'] = sorted(filter(mate1, fastqs))
    FILES[sample]['R2'] = sorted(filter(mate2, fastqs))

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)
