import io
import os
import pandas as pd

def read_vcf(path):
    
    """
    source code: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
    """
    
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

    return df