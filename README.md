[![PyPI](https://img.shields.io/pypi/v/simplesam.svg?)](https://pypi.python.org/pypi/simplesam)
[![Tests](https://github.com/mdshw5/simplesam/actions/workflows/tests.yml/badge.svg)](https://github.com/mdshw5/simplesam/actions/workflows/tests.yml)
[![Package Builds](https://github.com/mdshw5/simplesam/actions/workflows/deploy.yml/badge.svg)](https://github.com/mdshw5/simplesam/actions/workflows/deploy.yml)
[![Documentation Status](https://readthedocs.org/projects/simplesam/badge/?version=latest)](http://simplesam.readthedocs.io/en/latest/?badge=latest)

# Simple SAM parsing
Requiring no external dependencies (except a samtools installation for BAM reading)

# Installation
`pip install simplesam`

# Usage
For complete module documentation visit [ReadTheDocs](http://simplesam.readthedocs.io/en/latest/).

## Quickstart

```python
>>> from simplesam import Reader, Writer
```

Read from SAM/BAM files
```python
# can also read BAM
>>> in_file = open('data/NA18510.sam', 'r')
>>> in_sam = Reader(in_file)
```

Access alignments using an iterator interface
```python
>>> x = next(in_sam)
>>> type(x)
<class 'simplesam.Sam'>
>>> x
Sam(1:2:SRR011051.1022326)
>>> x.qname
'SRR011051.1022326'
>>> x.rname
'1'
>>> x.pos
2
>>> x.seq
'AACCCTAACCCCTAACCCTAACCCTAACCCTACCCCTAACCCTACCCCTCC'
>>> x.qual
'?<:;;=;>;<<<>96;<;;99;<=3;4<<:(;,<;;/;57<;%6,=:,((3'
>>> x.cigar
'8M1I42M'
>>> x.cigars
((8, 'M'), (1, 'I'), (42, 'M'))
>>> x.gapped('seq')
'AACCCTAACCCTAACCCTAACCCTAACCCTACCCCTAACCCTACCCCTCC'
>>> len(x)
50
>>> x.flag
35
>>> x.mapped
True
>>> x.paired
True
>>> x.duplicate
False
>>> x.secondary
False
>>> x.coords
[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51]
>>> x.tags
{'H1': 0, 'UQ': 33, 'RG': 'SRR011051', 'H0': 0, 'MF': 130, 'Aq': 25, 'NM': 2}
>>> str(x)
'SRR011051.1022326\t35\t1\t2\t255\t8M1I42M\t*\t0\t0\tAACCCTAACCCCTAACCCTAACCCTAACCCTACCCCTAACCCTACCCCTCC\t?<:;;=;>;<<<>96;<;;99;<=3;4<<:(;,<;;/;57<;%6,=:,((3\tAq:i:25\tH0:i:0\tH1:i:0\tMF:i:130\tNM:i:2\tRG:Z:SRR011051\tUQ:i:33\n'
```

Read the SAM sequence header structure
```python
>>> from pprint import pprint
>>> pprint(in_sam.header)
{'@HD': OrderedDict([('VN:1.0', ['GO:none', 'SO:coordinate'])]),
 '@SQ': {'SN:1': ['LN:247249719'],
         'SN:2': ['LN:242951149'],
         'SN:3': ['LN:199501827'],
         'SN:4': ['LN:191273063'],
         'SN:5': ['LN:180857866'],
         'SN:6': ['LN:170899992'],
         'SN:7': ['LN:158821424'],
         'SN:8': ['LN:146274826'],
         ...
 '@RG': {'ID:SRR011049': ['PL:ILLUMINA',
                           'PU:BI.PE.080626_SL-XAN_0002_FC304CDAAXX.080630_SL-XAN_0007_FC304CDAAXX.5',
                           'LB:Solexa-5112',
                           'PI:330',
                           'SM:NA18510',
                           'CN:BI'],
          'ID:SRR011050': ['PL:ILLUMINA',
                           'PU:BI.PE.080626_SL-XAN_0002_FC304CDAAXX.080630_SL-XAN_0007_FC304CDAAXX.6',
                           'LB:Solexa-5112',
                           'PI:330',
                           'SM:NA18510',
                           'CN:BI'],
         ...}
         }
}
```

Write SAM files from `Sam` objects
```python
# Reader and Writer can also use the context handler (with: statement)
>>> out_file = open('test.sam', 'w')
>>> out_sam = Writer(out_file, in_sam.header)
>>> out_sam.write(x)
>>> out_sam.close()
```

Write SAM files from `Sam` objects to stdout (allows `samtools view` compression)
```python
>>> from sys import stdout
>>> stdout_sam = Writer(stdout, in_sam.header)
>>> stdout_sam.write(x)
>>> stdout_sam.close()
```
```bash
$ python my_script_that_uses_simplesam.py | samtools view -hbo test.bam
```

# Example scripts
An example script [`pileup.py`](https://github.com/mdshw5/simplesam/blob/master/scripts/pileup.py) is installed with this module.
This script will generate an output that is similar to `samtools pileup` with the addition of several optional columns that summarize
counts for individual nucleotides (ACTGN) and deletions with respect to the reference (-). This script leverages the `Sam.gapped()` and
`Sam.parse_md()` methods to reconstruct position-specific counts from SAM alignment records.

```bash
$ pileup.py -h
usage: pileup [-h] [--version] [-c] [-i STATS] bam pileup

generate a simple pileup-like file from a sorted/indexed BAM file

positional arguments:
  bam                   sorted/indexed BAM file
  pileup                pileup output file

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -c, --counts          display counts for A/C/T/G/N/- separately (default: False)
  -i STATS, --stats STATS
                        tabulate mismatches to output file
```