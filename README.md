[![PyPI](https://img.shields.io/pypi/v/simplesam.svg?)](https://pypi.python.org/pypi/simplesam)
[![Build Status](https://travis-ci.org/mdshw5/simplesam.svg?branch=master)](https://travis-ci.org/mdshw5/simplesam)
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
>>> x.gapped('seq')
'AACCCTAACCCTAACCCTAACCCTAACCCTACCCCTAACCCTACCCCTCC'
>>> x.flag
35
>>> x.mapped
True
>>> x.duplicate
False
>>> x.secondary
False
>>> x.tags
{'H1': 0, 'UQ': 33, 'RG': 'SRR011051', 'H0': 0, 'MF': 130, 'Aq': 25, 'NM': 2}
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
