[![PyPI](https://img.shields.io/pypi/v/simplesam.svg?)](https://pypi.python.org/pypi/simplesam)

# Simple SAM parsing
Requiring no external dependencies (except a samtools installation for BAM reading)

# Installation
`pip install simplesam`

# Usage

```python
>>> from simplesam import Reader, Writer

# can also read BAM
>>> in_file = open('data/NA18510.sam', 'r')
>>> in_sam = Reader(in_file)

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

>>> out_file = open('test.sam', 'w')
>>> out_sam = Writer(out_file, in_sam.header)
>>> out_sam.write(x)
>>> out_sam.close()
```
