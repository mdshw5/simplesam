# pylint: disable=E1101, dangerous-default-value
"""
Classes to handle alignments in the SAM format.

Reader -> Sam
Writer <- Sam
"""
try:
    from collections import OrderedDict
except ImportError: #python 2.6
    from ordereddict import OrderedDict

import sys
import os
from itertools import groupby
from subprocess import Popen, PIPE
from io import TextIOWrapper
import re
from six import PY3, string_types

try:
    from multiprocessing.dummy.connection import Connection
except ImportError: #python2
    from _multiprocessing import Connection

__version__ = '0.0.4'

class DefaultOrderedDict(OrderedDict):
    def __init__(self, default, items=[]):
        super(DefaultOrderedDict, self).__init__(items)
        self._default = default

    def __missing__(self, key):
        self[key] = value = self._default()
        return value


class GenomicOrder(object):
    def __gt__(self, other):
        if self.rname != other.rname:
            return self.rname > other.rname
        elif self.pos != other.pos:
            return self.pos > other.pos
        else:
            return str(self) > str(other)

    def __lt__(self, other):
        if self.rname != other.rname:
            return self.rname < other.rname
        elif self.pos != other.pos:
            return self.pos < other.pos
        else:
            return str(self) < str(other)

    def __eq__(self, other):
        if self.rname == other.rname and self.pos == other.pos:
            return True
        else:
            return False


class Reader(object):
    """ Read SAM/BAM format file using iterator. """
    def __init__(self, f, regions=False, kind=None):
        ext = None
        self.spool = None  # use this to catch alignment during reader scraping
        try:
            _, ext = os.path.splitext(f.name)
            if f.name == '<stdin>':  # stdin stream
                self._sam_init(f)
            elif ext == '.bam' or kind == 'bam':
                self._bam_init(f, regions)
            elif ext == '.sam' or kind == 'sam':
                self._sam_init(f)
            else:
                self._sam_init(f)
            if (regions and ext != '.bam' and kind is None) or (regions and kind is not None and kind != 'bam'):
                self.__exit__()
                raise ValueError("Region support requires bam file.")
        except AttributeError:
            if isinstance(f, Connection):
                self._pipe_init(f)
            else:
                self._sam_init(f)

    def _pipe_init(self, f):
        """ Initializer for an object that is a pipe """
        header = []
        for line in iter(f.recv, ''):
            if line[0] == '@':
                header.append(line.rstrip('\n\r'))
            else:
                self.spool = line
                break
        self.header_as_dict(header)
        self.f = iter(f.recv, '')
        self._conn = 'pipe'

    def _sam_init(self, f):
        header = []
        self.f = f
        for line in self.f:
            if line[0] == '@':
                header.append(line.rstrip('\n\r'))
            else:
                self.spool = line
                break
        self.header_as_dict(header)
        self._conn = 'file'

    def _bam_init(self, f, regions):
        pline = ['samtools', 'view', '-H', f.name]
        try:
            p = Popen(pline, bufsize=-1, stdout=PIPE,
                      stderr=PIPE)
        except OSError:
            raise OSError('Samtools must be installed for BAM file support!\n')
        self.header_as_dict([line.decode('utf-8').rstrip('\n\r') for line in p.stdout])
        p.wait()
        if regions:
            try:
                open(''.join([f.name, '.bai']))
            except EnvironmentError:
                sys.stderr.write("BAM index not found. Attempting to index file.\n")
                index_p = Popen(['samtools', 'index', f.name], stdout=PIPE, stderr=PIPE)

                _, err = index_p.communicate()
                if index_p.returncode > 0 or re.search("fail", str(err)):
                    raise OSError("Indexing failed. Is the BAM file sorted?\n")
                else:
                    sys.stderr.write("Index created successfully.\n")
            pline = ['samtools', 'view', f.name, regions]
        else:
            pline = ['samtools', 'view', f.name]
        self.p = Popen(pline, bufsize=-1, stdout=PIPE,
                  stderr=PIPE)
        if PY3:
            self.f = TextIOWrapper(self.p.stdout)
        else:
            self.f = self.p.stdout

        self._conn = 'proc'

    def next(self):
        try:
            if self.spool:  # this will be the first alignment in a SAM file or stream
                line = self.spool.rstrip('\n\r')
                self.spool = None
            else:
                line = next(self.f).rstrip('\n\r')
            if line == '':
                raise StopIteration
            fields = line.split('\t')
            required = fields[:11]
            tags = fields[11:]
            return Sam(*required, tags=tags)
        except StopIteration:
            raise StopIteration

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self

    def subsample(self, n):
        """ Draws every nth read from self. Returns Sam. """
        for i, line in enumerate(self.f):
            if i % n == 0:
                fields = line.split('\t')
                required = fields[:11]
                tags = fields[11:]
                yield Sam(*required, tags=tags)

    def header_as_dict(self, header):
        """ Parse the header list and return a nested dictionary. """
        self.header = DefaultOrderedDict(OrderedDict)
        for line in header:
            line = line.split('\t')
            key, fields = (line[0], line[1:])
            try:
                self.header[key][fields[0]] = fields[1:]
            except IndexError:
                self.header[key][fields[0]] = ['']

    def merge_header(self, header):
        for key, values in header.items():
            for k, v in values.items():
                self.header[key][k] = v

    @property
    def seqs(self):
        """ Return just the sequence names from the @SQ library """
        for key in self.header['@SQ'].keys():
            yield key.split(':')[1]

    def tile_genome(self, width):
        """ Return UCSC-style regions tiling 'width' """
        assert isinstance(width, int)
        for k, v in self.header['@SQ'].items():
            rname = k.split(':')[1]
            seqlength = v[0].split(':')[1]
            for region in tile_region(rname, 1, int(seqlength), width):
                yield region

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self._conn == 'file':
            self.f.close()
        if self._conn == 'proc':
            self.f.close()
            self.p.terminate()


class Writer(object):
    """ Write SAM/BAM format file """
    def __init__(self, f, header=None):
        try:
            _, ext = os.path.splitext(f.name)
            if ext == '.bam':
                raise NotImplementedError('Bam writing support is not implemented.\n')
        except AttributeError:
            pass
        self.file = f
        if header is not None:
            self.header = DefaultOrderedDict(OrderedDict)
            self.merge_header(header)
        else:
            self.header = DefaultOrderedDict(OrderedDict)
            self.header['@HD']['VN:1.0'] = ['SO:unknown']
        self.header_dict_format()

    def merge_header(self, header):
        for key, values in header.items():
            for k, v in values.items():
                self.header[key][k] = v

    def header_dict_format(self):
        for key, value in self.header.items():
            for k, v in value.items():
                tags = '\t'.join(v)
                self.file.write('{key}\t{k}\t{tags}\n'.format(**locals()))

    def write(self, sread):
        self.file.write(str(sread))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class Sam(GenomicOrder):
    """ Store fields in each line of a SAM file, provided as a tuple. """
    cigar_nochange = set(("M", "N", "EQ", "X", "P"))
    cigar_len = set(("M", "D", "N", "EQ", "X", "P"))

    def __init__(self, qname='@', flag=4, rname='*', pos=0, mapq=255, cigar='*', rnext='*', pnext=0, tlen=0, seq='*', qual='*', tags=[]):
        self.qname = qname
        self.flag = int(flag)
        self.rname = rname
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = int(pnext)
        self.tlen = int(tlen)
        self.seq = seq
        self.qual = qual
        self._tags = tags
        self._cache = dict()
        self.tags = None

    def __str__(self):
        if self.tags:
            tag_fields = '\t'.join([encode_tag(tag, self.tags[tag]) for tag in sorted(self.tags.keys())])
        else:
            tag_fields = '\t'.join(self._tags)
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.format(self.qname,
                                                                                     str(self.flag),
                                                                                     self.rname,
                                                                                     str(self.pos),
                                                                                     str(self.mapq),
                                                                                     self.cigar,
                                                                                     self.rnext,
                                                                                     str(self.pnext),
                                                                                     str(self.tlen),
                                                                                     self.seq,
                                                                                     self.qual,
                                                                                     tag_fields)
    def __repr__(self):
        return "Sam({0}:{1}:{2})".format(self.rname, self.pos, self.qname)

    def __len__(self):
        return sum(c[0] for c in self.cigars if c[1] in self.cigar_len)

    def __getitem__(self, key):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        return self.tags[key]

    def __setitem__(self, key, value):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        self.tags[key] = value

    def index_of(self, pos):
        """ Return the relative index from genomic position """
        i = pos - self.pos
        if i >= 0:
            return i
        else:
            raise IndexError("Position {0:n} not in {1}.".format(pos, self.qname))

    def get(self, key, default_value):
        try:
            return self[key]
        except KeyError:
            return default_value

    def cigar_split(self):
        """ CIGAR grouping function modified from:
        https://github.com/brentp/bwa-meth
        """
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    def parse_md(self):
        """ Return the ungapped reference sequence from the MD tag, if present.
        """
        try:
            return self._cache['parse_md']
        except KeyError:
            pass
        try:
            md = self['MD']
        except KeyError:
            raise KeyError('MD tag not found in SAM record.')
        ref_seq = list(self.gapped('seq'))
        md_match = re.findall(r"([0-9]+)\^?([A-Z]+)?", md)
        ref_seq_i = 0
        for i, b in md_match:
            ref_seq_i += int(i)
            for mismatch in b:
                ref_seq[ref_seq_i] = mismatch
                ref_seq_i += 1
        self._cache['parse_md'] = ref_seq
        return ref_seq

    @property
    def cigars(self):
        try:
            return self._cache['cigars']
        except KeyError:
            self._cache['cigars'] = tuple(self.cigar_split())
            return self._cache['cigars']

    @property
    def paired(self):
        return bool(self.flag & 0x2)

    @property
    def mapped(self):
        return not (self.flag & 0x4)

    @property
    def secondary(self):
        return bool(self.flag & 0x100)

    @property
    def reverse(self):
        return bool(self.flag & 0x10)

    @property
    def passing(self):
        return not bool(self.flag & 0x200)

    @property
    def duplicate(self):
        return bool(self.flag & 0x400)

    def gapped(self, attr, gap_char='-'):
        """ Return string with all deletions wrt reference
         represented as gaps '-' and all insertions wrt reference
         removed.
        """
        try:
            ungapped = getattr(self, attr)
        except AttributeError:
            ungapped = self[attr]  # get dictionary key (tag) if attribute is missing
        gapped = []
        i = 0
        for n, t in self.cigars:
            if t in self.cigar_nochange:
                gapped.extend(ungapped[i:i + n])
                i += n
            elif t in ("D",):
                gapped.extend([gap_char] * n)
            elif t in ("I", "S"):
                i += n
            elif t in ("H",):
                pass
        return ''.join(gapped)

    @property
    def coords(self):
        return range(self.pos, self.pos + len(self))

    @property
    def safename(self):
        """Return self.qname without paired-end identifier if it exists"""
        if self.qname[-2] == '/':
            return self.qname[:-2]
        else:
            return self.qname


def parse_sam_tags(tagfields):
    """ Return a dictionary containing the tags """
    return dict([(tag, data) for tag, dtype, data in [decode_tag(x) for x in tagfields]])


def encode_tag(tag, data):
    """ Write a SAM tag in the format ``TAG:TYPE:data``
    >>> encode_tag('YM', '#""9O"1@!J')
    'YM:Z:#""9O"1@!J'
    """
    if isinstance(data, string_types):
        data_type = 'Z'
    elif isinstance(data, int):
        data_type = 'i'
    elif isinstance(data, float):
        data_type = 'f'
    else:
        raise NotImplementedError("Data {0} cannot be encoded as string, integer, or float tag.".format(data))
    value = ':'.join((tag, data_type, str(data)))
    return value


def decode_tag(tag_string):
    """ Parse a SAM format tag to a (TAG, TYPE, data) tuple.

    TYPE in A, i, f, Z, H, B

    >>> decode_tag('YM:Z:#""9O"1@!J')
    ('YM', 'Z', '#""9O"1@!J')
    >>> decode_tag('XS:i:5')
    ('XS', 'i', 5)
    >>> decode_tag('XF:f:100.5')
    ('XF', 'f', 100.5)
    """
    try:
        tag, data_type, data = tag_string.split(':')
    except ValueError:
        match = re.match(r'([A-Z]{2}):([iZfHB]):(\S+)', tag_string)
        tag = match.group(1)
        data_type = match.group(2)
        data = match.group(3)
    if data_type == 'i':
        return (tag, data_type, int(data))
    elif data_type == 'Z':
        return (tag, data_type, data)
    elif data_type == 'f':
        return (tag, data_type, float(data))
    elif data_type == 'A':  # this is just a special case of a character
        return (tag, data_type, data)
    elif data_type == 'H':
        raise NotImplementedError("Hex array SAM tags are currently not parsed.")
    elif data_type == 'B':
        raise NotImplementedError("Byte array SAM tags are currently not parsed.")
    else:
        raise NotImplementedError("Tag {0} cannot be parsed.".format(tag_string))


def tile_region(rname, start, end, step):
    """ Make non-overlapping tiled windows from the specified region
    >>> list(tile_region('chr1', 1, 250, 100))
    ['chr1:1-100', 'chr1:101-200', 'chr1:201-250']
    >>> list(tile_region('chr1', 1, 200, 100))
    ['chr1:1-100', 'chr1:101-200']
    """
    while start + step <= end:
        yield '%s:%d-%d' % (rname, start, start + step - 1)
        start += step
    if start < end:
        yield '%s:%d-%d' % (rname, start, end)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
