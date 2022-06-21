#!/usr/bin/env python

import sys
import argparse
import pkg_resources
from collections import deque
from collections import Counter
from collections import OrderedDict

class Extractor(OrderedDict):

	def __missing__(self, pos):
		value = (self.sam, deque(), deque())
		self[pos] = value
		return value

	def update(self, sam, min_qual):
		self.sam = sam
		g = sam.pos
		for i, (seq, qual) in enumerate(zip(sam.gapped('seq', '-'), sam.gapped('qual', '~'))):
			if qual > min_qual:
				self[g + i][1].append(seq)
				self[g + i][2].append(qual)

	def clear(self, stop=None):
		for pos in sorted(k for k in self.keys() if stop is not None and k < stop):
			sam, seq, qual = self.pop(pos)
			ref = sam.parse_md()[sam.index_of(pos)]
			try:
				s, q = zip(*[(s, q) for s, q in zip(seq, qual) if s != ref])
			except ValueError:
				s, q = ([], [])
			yield (sam.rname, str(pos), ref, str(len(seq) - len(s)), str(len(s)), ''.join(s), ''.join(q))


def pileup(args):
	from simplesam import Reader

	e = Extractor()
	bases = ('A', 'C', 'T', 'G', 'N', '-')
	stats = dict([('A', Counter()), ('C', Counter()), ('T', Counter()), ('G', Counter()), ('N', Counter()), ('-', Counter())])
	with Reader(args.bam) as bam:
		for read in bam:
			if not read.mapped:
				continue
			if read.duplicate:
				continue
			if read.secondary:
				continue
			for line in e.clear(stop=read.pos):
				if args.counts:
					args.pileup.write('\t'.join(line) + '\t' + '\t'.join([str(line[5].count(c)) for c in bases]) + '\n')
				else:
					args.pileup.write('\t'.join(line) + '\n')
				if args.stats:
					stats[line[2]][line[2]] += int(line[3])
					stats[line[2]].update(line[5])
			e.update(read, '!')

		for line in e.clear():
			if args.counts:
				args.pileup.write('\t'.join(line) + '\t' + '\t'.join([str(line[5].count(c)) for c in bases]) + '\n')
			else:
				args.pileup.write('\t'.join(line) + '\n')
			if args.stats:
				stats[line[2]][line[2]] += int(line[3])
				stats[line[2]].update(line[5])

		if args.stats:
			for ref, counts in sorted(stats.items()):
				for base in bases:
					ref_count = counts[ref]
					base_count = counts[base]
					try:
						percent_base = base_count / sum(counts.values())
					except ZeroDivisionError:
						percent_base = 0.
					args.stats.write('{ref}\t{base}\t{ref_count}\t{base_count}\t{percent_base:.4%}\n'.format(**locals()))


def main():
	parser = argparse.ArgumentParser(prog='pileup', description="generate a simple pileup-like file from a sorted/indexed BAM file")
	parser.add_argument('--version', action='version', version="%(prog)s version {0}".format(pkg_resources.get_distribution("simplesam").version))

	parser.add_argument('bam', type=argparse.FileType('r'), help="sorted/indexed BAM file ")
	parser.add_argument('pileup', type=argparse.FileType('w'), help="pileup output file")
	parser.add_argument('-c', '--counts', action='store_true', default=False, help="display counts for A/C/T/G/N/- separately (default: %(default)s)")
	parser.add_argument('-i', '--stats', type=argparse.FileType('w'), help="tabulate mismatches to output file")
	parser.set_defaults(func=pileup)

	args = parser.parse_args()
	args.func(args)

if __name__ == "__main__":
	main()