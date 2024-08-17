import sys

from typing import List

import numpy as np

from jcvi.apps.base import OptionParser, logger
from jcvi.formats.bed import Bed, BedLine
from jcvi.formats.fasta import Fasta


def main(args):
    p = OptionParser(main.__doc__)
    _, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, fastafile = args
    fasta = Fasta(fastafile)
    bed = Bed(bedfile)
    for read, sb in bed.sub_beds():
        get_breakpoint(read, sb, fasta)


def get_breakpoint(read: str, sb: List[BedLine], fasta: Fasta):
    """
    Given a read and its sub-beds, determine the breakpoint
    """
    ra, rb, read_id = read.split("_", 2)
    sids = [b.accn.split(":", 1)[0] for b in sb]

    n = len(sids)

    prefix_a = np.zeros(n, dtype="int")
    prefix_b = np.zeros(n, dtype="int")
    for i in range(n):
        if i == 0:
            prefix_a[i] = sids[i] == ra
            prefix_b[i] = sids[i] == rb
        else:
            prefix_a[i] = prefix_a[i - 1] + (sids[i] == ra)
            prefix_b[i] = prefix_b[i - 1] + (sids[i] == rb)
    suffix_a = [0] * n
    suffix_b = [0] * n
    for i in range(n - 1, 0, -1):
        suffix_a[i - 1] = suffix_a[i] + (sids[i] == ra)
        suffix_b[i - 1] = suffix_b[i] + (sids[i] == rb)
    ab = prefix_a + suffix_b
    ba = prefix_b + suffix_a

    if ab.max() > ba.max():
        idx = np.argmax(ab)
    else:
        idx = np.argmax(ba)
        ra, rb = rb, ra

    left = sb[idx].end
    print(sb[idx])
    if idx == n - 1:
        logger.info("Breakpoint at the end of `%s`. Skipped", read)
        return
    right = sb[idx + 1].start
    print(sb[idx + 1])
    mid = (left + right) // 2
    seq = fasta[read_id]
    end = len(seq)
    print(mid, f"{ra}=0-{mid}", f"{rb}={mid}-{end}")


if __name__ == "__main__":
    main(sys.argv[1:])
