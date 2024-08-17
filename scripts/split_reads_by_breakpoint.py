import sys

from collections import Counter
from io import TextIOWrapper
from typing import List

import numpy as np

from jcvi.apps.base import OptionParser, logger
from jcvi.formats.bed import Bed, BedLine
from jcvi.formats.fasta import Fasta


SUCCESS, FAIL = "SUCCESS", "FAIL"


def main(args):
    p = OptionParser(main.__doc__)
    _, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, fastafile = args
    fasta = Fasta(fastafile)
    bed = Bed(bedfile)
    new_fastafile = bedfile.replace(".gz", "").rsplit(".", 1)[0] + ".split.fasta"
    counter = Counter()
    with open(new_fastafile, "w", encoding="utf-8") as new_fasta:
        for read, sb in bed.sub_beds():
            ret = get_breakpoint(read, sb, fasta, new_fasta)
            counter[ret] += 1
    logger.info("Summary: %s", counter)
    logger.info("Split fasta file written to `{}`".format(new_fastafile))


def get_breakpoint(
    read: str, sb: List[BedLine], fasta: Fasta, new_fasta: TextIOWrapper
) -> str:
    """
    Given a read and its sub-beds, determine the breakpoint
    """
    ra, rb, orig_read = read.split("_", 2)
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
        return FAIL
    right = sb[idx + 1].start
    print(sb[idx + 1])
    mid = (left + right) // 2
    seq = fasta[read]
    end = len(seq)
    left_subread_id = f"{orig_read}|{ra}|0-{mid}"
    left_subread = seq.seq[:mid]
    right_subread_id = f"{orig_read}|{rb}|{mid}-{end}"
    right_subread = seq.seq[mid:]
    print(left_subread_id, right_subread_id)
    new_fasta.write(f">{left_subread_id}\n{left_subread}\n")
    new_fasta.write(f">{right_subread_id}\n{right_subread}\n")
    return SUCCESS


if __name__ == "__main__":
    main(sys.argv[1:])
