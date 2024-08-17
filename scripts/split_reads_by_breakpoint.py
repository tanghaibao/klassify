import sys

from typing import List

import numpy as np

from jcvi.apps.base import OptionParser
from jcvi.formats.bed import Bed, BedLine


def main(args):
    p = OptionParser(main.__doc__)
    _, args = p.parse_args(args)

    bed = Bed(args[0])
    for read, sb in bed.sub_beds():
        get_breakpoint(read, sb)


def get_breakpoint(read: str, sb: List[BedLine]):
    """
    Given a read and its sub-beds, determine the breakpoint
    """
    ra, rb = read.split("_", 2)[:2]
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
    if idx == n - 1:
        logger.info(f"Read {read} is fully contained in {sb[idx]}")
        mid = left
    else:
        right = sb[idx + 1].start
        mid = (left + right) // 2
    print(sb[idx])
    print(sb[idx + 1])
    print(mid, f"{ra}=0-{mid}", f"{rb}={mid}-end")


if __name__ == "__main__":
    main(sys.argv[1:])
