import sys

from collections import Counter, defaultdict
from typing import List

import numpy as np
import networkx as nx

from jcvi.formats.bed import Bed
from jcvi.apps.base import OptionParser, logger


SEQID_MISMATCH, SECONDARY_MATCH = "SEQID_MISMATCH", "SECONDARY_MATCH"
MIN_READ_SUPPORT = 3


def main(args: List[str]):
    p = OptionParser(main.__doc__)
    _, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(not p.print_help())

    (bedfile,) = args
    bed = Bed(bedfile)
    filtered = []
    seen = {}
    counter = Counter()
    # Store the primary match length
    for b in bed:
        length = b.end - b.start
        seen[b.accn] = max(seen.get(b.accn, 0), length)
    # Second pass to filter out the secondary matches
    for b in bed:
        expected_seqid = b.accn.split("|")[1]
        if expected_seqid != b.seqid:
            counter[SEQID_MISMATCH] += 1
            continue
        if seen[b.accn] != b.end - b.start:
            counter[SECONDARY_MATCH] += 1
            continue
        filtered.append(b)
    logger.info(counter)
    logger.info("Total filtered: %d", len(filtered))

    # Merge the regions that are clustered
    merged = []
    for b in filtered:
        if merged and (
            merged[-1][-1].seqid == b.seqid and merged[-1][-1].end >= b.start
        ):
            merged[-1].append(b)
        else:
            merged.append([b])

    region_to_reads = defaultdict(list)
    read_to_regions = defaultdict(list)
    for mb in merged:
        if len(mb) < MIN_READ_SUPPORT:
            continue
        mb_start = int(np.median([b.start for b in mb]))
        mb_end = int(np.median([b.end for b in mb]))
        region_name = f"{mb[0].seqid}:{mb_start}-{mb_end}"
        for b in mb:
            region_to_reads[region_name].append(b.accn)
            read_name = b.accn.split("|", 1)[0]
            read_to_regions[read_name].append(region_name)

    pair_to_reads = defaultdict(list)
    for read, regions in read_to_regions.items():
        if len(regions) == 2:
            pair_to_reads[tuple(sorted(regions))].append(read)

    G = nx.Graph()
    for pair, reads in pair_to_reads.items():
        G.add_edge(*pair, weight=len(reads))
    valid_pairs = sorted(tuple(sorted(x)) for x in nx.max_weight_matching(G))

    for pair in valid_pairs:
        print(pair, len(pair_to_reads[pair]))


if __name__ == "__main__":
    main(sys.argv[1:])
