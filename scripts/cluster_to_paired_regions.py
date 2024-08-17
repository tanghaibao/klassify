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

    read_to_regions = defaultdict(list)
    read_to_subreads = defaultdict(set)
    read_to_match = {}
    for mb in merged:
        if len(mb) < MIN_READ_SUPPORT:
            continue
        mb_start = int(np.median([b.start for b in mb]))
        mb_end = int(np.median([b.end for b in mb]))
        region_name = f"{mb[0].seqid}:{mb_start}-{mb_end}"
        for b in mb:
            read_name = b.accn.split("|", 1)[0]
            read_to_regions[read_name].append(region_name)
            read_to_subreads[read_name].add((b.accn, b.strand))
            read_to_match[b.accn] = f"{b.seqid}:{b.start}-{b.end}:{b.strand}"

    pair_to_reads = defaultdict(list)
    for read, regions in read_to_regions.items():
        if len(regions) == 2:
            pair_to_reads[tuple(sorted(regions))].append(read)

    G = nx.Graph()
    for pair, reads in pair_to_reads.items():
        G.add_edge(*pair, weight=len(reads))
    valid_pairs = sorted(tuple(sorted(x)) for x in nx.max_weight_matching(G))

    filtered_pair_to_reads = {}
    for pair in valid_pairs:
        n_reads = len(pair_to_reads[pair])
        if n_reads < MIN_READ_SUPPORT:
            continue
        counter = Counter()
        filtered_reads = []
        for read in pair_to_reads[pair]:
            (fa, _), (fb, _) = read_to_subreads[read]
            # Check if the subread is out of order
            fa_read, fa_seqid, fa_range = fa.split("|")
            fb_read, fb_seqid, fb_range = fb.split("|")
            assert fa_read == fb_read
            if fb_range.startswith("0-"):
                assert not fa_range.endswith("0-")
                fa, fb = fb, fa
                fa_seqid, fb_seqid = fb_seqid, fa_seqid
            counter[(fa_seqid, fb_seqid)] += 1
            filtered_reads.append((fa, fb))
        (ra_reordered, rb_reordered), _ = tuple(counter.most_common(1)[0])
        ra, rb = pair
        if not ra.startswith(ra_reordered):
            ra, rb = rb, ra
        assert ra.startswith(ra_reordered) and rb.startswith(rb_reordered)
        filtered_pair_to_reads[(ra, rb)] = filtered_reads

    # Finally, write out the paired regions
    header = [
        "Crossover ID",
        "Left",
        "Right",
        "Read Count",
        "Read ID",
        "Read Left",
        "Read Left Match",
        "Read Right",
        "Read Right Match",
    ]
    print("\t".join(header))
    for cid, ((ra, rb), reads) in enumerate(sorted(filtered_pair_to_reads.items())):
        cid += 1
        for i, read in enumerate(reads):
            row = [cid, "", "", ""] if i > 0 else [cid, ra, rb, len(reads)]
            row += [
                read[0].split("|")[0],
                read[0].split("|")[2],
                read_to_match[read[0]],
                read[1].split("|")[2],
                read_to_match[read[1]],
            ]
            print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    main(sys.argv[1:])
