#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("samfile", help="input sam or bam file")
# parser.add_argument("",help="",default=None)
parser.add_argument("-d", action="store_true", default=False)
parser.add_argument(
    "-m", "--matches", help="min num of matches to return", type=int, default=0
)
parser.add_argument("-t", "--threads", help="threads", type=int, default=8)
parser.add_argument("--header", action="store_true", default=False)
args = parser.parse_args()
DEBUG = args.d

import glob
import os
import sys
import re
import pysam


def makeHeader():
    rtn = ""
    rtn += "{}\t{}\t{}\t{}\t".format(
        "query_name", "query_start", "query_end", "query_length"
    )
    rtn += "{}\t{}\t{}\t".format("reference_name", "reference_start", "reference_end")
    rtn += "{}\t{}\t{}\t".format("perID_by_matches", "perID_by_events", "perID_by_all")
    rtn += "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
        "matches",
        "mismatches",
        "insertions",
        "deletions",
        "insertion_events",
        "deletion_events",
        "strand",
    )
    return rtn


def perId(matches, mismatch, ins, dele, insEvent, delEvent):
    if (matches + mismatch) == 0:
        return (0, 0, 0)

    bymatches = (100.0 * matches) / (matches + mismatch)
    byevents = (100.0 * matches) / (matches + mismatch + insEvent + delEvent)
    byall = (100.0 * matches) / (matches + mismatch + ins + dele)
    return (bymatches, byevents, byall)


# /(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)+/
# :[0-9]+ represents an identical block,
# -ata represents a deltion,i
# +gtc an insertion and
# *at indicates reference base a is substituted with a query base t.
# It is similar to the MD SAM tag but is standalone and easier to parse.


def formatRead(read):
    rtn = ""
    rtn += "{}\t{}\t{}\t{}\t".format(
        read.query_name,
        read.query_alignment_start,
        read.query_alignment_end,
        read.infer_query_length(),
    )
    rtn += "{}\t{}\t{}\t".format(
        read.reference_name, read.reference_start, read.reference_end
    )
    if read.flag == 4:
        return ""

    tags = read.get_tags()

    if "cs" in tags:
        cs = read.get_tag("cs")
        pattern = ":[0-9]+"
        match = re.findall(pattern, cs)
        if match is not None:
            match = [int(x[1:]) for x in match]
            matches = sum(match)

        pattern = "-[a-z]+"
        match = re.findall(pattern, cs)
        if match is not None:
            match = [len((x[1:])) for x in match]
            delEvent = len(match)
            dele = sum(match)

        pattern = "\+[a-z]+"
        match = re.findall(pattern, cs)
        if match is not None:
            match = [len((x[1:])) for x in match]
            insEvent = len(match)
            ins = sum(match)

        mismatch = cs.count("*")
    else:
        counts, events = read.get_cigar_stats()
        if "MD" in tags:
            mismatch = sum(
                [
                    bp.islower()
                    for _r, _q, bp in read.get_aligned_pairs(with_seq=True)
                    if bp is not None
                ]
            )
            matches = counts[0] - mismatch
        else:
            matches = counts[7]
            mismatch = counts[8]
        ins = counts[1]
        dele = counts[2]
        insEvent = events[1]
        delEvent = events[2]

    strand = "+"
    if read.is_reverse:
        strand = "-"
    bymatches, byevents, byall = perId(matches, mismatch, ins, dele, insEvent, delEvent)
    rtn += "{}\t{}\t{}\t".format(bymatches, byevents, byall)
    rtn += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        matches, mismatch, ins, dele, insEvent, delEvent, strand
    )
    if matches >= args.matches:
        return rtn
    return ""


if args.header:
    print(makeHeader())
samfile = pysam.AlignmentFile(args.samfile, threads=args.threads)
for read in samfile.fetch(until_eof=True):
    sys.stdout.write(formatRead(read))
samfile.close()
