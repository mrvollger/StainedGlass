#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import os
import sys
import argparse
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="input bed file")
    parser.add_argument(
        "--outputs", nargs="+", help="list of output files", required=True
    )
    args = parser.parse_args()

    bed = open(args.infile).readlines()
    N_WINDOWS = len(bed)
    N = len(args.outputs)
    splits = np.array_split(np.arange(N_WINDOWS), N, axis=0)
    for i, split in enumerate(splits):
        f = open(args.outputs[i], "w+")
        for pos in split:
            f.write(bed[pos])
        f.close()
