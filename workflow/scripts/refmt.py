#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import os
import sys
import argparse
import pandas as pd


# from dask.distributed import Client, progress
# import dask
# import dask.dataframe as dd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="samid table")
    parser.add_argument("outfile", help="output bed file")
    parser.add_argument("--window", help="window size", type=int, required=True)
    parser.add_argument(
        "--full",
        help="place to drop the full bed file with the extra columns",
        default=None,
    )
    parser.add_argument(
        "--ncolors", help="number of color qauntiles", type=int, default=10
    )
    parser.add_argument("--threads", help="nthreads", type=int, default=10)
    parser.add_argument("--fai", help="fai for the genome", required=True)
    parser.add_argument(
        "-d", help="store args.d as true if -d", action="store_true", default=False
    )
    parser.add_argument(
        "--one",
        help="filter for best 1-1 alignment",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()
    # client = Client(n_workers=1, threads_per_worker=args.threads, memory_limit="64GB")

    # df = dd.read_csv(args.infile, sep="\t")
    df = pd.read_csv(args.infile, sep="\t")

    if not args.one:
        # determine duplicates by these columns
        det_dup_cols = ["query_name", "reference_name"]
        # only keep the alignment with the most matching bases
        df.sort_values(by=det_dup_cols + ["matches"], inplace=True, ascending=False)
        df.drop_duplicates(subset=det_dup_cols, inplace=True)

    # extract the position from the name and add it to the start and end
    q_fix = df.query_name.str.extract(r"(.+):(\d+)-\d+", expand=True)
    df.query_name = q_fix[0]
    df.query_start = q_fix[1].astype(int)
    df.query_end = q_fix[1].astype(int) + args.window

    if not args.one:
        # extract ref coords
        r_fix = df.reference_name.str.extract(r"(.+):(\d+)-\d+", expand=True)
        df.reference_name = r_fix[0]
        df.reference_start = r_fix[1].astype(int)
        df.reference_end = r_fix[1].astype(int) + args.window

    # limit the size if we go over the chr end
    fai = pd.read_csv(args.fai, sep="\t", names=["chr", "length", "x", "y", "z"])[
        ["chr", "length"]
    ]
    df = pd.merge(df, fai, left_on="query_name", right_on="chr")
    df = pd.merge(df, fai, left_on="reference_name", right_on="chr")
    df.loc[df.query_end >= df.length_x, "length_x"] = (
        df.loc[df.query_end >= df.length_x, "length_x"] - 1
    )
    df.loc[df.reference_end >= df.length_y, "length_y"] = (
        df.loc[df.reference_end >= df.length_y, "length_y"] - 1
    )

    # columns to keep in the defualt output
    keep_cols = [
        "query_name",
        "query_start",
        "query_end",
        "reference_name",
        "reference_start",
        "reference_end",
        "perID_by_events",
        "strand",
    ]
    out = df
    if not args.one:
        out2 = out.copy()
        out2.reference_name = out.query_name
        out2.reference_start = out.query_start
        out2.reference_end = out.query_end
        out2.query_name = out.reference_name
        out2.query_start = out.reference_start
        out2.query_end = out.reference_end

        out = (
            pd.concat([out, out2])
            .sort_values(
                by=[
                    "query_name",
                    "query_start",
                    "reference_name",
                    "reference_start",
                    "perID_by_events",
                ],
                ascending=False,
            )
            .drop_duplicates(
                subset=[
                    "query_name",
                    "query_start",
                    "reference_name",
                    "reference_start",
                ]
            )
            .sort_values(
                by=[
                    "query_name",
                    "query_start",
                    "reference_name",
                    "reference_start",
                    "query_end",
                ]
            )
        )

    out["qcut"] = pd.qcut(
        out["perID_by_events"],
        min(args.ncolors, out.shape[0]),
        duplicates="drop",
        labels=False,
    )

    # sys.stdout.write("#"+"\t".join(keep_cols)+"\n")
    to_out = out[keep_cols]
    out_header = keep_cols
    out_header[0] = "#" + out_header[0]
    to_out.to_csv(
        args.outfile, index=False, header=out_header, sep="\t", compression="gzip"
    )

    if args.full is not None:
        out.to_csv(args.full, index=False, sep="\t", compression="gzip")
