#!/bin/python3
# Tables of the results of run_part_tese.sh (the COST files of out_folder).
#
# Usage: python3 table_part_tese.py [out_folder] [csv_file]
#   out_folder  out_folder of run_part_tese.sh (default ../TESE_OUT/ANS_TESE_PART)
#   csv_file    if given, the data of all instances is also written to it
#
# For each size and operation, and for each parameter (op, I, L and F), the
# table has, per value of the parameter (averaging over the other parameters),
# the mean cost (number of breakpoints of the two partitions) and the mean time
# (seconds) of each algorithm. With two or more algorithms, it also has, for
# each algorithm other than the first, the mean cost relative to the first one
# and the percentage of instances in which it is better/worse than the first.
# Instances without a partition (lines NA) are counted in the column NA and
# left out of the other columns.

import os
import re
import sys

import pandas as pd

out_folder = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "TESE_OUT", "ANS_TESE_PART")
csv_file = sys.argv[2] if len(sys.argv) > 2 else None

cost_folder = os.path.join(out_folder, "COST")
params = ["op", "I", "L", "F"]
name_re = re.compile(r"^N(\d+)_(.+?)(\d+)_I(\d+)_L(\d+)_F(\d+)\.txt$")

rows = []
for algo in sorted(os.listdir(cost_folder)):
    for root, _, files in os.walk(os.path.join(cost_folder, algo)):
        for file in files:
            m = name_re.match(file)
            if not m:
                continue
            n, mod, op, i, l, f = m.groups()
            with open(os.path.join(root, file)) as fin:
                for inst, line in enumerate(fin):
                    fields = line.split()
                    na = not fields or fields[0] == "NA"
                    rows.append({
                        "algo": algo, "N": int(n), "mod": mod, "op": int(op),
                        "I": int(i), "L": int(l), "F": int(f), "inst": inst,
                        "cost": None if na else int(fields[2]),
                        "time": None if na else float(fields[7])})
data = pd.DataFrame(rows)
if csv_file:
    data.to_csv(csv_file, index=False)

# the first algorithm is the reference (Greedy if it is there)
algos = sorted(data["algo"].unique(), key=lambda a: (a != "Greedy", a))
ref = algos[0]

def table(group):
    """Table of group (the instances of a size and operation) by each param."""
    key = ["N", "mod", "op", "I", "L", "F", "inst"]
    wide = group.pivot_table(index=key, columns="algo",
                             values=["cost", "time"], dropna=False)
    na = group.assign(na=group["cost"].isna()).pivot_table(
        index=key, columns="algo", values="na", aggfunc="sum")
    for param in params:
        by = wide.groupby(level=param)
        out = pd.DataFrame(index=by.size().index)
        for algo in algos:
            out["cost " + algo] = wide[("cost", algo)].groupby(level=param).mean()
        for algo in algos[1:]:
            cost, cost_ref = wide[("cost", algo)], wide[("cost", ref)]
            both = cost.notna() & cost_ref.notna()
            rel = (cost / cost_ref)[both].groupby(level=param).mean()
            better = (cost < cost_ref)[both].groupby(level=param).mean()
            worse = (cost > cost_ref)[both].groupby(level=param).mean()
            out[algo + "/" + ref] = rel
            out["%better " + algo] = 100 * better
            out["%worse " + algo] = 100 * worse
        for algo in algos:
            out["time " + algo] = wide[("time", algo)].groupby(level=param).mean()
        n_na = na.groupby(level=param).sum()
        if n_na.to_numpy().sum() > 0:
            for algo in algos:
                out["NA " + algo] = n_na[algo]
        out.loc["all"] = [
            out[c].mean() if not c.startswith("NA ") else out[c].sum()
            for c in out.columns]
        print(out.round(3).to_string())
        print()

pd.set_option("display.width", 250)
for (n, mod), group in data.groupby(["N", "mod"]):
    print("=" * 30, "N{} {}".format(n, mod), "=" * 30)
    table(group)
