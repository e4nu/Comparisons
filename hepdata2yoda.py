#!/usr/bin/env python3
"""Reference data of the Rivet analysis from the HEPData submission.

Reads the differential cross-section tables of the HEPData submission (YAML) and writes them as
YODA Estimate1D objects, named as HEPData names them in its own YODA export: table N of
submission.yaml becomes d<N>-x01-y01. Every uncertainty source is kept separately (stat, sys,...),
so Rivet shows their sum in quadrature. Only the cross-section tables are converted; the
acceptance and radiative correction tables are not something a generator can be compared to.

Once the record is public, the YODA file from HEPData (rivet-mkhepdata / "Download YODA")
gives the same objects and can replace this one.

Usage:
  python3 hepdata2yoda.py <HEPData submission folder> CLAS_2026_I0000000.yoda.gz
Needs the yoda Python module (from the YODA installation) and PyYAML.
"""
import os
import sys

import yaml
import yoda

ANALYSIS = "CLAS_2026_I0000000"


def tables(submission_dir):
    """(table number, data file) of every table, in the order of submission.yaml (1-based)."""
    with open(os.path.join(submission_dir, "submission.yaml")) as f:
        docs = [d for d in yaml.safe_load_all(f) if d and "data_file" in d]
    return [(i + 1, d["data_file"], d.get("name", "")) for i, d in enumerate(docs)]


def to_estimate(path, table, title):
    x = table["independent_variables"][0]
    y = table["dependent_variables"][0]
    edges = [b["low"] for b in x["values"]] + [x["values"][-1]["high"]]
    est = yoda.Estimate1D(edges, path)
    est.setTitle(title)
    est.setAnnotation("IsRef", 1)
    for i, v in enumerate(y["values"], start=1):
        b = est.bin(i)
        b.setVal(v["value"])
        for e in v.get("errors", []):
            b.setErr([-e["symerror"], e["symerror"]], e["label"])
    return est


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    submission_dir, out_file = sys.argv[1], sys.argv[2]
    out = {}
    for n, data_file, name in tables(submission_dir):
        if not data_file.startswith("differential_cross_section"):
            continue
        with open(os.path.join(submission_dir, data_file)) as f:
            table = yaml.safe_load(f)
        path = "/REF/%s/d%02d-x01-y01" % (ANALYSIS, n)
        out[path] = to_estimate(path, table, name)
        print("%s  <-  %s (%d bins)" % (path, data_file, out[path].numBins()))
    yoda.write(list(out.values()), out_file)
    print("Wrote", out_file)


if __name__ == "__main__":
    main()
