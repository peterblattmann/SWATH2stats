#!/usr/bin/env python

from __future__ import print_function
from mimetypes import guess_type
import sys
import csv
import gzip
import re
import gflags


gflags.DEFINE_string("input_file", None, "input file name (must be csv).")
gflags.DEFINE_string("output_file", None, "output file name")
FLAGS = gflags.FLAGS


def check_input_file_format(input_file):
    """Checks csv type and compression state of input file.

    Returns a tuple with the dialect to use in a csv reader object and
    an appropriate function to open the file.
    """
    csv_type, compression_type = guess_type(input_file)

    if compression_type == 'gzip':
        smart_open = gzip.open
    else:
        smart_open = open

    if csv_type == 'text/csv':
        with smart_open(input_file, 'rt') as fin:
            headers = fin.readline()
            if re.search(',', headers):
                file_dialect = csv.excel
            elif re.search('\t', headers):
                file_dialect = csv.excel_tab
            else:
                print('Unsupported input file type!')
                print('Supported file types are tab and comma separated csv files.')
                exit(1)
    else:
        print('Unsupported input file type!')
        print('Supported file types are tab and comma separated csv files.')
        exit(1)

    return (file_dialect, smart_open)


def peptide_to_transition(input_file, output_file, file_dialect, smart_open):
    """This function is responsible for the actual conversion from one
    peptide per row to one transition per tow
    """
    with smart_open(input_file, 'rt') as f, open(output_file, 'w') as out:
        csvf = csv.DictReader(f, dialect=file_dialect)
        out.write(
            "ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,"
            "IsotopeLabelType,Intensity,RT,Condition,BioReplicate,Run\n")

        for line in csvf:
            name = line["align_origfilename"].split("/")[-1]
            fas = line["aggr_Fragment_Annotation"].split(";")
            pas = line["aggr_Peak_Area"].split(";")
            for fa, pa in zip(fas, pas):
                out.write("%s,%s,%s,%s,NA,light,%s,%s,%s,%s,%s\n" %
                          (line["ProteinName"], line["FullPeptideName"],
                           line["Charge"], fa, pa, line.get("RT", ""),
                           line.get("Condition", ""), line.get("BioReplicate", ""),
                           name))


def main():
    "Converts from one row per peptide to one row per transition."
    (file_dialect, smart_open) = check_input_file_format(FLAGS.input_file)
    peptide_to_transition(FLAGS.input_file, FLAGS.output_file, file_dialect, smart_open)


if __name__ == "__main__":
    gflags.MarkFlagsAsRequired(["input_file", "output_file"])
    argv = FLAGS(sys.argv)
    main()
