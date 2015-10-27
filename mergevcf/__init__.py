"""
Entry point for mergevcf, a program to merge VCF files from multiple
callers into a single VCF, paying careful attention to normalise
structural variant calls.

https://github.com/ljdursi/mergevcf
"""
import mergedfile
import argparse
import os
import sys

def main():
    """Merge VCF files, output to stdout or file"""
    defsvwindow = 100

    parser = argparse.ArgumentParser(description='Merge calls in VCF files')
    parser.add_argument('input_files', nargs='+', help='Input VCF files')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="Specify output file (default:stdout)")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Specify verbose output")
    parser.add_argument('-l', '--labels', nargs='+',
                        help='Optional labels for each input VCF file (default:basenames)')
    parser.add_argument('-s', '--sv', action='store_true',
                        help='Force interpretation as SV (default:false)')
    parser.add_argument('-f', '--filtered', action='store_true',
                        help='Include records that have failed one or more filters (default:false)')
    parser.add_argument('-w', '--svwindow', default=defsvwindow, type=int,
                        help='Window for comparing breakpoint positions for SVs (default:'+str(defsvwindow)+')')

    args = parser.parse_args()
    input_files = args.input_files
    if args.labels is None:
        labels = [os.path.splitext(os.path.basename(f))[0] for f in input_files]
    else:
        labels = args.labels

    mergedfile.merge(input_files, labels, args.sv, args.output, 
                     slop=args.svwindow, verbose=args.verbose, 
                     filter_chrom=True, no_filter_variants=args.filtered)
