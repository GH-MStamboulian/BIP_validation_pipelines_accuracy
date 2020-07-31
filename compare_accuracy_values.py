#!/cm/shared/apps/python/3.6.5/bin/python
import argparse
import shutil
import os
import subprocess

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

class CustomFormatter(argparse.RawDescriptionHelpFormatter):
    pass

epilog ="""
    Steps performed by this script:
    1. check to see if BIP outputs exist, or else schedule bip runs
    2. Run LoD scripts over the different samples and different bip versions
    3. create an LoD summary table for the variants of interest in the study
    4. create a reportable table comparing the LoD values between the two versions of bip
    """

parser = argparse.ArgumentParser(description="BIP accuracy comparison runner",epilog=epilog, formatter_class=CustomFormatter)
parser.add_argument("-o", "--output_dir",
                    action="store",
                    dest="output_dir",
                    required=True,
                    help=("directory to store the tables"))
parser.add_argument("-d1", "--_data_tables_dir1",
                    action="store",
                    dest="_data_tables_dir1",
                    required=True,
                    help=("directory of the sample and variant tables for the study"))
parser.add_argument("-v1", "--version1",
                    action="store",
                    dest="version1",
                    required=True,
                    help=("specify the version of the BIP used"))
parser.add_argument("-d2", "--_data_tables_dir2",
                    action="store",
                    dest="_data_tables_dir2",
                    required=True,
                    help=("directory of the sample and variant tables for the study"))
parser.add_argument("-v2", "--version2",
                    action="store",
                    dest="version2",
                    required=True,
                    help=("specify the version of the BIP used"))

args = parser.parse_args()

if not os.path.isdir(args._data_tables_dir1) and not os.path.isdir(args._data_tables_dir2):
    raise RuntimeError("Nonexistant or invalid input directory '{}'".format(args.in_dir))

abs_dir = subprocess.Popen(["pwd"], stdout=subprocess.PIPE, shell=True).communicate()[0].decode('utf-8').replace('\n','/')

data_tables_dir=args._data_tables_dir1
out_dir=args.output_dir
version1=args.version1

out_dir1 = out_dir+'/out_1_bip_'+version1+'/'
if not os.path.isdir(out_dir1):
    os.makedirs(out_dir1)

os.system("python3 analysis.py " + data_tables_dir + "line_data.xlsx " + out_dir1)

data_tables_dir=args._data_tables_dir2
out_dir=args.output_dir
version2=args.version2

out_dir2 = out_dir+'/out_2_bip_'+version2+'/'
if not os.path.isdir(out_dir2):
    os.makedirs(out_dir2)

os.system("python3 analysis.py " + data_tables_dir + "line_data.xlsx " + out_dir2)

os.system("python3 combineAccuracyTables.py " + out_dir1 + " " + version1 + " " +out_dir2 + " " + version2 + " " + out_dir)
