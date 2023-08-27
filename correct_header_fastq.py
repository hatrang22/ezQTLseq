######################################################################
## This file is dedicated to correct the header format of fastq2 files 
## when running remove duplicate rule with SPET mode
## The default required format is ... 2:N:0
######################################################################

import argparse
import os
import gzip

def argparser():
    # Create the argument parser
    parser = argparse.ArgumentParser()
    
    # Add the required argument
    parser.add_argument("-fastq_files", nargs='+', help="List of fastq files")

    # Parse the command-line arguments
    args = parser.parse_args()
    
    return args

def correct_header(args):

    for fpath in args.fastq_files:
        with gzip.open(fpath, "rt") as f:
            header = f.readline()
            if not header.split(" ")[-1].startswith("2:N:0"):
                print(f"Correcting header of {os.path.basename(fpath)}...")
                fcontent = f.readlines()
                fcontent.insert(0, header)
                for i, line in enumerate(fcontent):
                    if line.startswith("@"):
                        fcontent[i] = line.split("/")[0] + " 2:N:0\n"
                with gzip.open(fpath, "wt") as f:   
                    f.writelines(fcontent)
                print(f"Finished correcting header of {os.path.basename(fpath)}.")
    

if __name__ == "__main__":
    
    correct_header(argparser())