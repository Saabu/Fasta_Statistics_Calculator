# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 10:56:45 2022

@author: ronit
"""

#import modules
import os, re, gzip, logging, argparse, io, csv, glob
from pathlib import Path
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor

#create custom formatter for argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass
parser = argparse.ArgumentParser(prog = "Fasta Quality Statistics Calculator",
                                 description= "To calculate the Total Length, GC%, N50 and L50 of given fasta files.", 
                                 usage = "%(prog)s [options]",
                                 epilog = '''\
                                     
                                     INFORMATION
                                     
1. Only SINGLE fasta file is supported. To calculate statistics for two or more files, place them in a directory. 
2. If the directory has folders with spaces in folder name, put the directory path within quotes.
3. BioPython needs to be installed before using the script.''', 
                                 formatter_class = CustomFormatter)

#add different arguments 
parser.add_argument("-f", "--logfile", default = ("logfile.log"), help = "Provide the file name or directory for the log file." "\n" "Example: -f logfile.txt")

parser.add_argument("-t", "--thread", type = int, required = True, help = "Provide the number of threads to be used for multiprocessing." "\n" "Example: -t 10")

parser.add_argument("-log", "--level", default = "warning", const = "warning", choices = ["debug", "info", "warning", "error", "critical"], nargs = "?", help = "Provide the logging level." "\n" "Example: -log debug")

parser.add_argument("-i", "--input", nargs = "?", type = Path, required = True, help = "Provide the input file or directory containing files." "\n" "Example: -i file.fasta")

parser.add_argument("-o", "--output", default = "output.csv", help = "Provide the file name or directory for the output file." "\n" "Example: -o output.csv")

args = parser.parse_args()

#create log file
logging.basicConfig(filename=args.logfile, level = args.level.upper(),
                    format="%(asctime)s - %(levelname)s:  %(message)s",
                    filemode="a+", force = True)

#check for directory or file from input
if os.path.isdir(args.input):
      files = glob.glob(os.path.join(str(args.input), "*"))
elif os.path.isfile(str(args.input)):
      files = (io.StringIO(str(args.input)))

#function for calculating the quality statistics
def calc(file):
       file_name = os.path.basename(file)   #get the file name
       
       try:
           inputfile = gzip.open(file, "rt")    #open gzipped file
           inputfile.read(1)
       except gzip.BadGzipFile:
           inputfile = open(file, "rt")     #open plain file
       seq_lengths = []; gc_list = []; temp = 0; n50 = 0; l50 = 0
       try:
           for rec in SeqIO.parse(inputfile, "fasta"):
               gc_list.append(len(re.findall("[gc]", str(rec.seq),
                                             flags=re.IGNORECASE)))     #count GC for each contig
               seq_lengths.append(len(re.findall("[atugc]", str(rec.seq),
                                                 flags=re.IGNORECASE)))     #count length for each contig
           sum_seqlen = sum(seq_lengths)
           gc_percentage = (sum(gc_list)/sum_seqlen)*100  #get GC%
           if sum(gc_list) == 0:
               logging.warning("%s has no GC content", file_name)   #log file name with no GC %
             
           gc_percentage = "{:.2f}".format(gc_percentage)
           
           seq_lengths.sort(reverse = True)
           for val in seq_lengths:      #get N50 and L50
               temp += val
               l50 += 1
               if temp >= sum_seqlen/2:
                   n50 = val
                   break

           return file_name, sum_seqlen, gc_percentage, n50, l50
       except:
           logging.error("%s is not a plain or gzipped fasta file", file_name)  #log unrecognised files
       finally:
           inputfile.close()

outlist = []
if __name__ == '__main__':
    try:
        #implement MultiProcessing
        with ProcessPoolExecutor(max_workers = args.thread) as executor:
            try:
                #create output csv file
                with open(args.output, "w", newline = '') as output:    
                    for file in files:
                        if os.path.isfile(file):
                            futures = executor.submit(calc, file)   #save output into "futures" variable
                            if futures.result() != None:
                                outlist.append(futures.result())
                            else:
                                continue
                        else:
                            continue
                    outputfile = csv.writer(output) #append the result
                    outputfile.writerow(['File Name','Total Length','GC Percentage','N50','L50'])
                    for row in outlist:
                        outputfile.writerow(row)
            except NameError:
                print(args.input, "does not exist")     #error message for wrong input file name or directory 
    except ValueError:
        print("Thread count must be higher than 0")
