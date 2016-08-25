'''
Created on Nov 14, 2014

@author: Igor Ruiz de los Mozos

Script will convert xlink file with counts in column 4 and creates a single line for each xlink cDNA. If using B file that have a decimal value in count (beacuse is a multimaper read) it assingns a 1 count.

Usage:
        python BED2BED_no_counts.py input_file.bed output_file.bed


'''


import sys

def BED2BED_no_counts(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
	count = float(col[3])

	if (count) >= 1:
        	
		
		for i in range(0,int(count)):
        	    fout.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + "" + '\t' + col[4] + "\t" + col[5] + '\n')
	
	elif (count) < 1:

        	#count = 1
		fout.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + "" + '\t' + col[4] + "\t" + col[5] + '\n')
        
	else:
		print col[3]

	line = fin.readline()
                    
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    BED2BED_no_counts(fname_in, fname_out)
else:
    print("python BED2BED_no_counts.py <input_file> <output_file>")
