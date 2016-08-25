'''
Created on 12 Feb 2016

@author: Igor Ruiz de los Mozos

Lift over Alu elements on Gorilla Resus 

'''
import sys
import os
import subprocess
from subprocess import call


def listdir_fullpath(d):
   return [os.path.join(d, f) for f in os.listdir(d)]    #example .... faFileLists = listdir_fullpath(outputdir + "/indexed")

def runUNIXCommands(commands):
    print commands
    try:
        retcode = subprocess.check_call(commands, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            print >>sys.stderr, "Child returned", retcode
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e



## Not needed function
def liftOver(bed_input, chain_file, bed_output, unMapped):
    cmd_list = ['/home/igor/Programs/bin/liftOver', bed_input, chain_file, bed_output, unMapped]
    p = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sys.stdout.flush()
    p.wait()
    err = p.stderr.read().strip()
    out = p.stdout.read().strip()
    if err:
            print err + "    ERROR"
    if out:
            print out + "    StandartOUT"
    sys.stdout.flush()
    print cmd_list
    

def read_in_line_save_to_file(fname_in, directory_name_out):    
    fin = open(fname_in, "rt")
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        file_name_out = directory_name_out + col[0]+"_"+ col[1]+"_"+col[2]+".bed"
        fout = open(file_name_out, "w")
        fout.write(col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t'  + col[4] + '\t' + col[5] +'\n')
        line = fin.readline()
        
    return()

read_in_line_save_to_file(fname_in, directory_name_out) 

outputdir = '/media/igor/DATA/UCL/Evolution_Alus/liftover_aproach/temp/'

#faFileLists = listdir_fullpath(outputdir + "/indexed")
bedFileLists = listdir_fullpath(outputdir)
bedFileLists = [each for each in bedFileLists if each.endswith('.bed')]
print bedFileLists ## list of files

runUNIXCommands("mkdir -p " + outputdir + "/Gor")
runUNIXCommands("mkdir -p " + outputdir + "/Join/Align")
runUNIXCommands("mkdir -p " + outputdir + "/Gor/fastas")
runUNIXCommands("mkdir -p " + outputdir + "/Human/fastas")

for bed in bedFileLists:
    tempOutName = outputdir + "Gor/" + bed.split("/")[-1].replace(".bed", "_Gor.bed")
    tempOutName_unmaped = outputdir + "Gor/" + bed.split("/")[-1].replace(".bed", "_Gor_unmaped.bed")
    tempOutName_join = outputdir + "Join/" + bed.split("/")[-1].replace(".bed", "_Join.fasta")
    print bed
    print tempOutName
    liftoverCommand = "/home/igor/Programs/bin/liftOver " + bed + " " + gorila_chain + " " + tempOutName + " " + tempOutName_unmaped
    runUNIXCommands(liftoverCommand)
    runUNIXCommands("rm -r " + tempOutName_unmaped)   ## Files that don't have an homologue in this organism
    filein = open(tempOutName, "r+")
    line = filein.readline()

    if os.path.getsize(tempOutName) == 0:   ## Remove empty files and continue (dont get fasta)
        runUNIXCommands("rm -r " + tempOutName) 
        continue
    else:
        ## get fasta en gorila
        tempOutName_gor_fasta = outputdir + "Gor/fastas/" + bed.split("/")[-1].replace(".bed", "_Gor.fasta")
        runUNIXCommands("bedtools getfasta -s -fi '/media/igor/DATA/UCL/Genomes/Primates/gorGor3.1/gorGor3.fa' -bed " + tempOutName + " -fo " + tempOutName_gor_fasta)
        ## get fasta en humano
        tempOutName_human_fasta = outputdir + "Human/fastas/" + bed.split("/")[-1].replace(".bed", "_Human.fasta")
        runUNIXCommands("bedtools getfasta -s -fi '/media/igor/DATA/UCL/Genomes/Human/ucsc.hg19.fasta' -bed " + bed + " -fo " + tempOutName_human_fasta)

        #runUNIXCommands("cat " + aluseq + " " + tempOutName_human_fasta + " " + tempOutName_gor_fasta + " > " + tempOutName_join)   ### Join with 
        runUNIXCommands("cat " + tempOutName_human_fasta + " " + tempOutName_gor_fasta + " > " + tempOutName_join)
        
        tempOutName_align = outputdir + "Join/Align/" + bed.split("/")[-1].replace(".bed", "_Align.fasta")
        runUNIXCommands("/home/igor/Programs/bin/muscle3.8.31_i86linux64 -in " + tempOutName_join + " -out " + tempOutName_align)
        
    filein.close()
    

    #runUNIXCommands("rm -r " + tempOutName_unmaped)

#def getfasta():
    
    


'''

aluseq = '/media/igor/DATA/UCL/Evolution_Alus/liftover_aproach/Alu_seq_model.fasta'
fname_in = '/media/igor/DATA/UCL/Evolution_Alus/liftover_aproach/alternative_Alu_exons_test.bed'
directory_name_out = '/media/igor/DATA/UCL/Evolution_Alus/liftover_aproach/temp/'
out_elememt= '/media/igor/DATA/UCL/Genomes/Human/Alus/Alu_element_name.txt'
gorila_chain= '/media/igor/DATA/UCL/Evolution_Alus/liftover_aproach/liftover_chains/hg19ToGorGor3.over.chain'

'''