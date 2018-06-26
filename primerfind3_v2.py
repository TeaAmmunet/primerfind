# -*- coding: utf-8 -*-
"""
Created on Thu Nov 09 10:13:13 2017

@author: Tea
"""
import sys
import itertools
import timeit
#import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
import os
#import csv
#%% Functions
def get_primername(fin,txtfile,error):
    names=pd.read_csv(txtfile,sep=',')
    g={}
    for i in range(0,len(names)):
        count=[]
        fw=names.loc[i][1] 
        fwbar=names.loc[i][2]
        rvbar=names.loc[i][3]
        rv=names.loc[i][4] #3 if reading with csv, 2 if reading with pandas
        count=countprimers2(fin,fw,fwbar,rv,rvbar,error)
        group=names.loc[i][0]
        #print count
        g[group]=count[0:2]
        #s[group]=count[2]
        sdf=pd.Series(count[2],name=group)
        cwd = os.getcwd()
        #print cwd
        #save_path = 'C:/example/'
        sdf.to_csv(cwd+'/results/'+group+'_seqs.csv',sep='\t')
        #print g
    return g
def countprimers2(filein,fwprimer,fwbar,revprimer,revbar, err):
    countbothfw_rcrev=0
    countbothrcfw_rev=0
    fw=Seq(fwprimer,generic_dna)
    fwbar=Seq(fwbar,generic_dna)
    rev=Seq(revprimer,generic_dna)
    revbar=Seq(revbar,generic_dna)
    rc_fw=fw.reverse_complement()
    rc_rev=rev.reverse_complement()
    rc_fwbar=fwbar.reverse_complement()
    rc_rvbar=revbar.reverse_complement()
    f= open(sys.argv[1],'r') if len(sys.argv) > 1 else sys.stdin
    n=1
    seqs=[]
    #for line in f: #When running with awk
    for line in itertools.islice(f, 1, None, 2):
        if n % 100000 == 0:
            print ('reading in line...', n)
        fullseq=Seq(line, generic_dna)
        #Cutting the sequence into 2 pieces
        #Trimming tag but taking +3 bp extra, since I saw there were some odd things in the sequence
        seq1=fullseq[:5+len(fwbar)+len(fw)+3] # first part of seq
        seq2=fullseq[len(fullseq)-(len(fw)+len(fwbar)+5)-3:] #end part of seq
        #only for fw+rcrev and rc_fw+rev
        #finding full sequences
        fwind1=seq1.find(fw)
        revind1=seq1.find(rev)
        rcfw_ind2=seq2.find(rc_fw)
        rcrv_ind2=seq2.find(rc_rev)
        #If none of the forward direction sequences is found  in full (find() outputs -1):
        if (fwind1 == (-1)) and (rcrv_ind2 ==(-1)):
            #try finding primers with errors
            countf1, posf1=find_with_error2(seq1,fw,err)
            #if primer found with errors, set index to pos (from find_with_errors)
            if countf1 != 0:
                fwind1 = 0
            countr2, posr2=find_with_error2(seq2,rc_rev,err)
            if countr2 != 0:
                rcrv_ind2 = 0
        #if none of the reverse direction sequences are found (find() outputs -1):
        if ((rcfw_ind2 == -1) and (revind1 == -1)):
            countf2, posf2=find_with_error2(seq2,rc_fw,err)
            if countf2 != 0:
                rcfw_ind2 = 0 
            countr1, posr1=find_with_error2(seq1,rev,err)
            if countr1 != 0:
                revind1=0
        #If one of the sequences is found:
        #fw in seq1 and rcrev in seq2
        if (fwind1 != (-1)) or (rcrv_ind2 !=(-1)):
            if (fwind1 != (-1)) and (rcrv_ind2 !=(-1)): #both found in full
                indf=fwind1
                #counting positions for a new seq: from position where primers start (ind) + length of primer
                indr=rcrv_ind2+(len(rc_rev)-1)
                #getting shorter seqs: first trimming GGTAG (5) then taking seq until beginning of primer
                #or until the end of reverse primer
                short_seqf=seq1[:indf]
                short_seqr=seq2[indr+1:len(seq2)-5]
                #Trying to find barcodes from short seqs with errors
                countf, posf=find_with_error2(short_seqf,fwbar,err)
                countr, posr=find_with_error2(short_seqr,rc_rvbar,err)
                if (countf != 0) and (countr!=0):
                    countbothfw_rcrev=countbothfw_rcrev+1
                    start=fwind1+17
                    stop=(len(fullseq)-len(seq2))+rcrv_ind2
                    seqs.append(line[start:stop])
                else:
                    countbothfw_rcrev=countbothfw_rcrev
                    seqs=seqs
            elif fwind1 != (-1): #fw found as whole in seq1
                err_rc_rev, posrcrev=find_with_error2(seq2,rc_rev,err)
                if err_rc_rev != 0: #rc_rev found with errors in seq 2
                    indf=fwind1
                    indr=rcrv_ind2+(len(rc_rev)-1)
                    short_seqf=seq1[:indf]
                    short_seqr=seq2[indr+1:len(seq2)-5]
                    countf, posf=find_with_error2(short_seqf,fwbar,err)
                    countr, posr=find_with_error2(short_seqr,rc_rvbar,err)
                    if (countf != 0) and (countr!=0):
                        countbothfw_rcrev=countbothfw_rcrev+1
                        #fw1 found in full, fw1 primer start positions gotten as is from find() -> seq start = primer start+17
                        #rev primer position is an end position from seq2->first position =end-17 ->len(full)-len(seq2)
                        start=fwind1+17
                        stop=(len(fullseq)-len(seq2))+(posrcrev-17)
                        seqs.append(line[start:stop])
                    else:
                        countbothfw_rcrev=countbothfw_rcrev
                        seqs=seqs
                else:
                    countbothfw_rcrev=countbothfw_rcrev
                    seqs=seqs
            elif rcrv_ind2 != (-1): #rc_rev found in seq2 in whole
                err_fw, posfw=find_with_error2(seq1,fw,err)
                if err_fw != 0: #fw found with errors from seq1!
                    indf=fwind1
                    indr=rcrv_ind2+(len(rc_rev)-1)
                    short_seqf=seq1[:indf]
                    short_seqr=seq2[indr+1:len(seq2)-5]
                    countf, posf=find_with_error2(short_seqf,fwbar,err)
                    countr, posr=find_with_error2(short_seqr,rc_rvbar,err)
                    if (countf != 0) and (countr!=0):
                        countbothfw_rcrev=countbothfw_rcrev+1
                        #reverse primer found in full -> start position of rev primer in seq2
                        #forward primer with errors -> stop position of fwprimer = start pos of seq 
                        start=fwind1
                        stop=(len(fullseq)-len(seq2))+rcrv_ind2
                        seqs.append(line[start:stop])
                    else:
                        countbothfw_rcrev=countbothfw_rcrev
                        seqs=seqs
                else:
                    countbothfw_rcrev=countbothfw_rcrev
                    seqs=seqs
        #rcfw in seq 2 or rev in seq1
        elif ((rcfw_ind2 != -1) or (revind1 != -1)): 
            if (rcfw_ind2 != -1) and (revind1 != -1):#both found in full
                indf=rcfw_ind2+(len(rc_fw)-1)
                indr=revind1
                short_seqf=seq2[indf+1:len(seq2)-5]
                short_seqr=seq1[:indr]
                countf, posf=find_with_error2(short_seqf,rc_fwbar,err)
                countr, posr=find_with_error2(short_seqr,revbar,err)
                if (countf != 0) and (countr!=0):
                    countbothrcfw_rev=countbothrcfw_rev+1
                    #both in full -> start positions from find
                    start=revind1+17
                    stop=(len(fullseq)-len(seq2))+rcfw_ind2
                    seqs.append(line[start:stop])
                else:
                    countbothrcfw_rev=countbothrcfw_rev
                    seqs=seqs
            elif (rcfw_ind2 != -1): #only rcfw found from seq2 in full
                err_rev, posrev=find_with_error2(seq1,rev,err)
                if err_rev != 0:
                    indf=rcfw_ind2+(len(rc_fw)-1)
                    indr=revind1
                    short_seqf=seq2[indf+1:len(seq2)-5]
                    short_seqr=seq1[5:indr]
                    countf, posf=find_with_error2(short_seqf,rc_fwbar,err)
                    countr, posr=find_with_error2(short_seqr,revbar,err)
                    if (countf != 0) and (countr!=0):
                        countbothrcfw_rev=countbothrcfw_rev+1
                        #rc_fw found in full from seq2
                        #revind found with errors -> end position
                        start=revind1
                        stop=(len(fullseq)-len(seq2))+rcrv_ind2
                        seqs.append(line)
                    else:
                        countbothrcfw_rev=countbothrcfw_rev
                        seqs=seqs
                else:
                    countbothrcfw_rev=countbothrcfw_rev
                    seqs=seqs
            elif (revind1 != -1): #only rev found from seq1 in full
                err_rc_fw, posrcfw=find_with_error2(seq2,rc_fw,err)
                if err_rc_fw !=0:
                    indf=rcfw_ind2+(len(rc_fw)-1)
                    indr=revind1
                    short_seqf=seq2[indf+1:len(seq2)-5]
                    short_seqr=seq1[5:indr]
                    countf, posf=find_with_error2(short_seqf,rc_fwbar,err)
                    countr, posr=find_with_error2(short_seqr,revbar,err)
                    if (countf != 0) and (countr!=0):
                        countbothrcfw_rev=countbothrcfw_rev+1
                        #rev found in full ->start position from find
                        #rcfw found with errors, stop position from find_with_errs -> stop pos - 17 = start pos in seq2
                        start=revind1+17
                        stop=(len(fullseq)-len(seq2))+(posrcfw-17)
                        seqs.append(line[start:stop])
                    else:
                        countbothrcfw_rev=countbothrcfw_rev
                        seqs=seqs
                else:
                    countbothrcfw_rev=countbothrcfw_rev
                    seqs=seqs
        else:
            countbothfw_rcrev=countbothfw_rcrev
            countbothrcfw_rev=countbothrcfw_rev
            seqs=seqs
        n=n+2
       
    return countbothfw_rcrev, countbothrcfw_rev, seqs

def find_with_error2(seq,primer,error):
    """find primer in seq allowing for error-number of mismatches.
    seq and primer should be given as strings"""
    count=0
    dist=0
    lenp=len(primer)
    comp_seq=list(str(seq))
    pos=0
    for i in range(0,len(seq)):
        start=i
        end=start+lenp
        if end <=len(seq):
            lseq=list(str(seq))
            comp_seq=lseq[start:end+1]
            lst=zip(comp_seq, list(str(primer)))
            #Distance by breaking if more than error nr of differences
            dist=0
            pos=i
            for k in lst:
                d=(lambda (x, y) : 0 if x == y else 1)(k)
                dist=dist+d
                pos=pos+1
                if dist > error:
                    count=count
                    pos=pos
                    break
            if dist <=error:
                count=count+1#count ONE match
                pos=pos
        else:
            break
    return count, pos
#%%
def main():
    start_time = timeit.default_timer()
    filein=sys.stdin
    txtf='Sample_primers_barcodes.csv'
    prims=get_primername(filein,txtf,2)
    #print(prims)
    res_df=pd.DataFrame(prims.items(),columns={'group','number'})#,'seqlines'})
    #print prims.viewitems()
    #print k,v in prims.iteritems() if not(v==(0,0,0,0))
    print (res_df[res_df['number']!=(0,0)])
    res_df.to_csv('primerfind3_results.csv', sep='\t')
    elapsed = timeit.default_timer() - start_time
    print('Calculation time %1.5f' % elapsed)
    #test=countprimers(filein,'TCCGTAGGTGAACCTGCGCGCACGCACTACAGA','',2)
    #print test
    
    return

main()

            
            
        
