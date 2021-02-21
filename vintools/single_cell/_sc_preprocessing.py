import pandas as pd
import time
import pysam
import os

def get_CB_UB(read):
    
    """Input is a pysam read generated using samfile.fetch"""
    
    cb_ub = read.qname[-17:].split("_")
    CB, UB = cb_ub[0], cb_ub[1]
    
    return CB, UB

def _add_barcode_tags(bcs, untagged_bam, out):
    
    """
    Add a cell barcode and UMI barcode as a tag from the read name. 
    Format of tags: CB:Z:[8 nt cell barcode] UB:Z:[8 nt UMI barcode]
    
    Parameters:
    -----------
    bcs
        array of barcodes. helps if it is already filtered. 

    untagged_bam
        unbarcoded bam
        
    out
        destination of the barcoded bam
    
    """
    
    untagged_bam = pysam.AlignmentFile(untagged_bam, "rb")
    tagged_bam = pysam.AlignmentFile(out, "wb", template=untagged_bam)
    
    barcoded_read_count = 0
    start_time = time.time()
    for i, read in enumerate(untagged_bam.fetch(until_eof=True)):
        
        CB_, UB_ = get_CB_UB(read)
        
        if CB_ in bcs:
            tags = read.get_tags()
            tags.append(("CB", CB_, "Z"))
            tags.append(("UB", UB_, "Z"))
            read.set_tags(tags)
            tagged_bam.write(read)
            
            barcoded_read_count += 1
            
        if i % 10000000 == 0:
            if i !=0:
                current_time = time.time()
                elapsed_time = (current_time - start_time) / 60
                percent_barcoded = (barcoded_read_count / i)*100
                print("Number of reads processed:", i)
                print("{:.2f}".format(elapsed_time), "minutes elapsed.")
                print("Number of reads barcoded:", barcoded_read_count, "Percentage barcoded:", "{:.2f}".format(percent_barcoded), "%")
            
            
    untagged_bam.close()
    tagged_bam.close()

    
def _split_bam(barcodes, unsplit_bam, outdir):
    
    """
    Function for splitting a bulk bam into sc-bams based on cell barcode ("CB" tag).
    
    Parameters:
    -----------
    barcodes
        single-column of a pandas df containing barcodes. order shouldn't matter.
    
    
    
    Notes:
    ------
    1. Be sure to first sort bam file by cell barcode 'CB'
        ex: samtools sort -t CB -@ 20 possorted_bam.dedup.bam > possorted_bam.dedup.sorted.bam
    
    Adapted from the original author, Dr. Huidong Chen
    """
    
    samfile = pysam.AlignmentFile(unsplit_bam, "rb")
        
    if(not os.path.exists(outdir)):
        os.makedirs(outdir)

    # variable to hold barcode index
    CB_hold = 'unset'
    itr = 0
    # read in upsplit file and process reads line by line

    start_time = time.time()
    for read in samfile.fetch(until_eof=True):
        if(len([ x for x in read.tags if x[0] == "CB"])>0):
            # barcode itr for current read
            CB_itr = read.get_tag( 'CB')
            if(CB_itr in barcodes.tolist()): # barcodes['cell']
                # if change in barcode or first line; open new file  
                if(CB_itr!=CB_hold or itr==0):
                    # close previous split file, only if not first read in file
                    if(itr!=0):
                        print(CB_itr)
                        split_file.close()
                    CB_hold = CB_itr
                    itr+=1
                    split_file = pysam.AlignmentFile(outdir + "{}.bam".format(CB_hold), "wb", template=samfile)
                # write reads with the same barcode to file
                split_file.write(read)    
    
    split_file.close()
    samfile.close()
    end_time = time.time()


    print(end_time - start_time)