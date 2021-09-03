def _deinterleaf_fastq(
    interleafed_file,
    outfile_prefix="./",
    silent=False,
    interleaf_script_path="~/scripts/deinterleaf.sh",
):

    """
    Parameters:
    -----------
    interleafed_file
        type: str
        
    outfile_prefix [ optional ]
        default: './'
        type: str
        
    silent [ optional ]
        default: False
        type: bool
        
    interleaf_script_path [ optional ]
        default: "~/scripts/deinterleaf.sh"
        type: str
    
    Returns:
    --------
    None
    
    deinterleaf.sh script contents:
    -------------------------------
    
    #!/bin/bash

    # deinterleaf.sh

    interleafed_file=$1
    outfile_prefix=$2

    paste - - - - - - - - < ${interleafed_file} \
        | tee >(cut -f 1-4 | tr "\t" "\n" > ${outfile_prefix}_read1.fastq) \
        |       cut -f 5-8 | tr "\t" "\n" > ${outfile_prefix}_read2.fastq

    """

    for file, out_prefix in zip(interleafed_file, outfile_prefix):
        command = " ".join(["bash", interleaf_script_path, file, out_prefix])
        os.system(command)