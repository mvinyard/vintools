import time

def get_sample_names(path):

    import os
    import glob

    fastqs = sorted(glob.glob(path + "*"))

    samples = [os.path.basename(file).split(".")[0].split("_")[0] for file in fastqs]
    samples = list(dict.fromkeys(samples))

    return samples

def assign_LR_reads(sample, base_path):

    read_one = base_path + sample + "_1.fastq"
    read_two = base_path + sample + "_2.fastq"

    return read_one, read_two

class STAR:
    def generate_genome(ref_genome, threads=10):

        import subprocess
        
        ref_genome_dir = ref_genome
        ref_genome_file = ref_genome + ".fa"

        command = (
            "STAR --runMode genomeGenerate --runThreadN "
            + str(threads)
            + " --genomeDir "
            + ref_genome_dir
            + " --genomeFastaFiles "
            + ref_genome
        )
        print("Running: ", command)

        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

        return output, error

    def align_reads(
        ref_genome,
        read_one,
        read_two,
        outfile_name,
        threads=10,
        outfile_format="BAM",
        zcat="zcat",
        outFilterMultimapNmax=2,
        outFilterMismatchNmax=3,
        alignIntronMax=1,
        alignEndsType="EndToEnd",
    ):
        import os
        import subprocess
        
        ref_genome_dir = ref_genome
        ref_genome_file = ref_genome + ".fa"
        
        if os.path.exists(ref_genome_dir + "genomeParameters.txt") !=True:
            generate_genome(ref_genome_dir, ref_genome_file, threads)
            
        command = (
            "STAR --runThreadN "
            + str(threads)
            + " --genomeDir "
            + ref_genome_dir
            + " --readFilesIn "
            + read_one
            + " "
            + read_two
            + " --outFileNamePrefix "
            + outfile_name
            + " --outSAMtype "
            + outfile_format
            + " SortedByCoordinate --outFilterMultimapNmax "
            + str(outFilterMultimapNmax)
            + " --outFilterMismatchNmax "
            + str(outFilterMismatchNmax)
            + " --alignIntronMax "
            + str(alignIntronMax)
            + " --alignEndsType "
            + alignEndsType
        )

        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

        return output, error

import time
from tqdm.notebook import tqdm as progress_bar

def _STAR_align(fastqs, ref_genome, alignment_outs, threads=10, viz_metrics=False):
    
    samples = get_sample_names(fastqs)
    
    alignment_out_messages = []
    alignment_err_messages = []

    for sample in progress_bar(samples, desc="Alignment progress"):
        
        sample_start_time = time.time()

        print("Now processing sample: ", sample)
        left_read, right_read = assign_LR_reads(sample, fastqs)
        out, err = STAR.align_reads(ref_genome, left_read, right_read,(alignment_outs + sample))
        alignment_out_messages.append(out)
        alignment_err_messages.append(err)
        
        sample_end_time = time.time()
        
        sample_time = (sample_end_time - sample_start_time) / 3600
        print("Sample ", sample, "was aligned in ", sample_time, " hours.")
        print("Aligned read files (.bam) are stored at: ", alignment_outs + sample)
        print("Metrics: ")
        
        if viz_metrics == True:
            metrics_file = alignment_outs + sample + "Log.final.out"

            metrics = open(metrics_file, "r")
            print(metrics.read())

        
    return alignment_out_messages, alignment_err_messages