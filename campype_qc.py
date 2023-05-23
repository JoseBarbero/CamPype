import os
import sys
from terminal_banner import Banner
from subprocess import call
from campype import read_input_files

def fastqc_call(input_file_fw, input_file_rv, output_dir):
    """
    FastQC call.
    
    Arguments:
        input_file_fw {string} -- Forward input file (and route).
        input_file_rv {string} -- Reverse input file (and route).
        output_dir {string} -- Output directory.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["fastqc", input_file_fw, input_file_rv, "-o", output_dir, "-q"]
    return call(arguments)

def multiqc_call(fastqc_dir, output_dir):
    """
    MultiQC call.
    
    Arguments:
        fastqc_dir {string} -- FastQC output directory.
        output_dir {string} -- MultiQC output directory.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["multiqc", fastqc_dir, "-o", output_dir]
    return call(arguments)



if __name__ == "__main__":
    
    output_folder = sys.argv[1]

    fastqc_dir = output_folder+"/fastq_quality_control"
    multiqc_dir = fastqc_dir+"/multiqc"

    # Create output directories
    os.mkdir(fastqc_dir)
    os.mkdir(multiqc_dir)
    
    # Run quality control (FastQC and MultiQC)  
    print(Banner(f"\nFastq reads quality control analysis: FastQC\n"), flush=True)

    for sample_basename, data in read_input_files("input_files.csv"):
        
        sample_fw = data["FW"]
        sample_rv = data["RV"]    
        
        # Fastqc needs the output directory to exist
        os.mkdir(fastqc_dir+"/"+sample_basename)

        # FastQC call
        fastqc_call(sample_fw, sample_rv, fastqc_dir+"/"+sample_basename)
        print(f"Analysis complete for {sample_fw}", flush=True)
        print(f"Analysis complete for {sample_rv}", flush=True)

    # MultiQC call
    print(Banner(f"\nMerge FastQC reports into a single report: MultiQC\n"), flush=True)
    multiqc_call(fastqc_dir, multiqc_dir)
    print(Banner(f"\nDONE\n"), flush=True)