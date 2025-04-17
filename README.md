### FARHA RNAseq Workflow

To run the pipeline follow the following steps:

1) Edit the `config.txt` file to set the paths to the reference genome, genome annotations in GTF format, directory for the STAR databseand full path to the directory containing the raw fastq reads.

The defaul location of the configuration file is the current working directory. You can also specify a different path and/or name using the parameter `-c`.

2) Ensure you have created a file named `samples.prefix` in your working directory containing the prefixes of the fastq read files. It is expected that the fastq files follow the format `<prefix>.R1/2.fastq.gz`.

If the samples.prefix file is in some other directory or has a different name, you can specify its location with the parameter `-p`.

3) Run the pipeline in biowulf using default values by typing:
```
./RNAseq_pipeline_v1.sh

```
If the patha to the prefix file and the config file are not the working directory, you can run the pipeline like so:
```
./RNAseq_pipeline_v1.sh -p /new/path/samples.prefix -c /my/new/config/path/config.txt
```  

The pipeline will generate the following directories:

**1-trimmed** containing the trimmed reads.

**2-alignments** with the STAR output files.

**3-deduplicated** with the deduplicated bam files (picard output).

**4-read_counts** with the read counts file.
 
