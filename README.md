# Transcriptomics Workflow

Instructions to take raw transcriptomics reads and progress to a counts file that can be used in `R` packages for statistical analysis of differential expression.

To follow this tutorial you will need to have some basic familiarity with working in a unix-like environment.  An excellent place to start would be to attend a [software carpentry](http://software-carpentry.org/) or [data carpentry](http://www.datacarpentry.org/) workshop.  The VLSCI also have several tutorials for [working with unix and/or a HPC environment](http://vlsci.github.io/lscc_docs/tutorials/)

## Connect to analysis server
This will usually mean opening a terminal and using `ssh` to connect. The general format of this command is;

```bash
	ssh username@my.hpc.server.edu.au
```

## Download your files
Often the sequencing center will provide a url to download raw sequencing data do that a tool like wget can be used to do the download. For example;

```bash
	wget -O SequencingData.tar https://data.sequencingcenter.edu.au/download?key
```

Once you have downloaded the data it is a good idea to keep the original unmodified files somewhere safe. 

## Unpack the files

How to do this depends on how the files were packaged up. The general tool to use is `tar`, possibly with options to unzip.  Assuming the downloaded data is in a file called `Project_ANDI1971.tar` we would use this command to unpack.

```bash
	tar -xvf Project_ANDI1971.tar
```

Alternatively for gzipped files (`.gz` extension) you can use the `-z` option to unzip and untar at the same time. 

```bash
	tar -zxvf Project_ANDI1971.tar.gz
```


## Organise Files

Now organise your raw read files. For the purposes of this tutorial we will assume that all your raw reads specific to this project have been placed in a folder called `corset`


## Possibly Clean Reads

Check the quality of your sequence data with fastqc.  If your reads contain adapter sequences you may need to trim those before going forward.  If you have a large fraction of unmapped reads and your read qualities are poor you may wish to consider trimming. 

## Obtain a reference transcriptome

This might involve a whole extra step of assembling the transcriptome from your own data. Here is some good advice on how to do that [http://oyster-river-protocol.readthedocs.io/en/master/](http://oyster-river-protocol.readthedocs.io/en/master/). Alternatively, a reference might already be available.  

We assume a reference transcriptome is already available and is called `transcriptome.fasta`.  Put this reference transcriptome inside the `corset` folder.



## Map Reads to the Reference

This is done with [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).  On many HPC systems this can be loaded by doing

```bash
	module load bowtie2
```

Check that it is installed by doing

```bash
	which bowtie2
```

This should return a path showing the location of the bowtie2 program

## Step Seven - Index the Transcriptome

From within the corset directory run the `bowtie2-build` command to index the transcriptome

```bash
	bowtie2-build transcriptome.fasta transcriptome
```

## Align Reads

Now make a file called `bowtie.sh` containing the following script.  Note that `bowtie2` options used in this script were derived from [this post](https://groups.google.com/forum/#!topic/corset-project/8Je6dPQ-BFk) on the corset users forum.

This script requires both bowtie2 and samtools to be installed and available. You might need to do `module load samtools`.  The reason we need samtools is to avoid creating a large intermediate file in sam format. Instead we pipe outputs from bowtie2 directly to samtools which converts to the more compacts `bam`.

```bash
for f in *R1*.fastq.gz; do 
	r1=$f
	r2=${f/R1/R2}
	r12name=${f/R1/both}
	outname=${r12name%.fastq.gz}

	echo "Running bowtie2 on $outname and writing outputs to $outputs.bam"

	bowtie2 --no-mixed --no-discordant --end-to-end --all \
	--min-score L,-0.1,-0.1 \
	--threads 32 \
	-x transcriptome -1 $r1 -2 $r2 | samtools view --threads 4 -b - > $outname.bam
done

```


## Run Corset

For this step you will need to have [corset](https://github.com/Oshlack/Corset) installed. Corset can use sample grouping information to aid clustering. This example assumes this information isn't available and that there is just one bam file per sample.  In that case just run corset like this;

```bash
	corset *.bam
```



