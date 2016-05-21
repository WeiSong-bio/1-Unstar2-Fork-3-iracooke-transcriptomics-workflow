# Transcriptomics Workflow

Instructions to take raw transcriptomics reads and progress to a counts file that can be used in `R` packages for statistical analysis of differential expression.

To follow this tutorial you will need to have some basic familiarity with working in a unix-like environment.  An excellent place to start would be to attend a [software carpentry](http://software-carpentry.org/) or [data carpentry](http://www.datacarpentry.org/) workshop.  The VLSCI also have several tutorials for [working with unix and/or a HPC environment](http://vlsci.github.io/lscc_docs/tutorials/)

## Step One - Connect to analysis server
This will usually mean opening a terminal and using `ssh` to connect. The general format of this command is;

```bash
	ssh username@my.hpc.server.edu.au
```

## Step Two - Download your files
Often the sequencing center will provide a url to download raw sequencing data do that a tool like wget can be used to do the download. For example;

```bash
	wget -O SequencingData.tar https://data.sequencingcenter.edu.au/download?key
```

Once you have downloaded the data it is a good idea to keep the original unmodified files somewhere. 

## Step Three - Unpack the files

How to do this depends on how the files were packaged up. The general tool to use is `tar`, possibly with options to unzip.  Assuming the downloaded data is in a file called `Project_ANDI1971.tar` we would use this command to unpack.

```bash
	tar -xvf Project_ANDI1971.tar
```

For gzipped files (`.gz` extension) you can use the `-z` option to unzip and untar at the same time. 

```bash
	tar -zxvf Project_ANDI1971.tar.gz
```

Now move all these raw reads into a folder called `raw_data`

```bash
	mkdir raw_data
	mv Project_ANDI1971 raw_data/
```

## Step Four - Possibly Clean Reads

Check the quality of your sequence data with fastqc.  If your reads contain adapter sequences you may need to trim those before going forward.  

## Step Five - Obtain a reference transcriptome

This might involve a whole extra step of assembling the transcriptome from your own data.  Alternatively, a reference might already be available.  We assume a reference transcriptome is already available and is called `Trinity.fasta`.  Put this reference transcriptome inside the `raw_data` folder.


## Step Six - Organise Files

This is done with [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).  On many HPC systems this can be loaded by doing

```bash
	module load bowtie2
```

Check that it is installed by doing

```bash
	which bowtie2
```

This should return a path showing the location of the bowtie2 program

Now make a directory where we are going to do all the alignment and read counting.  Change into this new directory

```bash
	mkdir corset
	cd corset
```

Run the following script to make symbolic links (aliases) for all your raw data files inside the `corset` directory


```bash
file_paths=`find ../raw_data/ -name *.fastq.gz`

for f in $file_paths; do 
	link_name=`basename $f`	
	ln -s $f $link_name
done
```

Also make a symbolic link for the reference transcriptome

```bash
	ln -s ../raw_data/Trinity.fasta transcriptome.fasta
```

## Step Seven - Index the Transcriptome

From within the corset directory run the `bowtie2-build` command to index the transcriptome

```bash
	bowtie2-build transcriptome.fasta transcriptome
```

## Step Seven - Align Reads

Now make a file called `bowtie.sh` containing the following script.  Note that `bowtie2` options used in this script were derived from [this post](https://groups.google.com/forum/#!topic/corset-project/8Je6dPQ-BFk) on the corset users forum.

This script requires both bowtie2 and samtools to be installed and available. You might need to do `module load samtools`.  The reason we need samtools is to avoid creating a large intermediate file in sam format. Instead we pipe outputs from bowtie2 directly to samtools which converts to the more compacts `bam`.

```bash
for f in `ls *R1*.fastq.gz`; do 
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


