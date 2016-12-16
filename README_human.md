# Transcriptomics Workflow

Instructions to take raw transcriptomics reads and progress to a counts file that can be used in `R` packages for statistical analysis of differential expression.

To follow this tutorial you will need to have some basic familiarity with working in a unix-like environment.  An excellent place to start would be to attend a [software carpentry](http://software-carpentry.org/) or [data carpentry](http://www.datacarpentry.org/) workshop, or complete the [Unix bootcamp](
http://rik.smith-unna.com/command_line_bootcamp/?id=9xnbkx6eaof).  The VLSCI also have several tutorials for [working with unix and/or a HPC environment](http://vlsci.github.io/lscc_docs/tutorials/)

## Connect to Analysis Server
This will usually mean opening a terminal and using `ssh` to connect. The general format of this command is;

```bash
	ssh username@my.hpc.server.edu.au
```
## Adding Files into a HPC Terminal
To begin, navigate to the directory containing the files you want to upload to the server. A new terminal should appear on the screen.
The next step is to transfer the files from the new terminal into the HPC terminal. This step is completed using the 'secure copy' or `scp` command and is as follows:

```bash
	scp -r /path/to/the/data/ username@my.hpc.server.edu.au:~/
```

If you are unsure of the path to the data, use the `pwd` command in the terminal, as this will display the path to your current file/directory.
*Note the space either side of the `-r` and before the username*


*After downloading, it is a good idea to keep the original, unmodified files somewhere safe.*

#### Combining Data from Different Files
*Depending on the location of your files, you may not need to complete this step*
It may be the case that your data files are located in separate files within one directory \(my_directory\). Start by finding all of the required files by using the `find` command. This will ensure that you have all the correct files when you attempt to link them. The command is as follows:

```bash
	find path/to/files/ -name '*common.file.name'
```

The above command tells the cluster to find, within a particular directory, all the files that contain that particular name. You may wish to change the search to find all files that end with a particular file extension. If this is the case you would change '\*common.file.name' with '\*.file.extension'.

The output for this command will show a list of the files that contain the common name that you specified above.

### Creating a Loop to Link Files
A loop is a tool that enables you to execute a set of commands repeatedly.
This step uses the loop to take the original files and create a link with them into another directory. This allows one version of the files to be modified, while the originals remain untouched. The command for this is as follows:

```bash
	for f in `find path/to/files -name '*common.file.name'`; do ln -s $f .; done
```

This command says that for all the files listed in the step above, link them to the current directory.  
This creates a link between one directory \(i.e. the directory that contains the original files\) and your current directory \(indicated by the period.\) which allows you to access all the files in the original directory in your current directory, without changing anything about the original files.

This process is much more efficient than using the `cp` command to copy all the files into a new directory. 

## Downloading the Human Genome
This step will be important when using the HISAT program to align the RNAseq reads onto the genome. This step will use the `wget` command, which allows the user to download particular files from the web. The general format of this command is:

```bash
	wget http://website.come/files
```
To download a version of the *human genome* with HiSat2 indexes, enter the following command into the terminal on HPC

```bash
	wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
```

### Unpacking the Genome File
The human genome file given above is a compressed file and therefore needs to be unpacked.  *This step will be different if you have downloaded a different version of the human genome*.
To unpack the version given above, enter the command as follows:

```bash
	tar -zxvf grch38.tar.gz
```

## Load the HISAT2 Program
This is done using the HISAT2 program, which can be loaded onto the cluster using the command:

```bash
	module load hisat2
```

You may check that the program has been installed by entering:

```bash
	which hisat2
```

This should return a path showing the location to the hisat2 program.

Using this same process load the \'samtools\' program.
