# MATAFILER
## Overview
- Pipeline created to assemble metagenomes, profile miTags, profile functions, profile taxonomy using a variety of approaches (MATAFILER.pl)
- Build a gene catalog based on these assemblies and predicted genes, build abundance matrices from these and annotate the genes functionally (geneCat.pl)
  
##INSTALLATION
Most importantly set up the config.txt (symlink to Mods/MATAFILERcfg.txt). A bunch of programs needs to be present on your system, and in this file you either need to link to the dir or the executable. Follow the examples set up on my system, to see what kind of file is needed. There are comments that define for which part of the pipeline you need what programs.  
In general, LCA and sdm can be obtained from the LotuS pipeline.  
programs under the header "ESSENTIAL" are essential to have to run the pipeline.
The pipeline is in alpha state, expect to have some headache getting it to run, though I will help where I can in the process. Note that some essential programs are already in the bin folder, check that all of them are executable (otherwise do 'chmod +x bin/*').  
Last, some Perl libraries need to be globally available. To do this add the following to your .bashrc:  
export PERL5LIB=absolute_path_to_METAFILER_dir:$PERL5LIB  
where 'absolute_path_to_METAFILER_dir' is simply the dir the README.md is in, that you are now reading.  
### temorary files
Since the pipeline is expected to run on a compute cluster, temporary directories are of enormous importance for a) performance and b) file exchange between compute nodes that are usually physically separated clusters (see below for more information on this).  
The pipeline expects a path to a storage that is globally available (given by argument "globalTmpDir") on all nodes and a tmp dir that is locally available on each node (given by arguemnt "nodeTmpDir" in config file).  

## Quickstart
### Create a mapping file for your experiment
Most importantly you need a mapping file to your files. See 'examples' dir for some map examples (also explaining how to do compound assemblies, compound mapping). These tags are important in the mapping file:
- #SmplID	Path - map always has to start with the *#SmplID* tag. This will be used in all subsequent tables, intermediary files etc to identify a sample and shoul be as short and descriptive as possible. *DO NOT USE SPECIAL CHARACTERS IN THE SMPLID, keep it basic*!  
*Path* is the relative path to fastq[.gz] files for each sample (see #DirPath, this needs to be set to the absolute path). All files ending with .fq or .fastq (can have .gz after) in the dir will be used for that specific samples. 1. or 2. indicates first or second read. E.g. al0-0_12s005629-2-1_lane3.2.fq.gz is the second read, here the pipeline expects to have al0-0_12s005629-2-1_lane3.1.fq.gz in the same dir.  
Further, you can add the following specifics for each single sample:   
*AssmblGrps* - set this to a number or string. all samples with the same tag will be assembled together (e.g. samples from the same patient at different time points).  
*MapGrps* - set a tag here as in AssmblGrps. All reads from these samples will be thrown together, when mapping against target sequences (only works with option "map2tar" and "map2DB").
*SupportReads* - in case you have additional reads, that are not normal illumina hiSeq, e.g. miSeq or hiSeq in mate pair sequence mode.
 
- #DirPath	Base directory where subdir with the fastqs can be found. You can insert this on several lines, if the base path changes for all samples afterwards.
- #OutPath	Where to write the output (can be massive, make sure you have enough space)
- #RunID	The directory below OutPath, where results are stored. Also serves as global identifier for this run
- #mocatFiltPath	If for some reason you are forced to use mocat filtered fastqs and not the original, unfiltered files (strongly recommended), than you can indicate in which subdir these mocat files can be found

After this follow the sample IDs and the relative path, where to find the input fastqs.  
See _examples/example_map_assemblies.map_ for a very complicated mapping file with several source dirs.

## Information about the pipeline
### I/O (Input/Output)
Analysing a shotgun metagenomic experiment can be a computationally extremely demanding task, as in some experiments several TB of data can be accumulating. MATAFILER was designed with the latter case in mind, but can of course also handle smaller experiment.  
In order to be able to cope with these data amounts, a lot of 'file juggling' is happening behind the scences. A lot of temporary files are being created that don't need to be saved on long term storage solution that are backed upand generally also slower. For this purpose big servers usually have a 'scratch' dir that is the global temporary storage on which all nodes in a cluster can write, but that is not backed up and might be cleaned infrequently. Further, usually each node has a local temp dir, to which only that specific node has access. Using these temporary solutions does make the whole cluster more stable and also enable other users to use a cluster more efficiently. To give you an example: if you have an IO heave process like searching with diamond through a lot of reads, you will use up the bandwith provided by your permanent storage very quickly. This could lead to situations where 500 cores on the cluster are busy with running in parallel diamond searches, but since the IO is so severely limited, only a small fraction of data trickles through to these jobs, effectively maybe giving 16 cores work. In this case the cluster would be unnecessarily blocked and the 500 core job would also take much longer than needed. That is the reason why file juggling is so important and why so much development effort went into optimizing this for MATAFILER.  
To take advantage of this, I strongly recommend to ask your sysadmin where the local and global temp storage on your cluster are and set in the MATAFILER config the variables 'globalTmpDir' and 'nodeTmpDir' variables correspondingly. 
