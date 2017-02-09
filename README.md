# MATAFILER
## Overview
- Pipeline created to assemble metagenomes, profile miTags, profile functions, profile taxonomy using a variety of approaches (MATAFILER.pl)
- Build a gene catalog based on these assemblies and predicted genes, build abundance matrices from these and annotate the genes functionally (geneCat.pl)

##INSTALL
- How to install:
Most importantly set up the config.txt (symlink to Mods/MATAFILERcfg.txt). A bunch of programs needs to be present on your system, and in this file you either need to link to the dir or the executable. Follow the examples set up on my system, to see what kind of file is needed. There are comments that define for which part of the pipeline you need what programs.
In general, LCA and sdm can be obtained from the LotuS pipeline. 
programs under the header "ESSENTIAL" are essential to have to run the pipeline.
The pipeline is in alpha state, expect to have some headache getting it to run, though I will help where I can in the process. Note that some essential programs are already in the bin folder, check that all of them are executable (otherwise do 'chmod +x bin/*').
Last, some Perl libraries need to be globally available. To do this add the following to your .bashrc:
export PERL5LIB=absolute_path_to_METAFILER_dir:$PERL5LIB
where 'absolute_path_to_METAFILER_dir' is simply the dir the README.md is in, that you are now reading.

##QUICKSTART
Most importantly you need a mapping file to your files. See 'examples' dir for some map examples (also explaining how to do compound assemblies, compound mapping). These tags are important in the mapping file:
- #SmplID	Path - map always has to start with the #SmplID tag, Path is the relative path to fastq[.gz] files for each sample
- #DirPath	Base directory where subdir with the fastqs can be found
- #OutPath	Where to write the output (can be massive, make sure you have enough space)
- #RunID	The directory below OutPath, where results are stored. Also serves as global identifier for this run
- #mocatFiltPath	If for some reason you are forced to use mocat filtered fastqs and not the original, unfiltered files (strongly recommended), than you can indicate in which subdir these mocat files can be found

After this follow the sample IDs and the relative path, where to find the input fastqs.
