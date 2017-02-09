# MATAFILER
-Pipeline created to assemble metagenomes, profile miTags, profile functions, profile taxonomy using a variety of approaches (MATAFILER.pl)
-Build a gene catalog based on these assemblies and predicted genes, build abundance matrices from these and annotate the genes functionally (geneCat.pl)

#Quickstart
- How to install:
Most importantly set up the config.txt (symlink to Mods/MATAFILERcfg.txt). A bunch of programs needs to be present on your system, and in this file you either need to link to the dir or the executable. Follow the examples set up on my system, to see what kind of file is needed. There are comments that define for which part of the pipeline you need what programs.
In general, LCA and sdm can be obtained from the LotuS pipeline. 
programs under the header "ESSENTIAL" are essential to have to run the pipeline.
The pipeline is in alpha state, expect to have some headache getting it to run, though I will help where I can in the process.
