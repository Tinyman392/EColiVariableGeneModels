# Accessory Genome Prediction

## Requirements

There are quite a few requirements for this repository.  Bash shell and Python scripts are used heavily in many of the pipelines.    You'll need two external tools along with various Python 2 tools (version 2.7.16 used).  The external tools are KMC and GenomicModelCreator.  

KMC can be downloaded [here](https://refresh-bio.github.io/software/) (version 2.3.0 used).  You'll want the KMC executables in your PATH variable in your active terminal.  

The GenomicModelCreator can be downloaded [here](https://github.com/Tinyman392/GenomicModelCreator/tree/edf2c58a616a2bdb7bbaa19aee87dd7b795afcba).  Although you can add the buildModel.py script to your PATH variable, it is not required.    GenomicModelCreator has it's own set of pre-requisites which should be installed:
- numpy ([website](https://numpy.org), [anaconda](https://anaconda.org/anaconda/numpy)) [version 1.16.5 used]
- sklearn ([website](https://scikit-learn.org/stable/), [anaconda](https://anaconda.org/anaconda/scikit-learn)) [version 0.20.3 used]
- scipy ([website](https://www.scipy.org), [anaconda](https://anaconda.org/anaconda/scipy)) [version 1.2.1 used]
- xgboost ([website](https://xgboost.readthedocs.io/en/latest/), [anaconda](https://anaconda.org/conda-forge/xgboost)) [version 0.90 used]

The kmc.sh script from this repository should be in your PATH variable.  

The models on the BV-BRC (formerly PATRIC) FTP are in a zip file, so you'll need the zip command installed on your machine.  Precomputed models can be downloaded <a href="ftp://ftp.bvbrc.org//datasets/Nguyen_et_al_2023.zip">here</a>.

The prediction script relies on finding the conserved genes used to train the model in an annotated fasta file.  This requires the use of an annotated genome on the [BV-BRC](https://www.bv-brc.org).  This can be done through the web browser or through the [BV-BRC's command line interface](https://www.bv-brc.org/docs/cli_tutorial/index.html).  

## Contents

- downLoadFTP.sh : Downloads genomes from the BV-BRC FTP to working directory.
- genomes.gid.lst : List of E. coli genomes to download
- getClusters.py : Script clusters genomes' KMC output to generate a diverse set of 500, 1000, 2000, and 4000 genomes.
- getSubsample.sh : Script will subsample the genomes from the *getClusters.py* to get a small diverse set of genomes.  
- kmc.sh : runs kmc and kmc_dump tools in tandem and outputs tab delimited file
- parseFTP.py : Parses downloaded genomes from FTP to produce conserved and accessory genes.  
- plf.con.lst : File containing a list of PLFs that the precomputed models are trained from.  
- predict.py : This script will make a prediction for a genome given it's annotation tabular file, a list of conserved PLFs, and directory of models.  
- README.md : This file
- runKMC.sh : Runs KMC on all fasta files in a directory
- trainPipeline.sh : Runs the training pipeline to build new models.

## Training Models

A pre-computed set of models can be downloaded from <a href="ftp://ftp.bvbrc.org//datasets/Nguyen_et_al_2023.zip">the BV-BRC FTP site]</a>.  This first pipeline is designed to allow you to retrain the same models or even train your own with some modification.  If you wish to just predict on pre-computed models, please go see the [Predicting with models](#predicting-with-Models) section of this README.  

To train your own models, you must download the genomes used for the paper, parse the downloaded genomes, run KMC on the set of conserved genes, cluster the genomes, then finally train on the dataset.  This is done with the following scripts:
- downloadFTP.sh
- parseFTP.py
- runKMC.sh
- getClusters.py
- buildModel.py (*external*)

#### downloadFTP.sh

This script will download a set of genomes from the BV-BRC FTP database.  This script takes as an argument a file containing a list of genome IDs to download from the FTP database.  This list of genome IDs can be found in the *genomes.gid.lst* file located in this repository.  

This script will output an ftp directory in the current working directory.  **If the directory exists already, that directory will be deleted and recreated.**

Note that the script will end up downloading all 34k E. coli genomes from the BV-BRC.  This results in 471GB of data being downloaded!

This script can be run as follows:

``` bash
downloadFTP.sh PATH/TO/REPOSITORY/genomes.gid.lst
```

#### parseFTP.py

This script will parse the BV-BRC ftp directory created by the downloadFTP.sh script.  It takes in various arguments as described below:
- -f | --ftp : Specify the location of *genomes* directory in the *ftp* directory that was created by the *downloadFTP.sh* script.  This is a required parameter.
- -g | --gid : Specify a file containing a list of genome IDs to filter off of.  This is an optional parameter
- -o | --out_pref : Specify the prefix for output files and directories.  The output will be *out_pref*.acc.plf/, *out_pref*.con.fasta/, *out_pref*.cnts.tab, *out_pref*.plf.con.lst.  This is an optional parameter; the default value is "outPref".
- -n | --n_plf : Specify the number of top PLFs to get per genome.  This is an optional parameter; the default value is 100.  

``` bash
parseFTP.py -f /PATH/TO/DOWNLOADED/ftp/genomes -o out_pref
```

#### runKMC.sh

This script will run the KMC tool on all fasta files contained within a specified directory.  This is to be used on the fasta directory specified by the *-o | --out_pref* option from the *parseFTP.py* script.  Specifically you'll want to use the *out_pref.con.fasta/* directory (replace "out_pref" with whatever you specified with the -o|--out_pref option).  It also takes an additional output directory for all the KMC outputs.  

``` bash
runKMC.sh /PATH/TO/out_pref.con.fasta/ output_directory
```

#### Step: run KMC on all genomes (concatenated)

To get clusters, you need to concatenate all the fasta files then run KMC on that concatenated file.  This allows us to get a master list of k-mers.  

``` bash
cat /PATH/TO/runKMC/OUTPUT/*.fasta > allFasta.fasta
kmc.sh 7 allFasta.fasta allFasta.fasta . 1
```

#### getClusters.py

This script will subsample out differ numbers of genomes based on a hierarchical clustering algorithm.  It will get clusters of 500, 1000, 2000, and 4000.  It takes in 3 parameters:
- KMC directory : Specify the KMC directory from *runKMC.sh*
- All fasta KMC : From the step above, we specified allFasta.fasta and this outputted a file named allFasta.fasta.7.kmrs, this is what you want to specify for this parameter.  
- Output prefix : Specify the output prefix.  This script will output multiple directories, one for each of the clusters created.  The directories will be named *out_pref.XXXX.clusts* where XXXX is the size of the cluster.

``` bash
python getClusters.py /PATH/TO/runKMC/OUTPUT/ PATH/TO/allFasta.fasta.7.kmrs out_pref
```

#### getSubsample.sh

From here will create individual directories for the cluster directories, one for fasta file and one for KMC files.  

``` bash
bash getSubsample.sh /PATH/TO/out_pref.con.fasta/ /PATH/TO/out_pref.con.kmc/ /PATH/TO/out_pref.XXXX.clusts
```

This will output two directories named *out_pref.con.kmc.XXXX* and *out_pref.con.fasta.XXXX*.  

#### buildModel.py

This script is part of the [GenomicModelCreator](https://github.com/Tinyman392/GenomicModelCreator/tree/edf2c58a616a2bdb7bbaa19aee87dd7b795afcba) github.  The output from the *parseFTP.py* script included a directory of accessory PLF tabular files named *out_pref*.acc.plf/.  We will need to build one model per file in that directory.  The loop below should be able to train all these models.  Note that each model (4000 genomes, 100 PLFs) took 4 minutes per model on a machine utilizing with 128 Intel Xeon Gold 6148 cores.  With a total of over 3000 models, this took us over a week to train.

You'll be specifying the following parameters for the model building:
- -f : Fasta directory, specify the one from the getSubsample.sh output.
- -t : Tabular file to train off of.  This is the tabular file that is in the *out_pref.acc.plf* directory.  
- -T : Temp directory.  Any directory can be specified, **note that this directory may be cleared prior to training!**
- -o : Output directory for the model.  It is useful to organize all the model outputs into a directory of their own.  
- -n : Number of threads.  
- -d : Depth of 16, as used in the paper.
- -k : k size of 7 for KMC and k-mer counting.  
- -K : KMC output directory.  Specify the one from the getSubsample.sh output.  This isn't required, but if building multiple models (like 3000+ models), not having to run KMC on all genomes 3000 times will save time in the long run.
- -c : Specify that the model is a classifier (vs regressor)
- -w : Weight the training inputs by class distribution to reduce the impact of class imbalance.
- -W : Specify that the data are not MIC models or listed in powers of 2 that need to be linearized.  The buildModel script was originally designed for MIC prediction, so this parameter is used to turn that off.  

``` bash
mkdir /PATH/TO/MODEL/OUTPUT/

for i in /PATH/TO/out_pref.acc.plf/*.tab; do
	python /PATH/TO/GENOMIC/MODEL/CREATOR/buildModel.py -f /PATH/TO/out_pref.fasta.XXXX/ -t $i -T temp/ -o /PATH/TO/MODEL/OUTPUT/$(basename $i) -n N_THREADS -d 16 -k 7 -K /PATH/TO/KMC/DIRECTORY/ -c True -w True -W False
done
```

This small loop should train all the models, one at a time.  It is quite computationally intensive depending on the number of PLFs being used and the number of genomes used.  

#### trainPipeline.sh

This script will run the pipeline from downloading the genomes from FTP all the way through training the models.  It takes two parameters:
- Path to the GenomicModelCreator directory
- Number of max cores to use during training

``` bash
bash trainPipeline.sh /PATH/TO/GenomicModelCreator/ 128
```

Going through the script you'll see the entire pipeline run in order.  This script is set up to only run one of the models instead of all 3000+.  Commenting out the *break* statement on will allow the script to run all 3000+ models.  

## Predicting with Models

#### Download models from BV-BRC FTP

The pre-computed models are stored on the BV-BRC FTP in a 5.4GB zip file.  You'll need to download the file then unzip it.  It is located [here](ftp://ftp.bvbrc.org//datasets/Nguyen_et_al_2023.zip).  You can also download it and unzip it using the curl and unzip commands below.

``` bash
curl ftp://ftp.bvbrc.org//datasets/Nguyen_et_al_2023.zip > models.100.4000.zip
unzip models.100.4000.zip
```

#### Annotate a Genome

Use the [comprehensive genome analysis](https://www.bv-brc.org/app/ComprehensiveGenomeAnalysis) or [genome annotation](https://www.bv-brc.org/app/Annotation) tool on the BV-BRC to annotate a genome.  This can be done using the web interface or the CLI's [p3-submit-genome-annotation](https://www.bv-brc.org/docs/cli_tutorial/command_list/p3-submit-genome-annotation.html) (the comprehensive genome analysis is not yet supported on the CLI).  Note that the genome annotation service requires an assembled genome.  BV-BRC offers a [genome assembly](https://www.bv-brc.org/app/Assembly2) service as well as a tool on the CLI to do it as well called [p3-submit-genome-assembly](https://www.bv-brc.org/docs/cli_tutorial/command_list/p3-submit-genome-assembly.html).

You'll need some files from the annotation.  

If you did a comprehensive genome analysis, you'll find one file named as:
- annotation/annotation.txt : this is the feature tabular file that the *predict.py* script will use.

If you did a genome annotation, you'll find one file named as:
- [genome name].txt : this is the feature tabular file that the *predict.py* script will use.  

Note, if you have an old annotation of the genome, the txt file may not contain the PLFam column and the annotation job may need to be rerun.

Once you've annotated the genome, you could use the [p3-cp](https://www.bv-brc.org/docs/cli_tutorial/command_list/p3-cp.html) command from the [BV-BRC Command Line Interface](https://www.bv-brc.org/docs/cli_tutorial/index.html) to download the annotation file or download it through the BV-BRC web interface.  

#### predict.py

Once your genome has been annotated and you've pulled the appropriate files, you can use our prediction software.  It takes the following parameters:
- -f | --feature_tab : Feature tabular file that you downloaded from your annotated genome (step above).  
- -p | --cons_plfs : file containing the list of PLFs used to train the model.  For our pre-computed models, it is named *plf.con.lst* in the root directory of this Github.
- -m | --models : directory containing all the models which will be predicted on.  Our precomputed models are on the BV-BRC FTP.  If you trained your own models using the steps above, you would have had to specify your models directory in the last step when training.  
- -t | --temp_dir : Temporary directory to hold fasta files, KMC output, etc.  This directory may be cleared when running the script!  The default value for this is *temp/*.  

``` bash
python predict.py -f [feature tab file] -p /PATH/TO/plf.con.lst -m /PATH/TO/MODELS/
```

This script will output, to standard output a tab-delimited table with two columns:
- PLFam
- Prediction for presence or absence of the PLFam

For the second column, there are 4 possible outcomes:
- Y  : At least 4 of the 5 folds for the model predicted that the PLF should exist in the genome.
- Y* : 3 of the 5 folds for the model predicted that the PLF should exist in the genome.  
- N  : At least 4 of the 5 folds for the model predicted that the PLF should not exist in the genome.
- N* : 3 of the 5 folds for the model predicted that the PLF should not exist in the genome.  
