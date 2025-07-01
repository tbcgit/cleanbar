# CleanBar

CleanBar is a flexible tool for demultiplexing reads tagged with sequentially ligated barcodes (split-and-pool barcoding). It searches for barcodes at both ends of the sequence, accommodates variation in barcode positions and linker lengths, handles diverse ligation errors, and minimizes misclassification of natural barcode-like sequences. It also provides statistics to help optimize laboratory protocols. CleanBar was originally developed to process reads from the Atrandi platform, but it is adaptable to a wide range of barcode configurations and linker types.

## How does CleanBar work?

### Detection of complete barcode strings
On the Atrandi platform, barcodes are ligated in the order A, B, C, D. CleanBar starts its search with barcode D, as this is the last barcode added to the string. After detection of the barcode D, CleanBar begins searching for a sequence corresponding to barcode C starting from the last nucleotide of the identified barcode D, followed by a search for barcode B from the end of barcode C, and finally for barcode A from the last nucleotide of barcode B. Within each barcode set, CleanBar selects the barcode detected at the closest distance from the preceding barcode. Additionally, the program calculates the distances between the detected barcodes, corresponding to the linker length.

The expected barcode string length produced by the Atrandi platform is 44 bp, corresponding to four 8 bp barcodes and three 4 bp linkers. Barcodes do not necessarily start at the beginning of the read, as sequencing adapters may not have been trimmed properly. Additionally, the total barcode length may deviate slightly from 44 bp due to (rare) ligation errors, such as incomplete barcode strings or linkers that are not exactly 4 bp long. CleanBar is specifically designed to handle such cases and ensure reliable barcode detection.

CleanBar searches for barcodes on both ends of each read in the input FASTQ file using exact string matching, without allowing mismatches or indels. When a barcode string is found at one or both ends, it is trimmed, and the read is saved to a FASTQ file named according to the identified barcodes. For example, a read containing barcode D with label E10, barcode C  with label E8, barcode B  with label C6, and barcode A  with label B2 will be saved in a file named ``ff_E10_E8_C6_B2.fq``. If a read contains D–C–B–A barcodes on both ends, and the barcode string is detected in the reverse complement, the read is still assigned to the FASTQ file corresponding to the barcodes identified in the forward (direct) sequence - each read appears in the output only once, regardless of whether the barcode string was found on one or both ends. All sequences with complete barcode strings are saved in ``res_4barcodes`` folder. 

### Detection of incomplete barcode strings
Although incomplete barcode strings do not offer single-cell resolution, some users may choose to retain reads containing 2 or 3 barcodes for downstream analysis, especially when optimizing lab protocols. If no barcode from set D is detected, CleanBar begins the search with set C. If no barcode C is found, it continues with set B. As a result, FASTQ files may be generated with only two or three barcode labels. Missing barcodes are indicated with a double underscore ``(__)``. For example, a file name  ``ff_E10_E8_C6__.fq`` would represent a read where barcode A was not detected. Complete D–C–B–A barcode strings always take priority. A read is saved based on an incomplete barcode string only if a complete barcode string is not found on the opposite end of the sequence. CleanBar can also detect barcode strings where one or two barcode sets were skipped, which may occur due to ligation errors. Sequences with 2 or 3 identical barcodes are stored in a single file which will be stored in the ``res_23barcodes`` folder.

Reads with only one detected barcode are not saved in separate FASTQ files. This is because any 8 bp sequence has a ~0.25% chance of occurring naturally in DNA, offering no strong evidence that the match resulted from the barcoding process. Instead, these reads are grouped together with reads that contain no barcodes in a separate file named ``<filename>_0_1.fastq``, where ``<filename>`` corresponds to the name of the input FASTQ file. While sequences in this file lack single-cell resolution, they are still valuable and should not be discarded. They can support bulk-level analyses, similar to those performed in conventional metagenomics.


## What do you need to get started?

To use CleanBar, you will need:

- FASTQ file containing your sequence data from a split-and-pool barcoding assay. This file must be quality-filtered prior to use
- ``barcodes.txt`` file listing the barcodes used in your assay

We provide a ``barcodes.txt`` file for the Atrandi platform, based on the SINGLE-MICROBE DNA BARCODING protocol (Doc. No. DGPM02323206001). If your data were generated using this system, you can use the file directly. If you're working with data from a different split-and-pool barcoding platform, please refer to the section "How to adapt CleanBar to other split-and-pool barcoding platforms" in this manual for instructions on how to generate a compatible ``barcodes.txt`` file.
  
## How instal CleanBar:

### Installation from source:

#### Use the "clone" command:
````
 git clone  https://github.com/tbcgit/cleanbar
````

### Installation in Windows:

If you're using Windows or macOS and have experience with C or C++ programming, you can manually download all the necessary files from this repository to your working directory and compile CleanBar.c using your preferred C compiler. Once compiled, you will be able to run CleanBar under the same conditions as if you had cloned the repository via GitHub. Simply follow the same usage instructions described in this manual.



## Setup:

Move to the cleanbar folder:

````
cd  cleanbar
````

Place your FASTQ file in the ``cleanbar`` folder. Alternatively, you can create your own working directory, but make sure that the following files are present in that folder:

- ``barcodes.txt`` (use the provided file or your own)
- ``CleanBar.c``
- ``prepare.sh``

Run the prepare command, which will compile the "CleanBar.c" program and create folders "res_4barcodes" and "res_23barcodes".  If these folders already exists, CleanBar will delete their content. 

````
bash prepare.sh
````

Each time you run CleanBar in the same directory, you must either delete these folders and create new ones manually using the commands below, or simply run the "prepare.sh" script again, which handles this automatically.

````
mkdir  res_4barcodes  ||  rm -r  res_4barcodes/*.fq
mkdir  res_23barcodes ||  rm -r  res_23barcodes/*.fq
````

## 
You can use one of our example files for testing, 
If we are going to use the Atrandi technology, we can download the file ‘barcodes.txt’ from this GitHub.
After launching the "prepare.sh" script, a CB executable file will be generated that will accept at least two arguments, the file with the barcodes and the FASTQ file to be processed.
If we don't have a FASTQ file, we can use the GitHub file "Atrandi_1k.fq" or "Atrandi_4k.fastq".

For example:

````
./CB barcodes.txt Arandi_1k.fq
````

Every time we run the CB command again, it is necessary to delete all the files in the res_4barcodes and res_32_barcodes folders, 
which we can do by re-running the prepare.sh script
 	








## An example of how to view the CleanBar options:

For information on what arguments the ``CleanBar`` program accepts, write:

````
./CB  --help
````

and the output of the programme will be:


````
===================================================
=     CleanBar : Single Cell data analysis        =
===================================================
=  Vicente Arnau & Maria Dzunkova . 19-VI-2025    =
===================================================

SINTAX: ./CB  <options>  BARCODES_File  FASTQ_File >  Screen_output_File.txt

<options>:
 -l     : Number of nt parsed at the start and the end of a read
 -s     : Number of reads showed on the screen
 -bn    : Number of barcodes for group
 -bs    : number of nucleotides per barcode (BARCODE_SIZE)
 -ls    : number of nucleotides per link (LINK_SIZE)
        : Screen_output_File.txt is optional
````
		
If you do not use any of the options provided by CleanBar, you should know that the default values are:

- The number of nucleotides analysed at the ends of the sequences is 88.
- The number of sequences that are displayed on the screen is 2000.
- Each group has 24 BARCODES.
- Each BARCODE has 8 nucleotides.
- The separation between each BARCODE is 4 nucleotides.

	

## An example of how to run the program on the ``Atrandi_1k.fq`` file:

If you do not have a FASTQ file to process, you can download for example the file ‘Atrandi_1k.fq’ from gitHub. This file has 1000 sequences with Atrandi barcodes at the ends of their sequences. At one or both ends.
To run the CleanBar program on this file with 1000 sequences we will write:

````
./CB  barcodes.txt Atrandi_1k.fq
````
 
In the ``res_4barcodes`` folder, we will have as many files as cells with 4 identical barcodes we have.
In the ``res_23barcodes`` folder, as many files will be generated as cells with 2 or 3 shared barcodes are obtained from using the ClenBar program.

We also generate 4 files with additional information about the analysis of each sequence (read) of the analyzed ``Atrandi_1k.fq`` file:
- ``Atrandi_1k_stats.txt`` --> provides the file name in which each sequence was stored, the position of the last detected barcode, and the lengths of linkers between barcodes
- ``Atrandi_1k_summary.txt``  --> lists the barcodes found in both, direct and reverse complement strings for each sequence
- ``Atrandi_1k_1_0_bar.fq``  --> Have the reads with only one or no barcodes of Atrandi_1k.fq.
- ``Atrandi_1k_3_2_bar.fq``  --> Have the reads with 3 or 3 barcodes of Atrandi_1k.fq.

The CleanBar program will display on the screen the barcodes read from the Barcodes.txt file. And then how it detects the barcodes in the first 1000 sequences (reads) of the Atrandi_1k.fq file (a user-adjustable parameter with -s flag).

## Another example of using CleanBar on the  ``Atrandi_4k.fastq`` file:


The file Atrandi_4k.fastq has 4k sequences with barcodes at their ends. 
We can process it with CleanBar, but first we must delete all the files in the "res_4barcodes" and "res_3barcodes" folders, because if we don't do it, the new sequences to process will be added to the files we already have created, if their barcodes match.
We also want to show how only the first 50 lines of Atrandi_4k.fastq are processed, so we will use the command ‘-s 50’.
That's why we will execute these two commands:

````
sh prepare.sh
./CB -s 50 barcodes.txt Atrandi_4k.fastq
````

If we want the CleanBar screen output, which includes the detail of how the barcodes are detected in the first 50 lines, we will run CleanBar as follows, generating a text file "screen.txt" where this information will be saved:


````
sh prepare.sh
./CB -s 50 barcodes.txt Atrandi_4k.fastq  >  screen.txt
````


## Example of using CleanBar with a file with 20 barcodes per set.

We can download from GitHub the file "barcodes_20.txt" which has only 20 barcodes in each set. 
If we want to process the example file "Atrandi_4k.fastq" using this barcodes file, we will have to tell the program on the run line that the barcodes file to be used now has only 20 barcodes. We will do it using the option "-bn 20" and the file "barcodes_20.txt.
Also, we only want to display 60 lines per screen. We will indicate it with the option "-s 60".
This way:


````
sh prepare.sh
./CB   -s 60 -bn 20 barcodes_20.txt Atrandi_4k.fastq
````
# How to adapt CleanBar to other split-and-pool barcoding platforms
## Example of using CleanBar with a file with only 3 sets of barcodes.


Let's imagine that we now have an experiment that has made use of only 3 sets of barcodes. 
And the third set, instead of having 24 barcodes, has only 20 barcodes.
We can use the file "barcodes_xx.txt" that has the last set with the barcodes defined as "XXXXXXXXXX", so they will not be found in the sequence extresms.
And as the third set of barcodes has only 20 barcodes defined, the remaining 4 barcodes will also be replaced by barcodes defined as "XXXXXXXXXX".
This will be the content of the file "barcodes_xx.txt":

````
#Barcode A
A1 GTAACCGA
A2 TCCTCAAC
A3 TGGTCTCA
B1 GACAGCAT
. . .
H2 ATGTCTGC
H3 ACAATCCG

#Barcode B
A4 TACAACCG
A5 GCTGGATA
A6 CATCGTTG
. . .
H4 GGAACTGT
H5 TATGCGAC
H6 TGTTGGAC

#Barcode C
A7 TACAGCAG
A8 TTCGGTAG
A9 GATACCGA
 . . .
G7 CGCTACTA
G8 AAGGTGAC
G9 XXXXXXXX
H7 XXXXXXXX
H8 XXXXXXXX
H9 XXXXXXXX

#Barcode D
A10 XXXXXXXX
A11 XXXXXXXX
A12 XXXXXXXX
 . . .
N11 XXXXXXXX
H12 XXXXXXXX
````


## Another example:

This is an extreme example of the use of CleanBar.
Now we have a barcode file with only 5 barcodes per set and the barcodes are only 6 nucleotides long.
And we want to analyse only the first 70 nucleotides of the test Fastq file. With link sequences of 6 nucleotides.
And displaying on screen only the first 40 lines processed.
We should run the following command:


````
sh prepare.sh
./CB   -s 40 -l 70 -bn 5 -bs -ls 6  barcodes_s6.txt Atrandi_4k.fastq
````

As we can see, CleanBar is a very flexible programme that can be adapted to a wide variety of use cases.

