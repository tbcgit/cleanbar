# CleanBar

CleanBar is a flexible tool for demultiplexing reads tagged with sequentially ligated barcodes (split-and-pool barcoding). It searches for barcodes at both ends of the sequence, accommodates variation in barcode positions and linker lengths, handles diverse ligation errors, and minimizes misclassification of natural barcode-like sequences. It also provides statistics to help optimize laboratory protocols. CleanBar was originally developed to process reads from the Atrandi platform, but it is adaptable to a wide range of barcode configurations and linker types.

## How does CleanBar work?

### Detection of complete barcode strings
On the Atrandi platform, barcodes are ligated in the order A, B, C, D. CleanBar starts its search with barcode D, as this is the last barcode added to the string. After detection of the barcode D, CleanBar begins searching for a sequence corresponding to barcode C starting from the last nucleotide of the identified barcode D, followed by a search for barcode B from the end of barcode C, and finally for barcode A from the last nucleotide of barcode B. Within each barcode set, CleanBar selects the barcode detected at the closest distance from the preceding barcode. Additionally, the program calculates the distances between the detected barcodes, corresponding to the linker length.

The expected barcode string length produced by the Atrandi platform is 44 bp, corresponding to four 8 bp barcodes and three 4 bp linkers. Barcodes do not necessarily start at the beginning of the read, as sequencing adapters may not have been trimmed properly. Additionally, the total barcode length may deviate slightly from 44 bp due to (rare) ligation errors, such as incomplete barcode strings or linkers that are not exactly 4 bp long. CleanBar is specifically designed to handle such cases and ensure reliable barcode detection.

CleanBar searches for barcodes on both ends of each read in the input FASTQ file using exact string matching, without allowing mismatches or indels. When a barcode string is found at one or both ends, it is trimmed, and the read is saved to a FASTQ file named according to the identified barcodes. For example, a read containing barcode D with label D11, barcode C  with label C9, barcode B  with label E5, and barcode A  with label C2 will be saved in a file named ``ff_D11_C9_E5_C2.fq``. If a read contains D–C–B–A barcodes on both ends, and the barcode string is detected in the reverse complement, the read is still assigned to the FASTQ file corresponding to the barcodes identified in the forward (direct) sequence - each read appears in the output only once, regardless of whether the barcode string was found on one or both ends. All sequences with complete barcode strings are saved in ``res_4barcodes`` folder. 

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

Use the "clone" command:
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

Run the ``prepare.sh`` command, which will compile the ``CleanBar.c`` program and create folders ``res_4barcodes`` and ``res_23barcodes``.  If these folders already exists, CleanBar will delete their content. 

````
bash prepare.sh
````

Each time you run CleanBar in the same directory, you must either delete these folders and create new ones manually using the commands below, or simply run the ``prepare.sh`` script again, which handles this automatically.

````
mkdir  res_4barcodes  ||  rm -r  res_4barcodes/*.fq
mkdir  res_23barcodes ||  rm -r  res_23barcodes/*.fq
````

## Running CleanBar

### How to run CleanBar:
Here we show how to run CleanBar using one of our example FASTQ files, Atrandi_1k.fq. The command is simple: it starts with ./CB, followed by the barcode list file barcodes.txt, and then your input FASTQ file:

````
./CB barcodes.txt Atrandi_4k.fastq
````

### Output files:

#### Screen output: 

By default, CleanBar displays the exact alignment of the first 2,000 sequences (this number can be adjusted using the ``-s`` flag — see the Options section below). Here we show an example of the read number 1900, which contains no barcodes in the forward (direct) sequence (dir), but four consecutive barcodes were detected in the reverse complement sequence (rc). Barcode D with label H11 was found at position 6, barcode C with label F9 at position 18, barcode B with label B5 at position 30, and barcode A with label C2 at position 42.

````
---------  1990 - dir ----------------
---------  1990 - rc  ----------------
6	 -->[BARCODE_D]
CGATCTCAACCAACAGGAATGGTGTGACTCCGAACTTGAAGGGTGACTCTTTGGTTTACCGCCGGGCTGGAGGGCAAAAATGCCTTAT
      CAACCAAC	 --> H11

18	 --> [BARCODE_C]
CGATCTCAACCAACAGGAATGGTGTGACTCCGAACTTGAAGGGTGACTCTTTGGTTTACCGCCGGGCTGGAGGGCAAAAATGCCTTAT
                  ATGGTGTG	 --> F9

30	 --> [BARCODE_B]
CGATCTCAACCAACAGGAATGGTGTGACTCCGAACTTGAAGGGTGACTCTTTGGTTTACCGCCGGGCTGGAGGGCAAAAATGCCTTAT
                              CGAACTTG	 --> B5

42	 --> [BARCODE_A]
CGATCTCAACCAACAGGAATGGTGTGACTCCGAACTTGAAGGGTGACTCTTTGGTTTACCGCCGGGCTGGAGGGCAAAAATGCCTTAT
                                          GTGACTCT	 --> C2

````


#### Folders:

In the ``res_4barcodes`` folder, you will find all demultiplexed FASTQ files wih complete barcodes. For example the FASTQ ``ff_D11_C9_E5_C2.fq`` contains three reads that contained barcod D with label D11, barcode C  with label C9, barcode B  with label E5, and barcode A  with label C2, which means that these three reads come from the same single-cell. 

In the ``res_23barcodes folder``, you will find all demultiplexed FASTQ files containing incomplete barcode strings. For example, the file ``ff_____A8_F4_B2.fq`` contains reads in which barcode D is missing, but the remaining barcodes, particularly barcode A, may still provide information, which is can be valuable in certain experimental setups. For instance, if multiple samples were combined on the same barcoding plate, the well position of the first round of barcode ligations (barcodes from the group A) may allow recovery of bulk-level data from these incomplete reads, even if single-cell resolution is lost.

#### Files: 

CleanBar also generates four additional output files, each named based on the input FASTQ file (in this example, Atrandi_4k):

``Atrandi_4k_stats.txt`` Contains the output file name for each sequence, the position of the last detected barcode, and the lengths of the linkers between barcodes.

``Atrandi_4k_summary.txt`` Lists the barcodes detected in both the forward and reverse complement sequences for each read.

``Atrandi_4k_1_0_bar.fq`` Contains reads with one or no barcodes detected.

``Atrandi_4k_3_2_bar.fq`` Contains reads with two or three barcodes detected.



## Default options
For information on what arguments the CleanBar accepts, write:

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

If you use the Atrandi platform, you do not need to modify these arguments. The default options are the following:

``-l 88``  Number of nucleotides analyzed at each end of the sequence.

``-s`` Number of reads showed on the screen by default is 2000.

``-bn 24`` Number of barcodes for group. Each group contains by default 24 barcodes.

``-bs 8`` Number of nucleotides per barcode. Each barcode is 8 nucleotides long.

``-ls 8`` Number of nucleotides per linker. The separation between barcodes (linker length) is 4 nucleotides. This is the expected value only, because CleanBar indicates the exact linker lengts in its output file. 

For example, if you want to see 50 sequences in the screen output, the command would be the following:

````
sh prepare.sh
./CB -s 50 barcodes.txt Atrandi_4k.fastq
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

