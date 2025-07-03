# CleanBar

CleanBar is a flexible tool for demultiplexing reads tagged with sequentially ligated barcodes (split-and-pool barcoding). It searches for barcodes at both ends of the sequence, accommodates variation in barcode positions and linker lengths, handles diverse ligation errors, and minimizes misclassification of natural barcode-like sequences. It also provides statistics to help optimize laboratory protocols. CleanBar was originally developed to process reads from the Atrandi platform, but it is adaptable to a wide range of barcode configurations and linker types.


## How does CleanBar work?

### Detection of complete barcode strings
On the Atrandi platform, barcodes are ligated in the order A, B, C, D. CleanBar starts its search with barcode D, as this is the last barcode added to the string. After detection of the barcode D, CleanBar begins searching for a sequence corresponding to barcode C starting from the last nucleotide of the identified barcode D, followed by a search for barcode B from the end of barcode C, and finally for barcode A from the last nucleotide of barcode B. Within each barcode set, CleanBar selects the barcode detected at the closest distance from the preceding barcode. Additionally, the program calculates the distances between the detected barcodes, corresponding to the linker length.

The expected barcode string length produced by the Atrandi platform is 44 bp, corresponding to four 8 bp barcodes and three 4 bp linkers. Barcodes do not necessarily start at the beginning of the read, as sequencing adapters may not have been trimmed properly. Additionally, the total barcode length may deviate slightly from 44 bp due to (rare) ligation errors, such as incomplete barcode strings or linkers that are not exactly 4 bp long. CleanBar is specifically designed to handle such cases and ensure reliable barcode detection.

CleanBar searches for barcodes on both ends of each read in the input FASTQ file using exact string matching, without allowing mismatches or indels. When a barcode string is found at one or both ends, it is trimmed, and the read is saved to a FASTQ file named according to the identified barcodes. For example, a read containing barcode D with label D11, barcode C  with label C9, barcode B  with label E5, and barcode A  with label C2 will be saved in a file named ``ff_D11_C9_E5_C2.fq``. If a read contains D–C–B–A barcodes on both ends, and the barcode string is detected in the reverse complement, the read is still assigned to the FASTQ file corresponding to the barcodes identified in the forward (direct) sequence - each read appears in the output only once, regardless of whether the barcode string was found on one or both ends. All sequences with complete barcode strings are saved in ``res_4barcodes`` folder. 

### Detection of incomplete barcode strings
Although incomplete barcode strings do not offer single-cell resolution, some users may choose to retain reads containing 2 or 3 barcodes for downstream analysis, especially when optimizing lab protocols. If no barcode from set D is detected, CleanBar begins the search with set C. If no barcode C is found, it continues with set B. As a result, FASTQ files may be generated with only two or three barcode labels. Missing barcodes are indicated with a double underscore ``(__)``. For example, a file name  ``ff_E10_E8_C6__.fq`` would represent a read where barcode A was not detected. Complete D–C–B–A barcode strings always take priority. A read is saved based on an incomplete barcode string only if a complete barcode string is not found on the opposite end of the sequence. CleanBar can also detect barcode strings where one or two barcode sets were skipped, which may occur due to ligation errors. Demultiplexed FASTQ files with 2 or 3 identical barcodes are stored in the ``res_23barcodes`` folder. Additionally, they are all merged into a large <filename>_3_2.fastq, containing all reads with 3 or 2 barcodes. 

Reads with only one detected barcode are not saved in separate FASTQ files. This is because any 8 bp sequence has a ~0.25% chance of occurring naturally in DNA, offering no strong evidence that the match resulted from the barcoding process. Instead, these reads are grouped together with reads that contain no barcodes in a separate file named ``<filename>_0_1.fastq``, where ``<filename>`` corresponds to the name of the input FASTQ file. While sequences in this file lack single-cell resolution, they are still valuable and should not be discarded. They can support bulk-level analyses, similar to those performed in conventional metagenomics.



## What do you need to get started?

To use CleanBar, you will need:

- FASTQ file containing your sequence data from a split-and-pool barcoding assay. This file must be quality-filtered prior to use.
- ``barcodes.txt`` file listing the barcodes used in your assay.

We provide a ``barcodes.txt`` file for the Atrandi platform, based on the SINGLE-MICROBE DNA BARCODING protocol (Doc. No. DGPM02323206001). If your data were generated using this system, you can use the file directly. If you're working with data from a different split-and-pool barcoding platform, please refer to the section "How to use CleanBar with other split-and-pool barcoding platforms?" in this manual for instructions on how to generate a compatible ``barcodes.txt`` file.

  
## How to install CleanBar?

### Installation from source:

Use the "clone" command:
````
 git clone  https://github.com/tbcgit/cleanbar
````

### Setup:

Move to the cleanbar folder:

````
cd  cleanbar
````

Place your FASTQ file in the ``cleanbar`` folder. Alternatively, you can create your own working directory, but make sure that the following files are present in that folder:

- ``barcodes.txt`` (use the provided file or your own)
- ``CleanBar.c``
- ``prepare.sh``

Run the ``prepare.sh`` command, which will compile the ``CleanBar.c`` program and create folders ``res_4barcodes`` and ``res_23barcodes``.

````
bash prepare.sh
````
The compiler program gcc is usually pre-installed on most Linux systems and may also be available on macOS if developer tools are installed. If you receive an error indicating that the gcc command is not recognized, you can install it on Linux distributions using: ``sudo apt install build-essential`` . On macOS, you can install it by running: ``xcode-select --install``.

Running the ``prepare.sh`` script will automatically create two folders: ``res_4barcodes`` and ``res_23barcodes``.  If these folders already exists, CleanBar will delete their content. Each time you run CleanBar in the same directory, you must either delete these folders and create new ones manually using the commands below, or simply run the ``prepare.sh`` script again, which handles this automatically.

````
mkdir  res_4barcodes  ||  rm -r  res_4barcodes/*.fq
mkdir  res_23barcodes ||  rm -r  res_23barcodes/*.fq
````

### Installation in Windows:

If you're using Windows, you can use Git Bash, which is a Windows application that provides a Unix-like command-line environment, allowing Windows users to use Git commands and other common Unix command-line tools.


## Running CleanBar

### How to run CleanBar:
It's so simple! Here we show how to run CleanBar using our example FASTQ file ``Atrandi_4k.fastq``, which is located in the folder EXAMPLES, together with the barcode list file ``barcodes.txt``. The easiest way to test the program is to copy ``Atrandi_4k.fastq`` into the main ``cleanbar`` working directory using the following command: ``cp EXAMPLES/Atrandi_4k.fastq``.

````
./CB barcodes.txt Atrandi_4k.fastq
````

### Output files:

#### Screen output: 

By default, CleanBar displays the exact alignment of the first 2,000 sequences (this number can be adjusted using the ``-s`` flag — see the Options section below). Here we show an example of the read number 1900, which contains no barcodes in the direct sequence (dir), but four consecutive barcodes were detected in the reverse complement sequence (rc). Barcode D with label H11 was found at position 6, barcode C with label F9 at position 18, barcode B with label B5 at position 30, and barcode A with label C2 at position 42.

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

- In the ``res_4barcodes`` folder, you will find all demultiplexed FASTQ files wih complete barcodes. For example the FASTQ ``ff_D11_C9_E5_C2.fq`` contains three reads that contained barcod D with label D11, barcode C  with label C9, barcode B  with label E5, and barcode A  with label C2, which means that these three reads come from the same single-cell. 

- In the ``res_23barcodes folder``, you will find all demultiplexed FASTQ files containing incomplete barcode strings. For example, the file ``ff_____A8_F4_B2.fq`` contains reads, in which barcode D is missing, but the remaining barcodes, particularly barcode A, may still provide information, which is can be valuable in certain experimental setups. For instance, if multiple samples were combined on the same barcoding plate, the well position of the first round of barcode ligations (barcodes from the group A) may allow recovery of bulk-level data from these incomplete reads, even if single-cell resolution is lost.

#### Files: 

CleanBar also generates four additional output files, each named based on the input FASTQ file (in this example, Atrandi_4k):

- ``Atrandi_4k_stats.txt`` Contains the output file name for each sequence, the position of the last detected barcode, and the lengths of the linkers between barcodes. 

- ``Atrandi_4k_summary.txt`` Lists the barcodes detected in both the forward and reverse complement sequences for each read.

- ``Atrandi_4k_1_0_bar.fq`` Contains reads with one or no barcodes detected.

- ``Atrandi_4k_3_2_bar.fq`` Contains reads with two or three barcodes detected.

Here is an example from the ``Atrandi_4k_stats.txt`` file, specifically for read number 1990 (also shown in the screen output example above). This line indicates that the read with ID ``@m64105_240622_213943/132618/ccs`` did not contain any barcodes in its direct sequence ``1990_dir	      ------      	   0	(#0)``, but CleanBar detected four barcodes ``(#4)`` in the reverse complement sequence ``1990_rc``. The last barcode A started at position ``42``, and the four barcodes were separated by three linkers, each 4 bp long ``4 4 4``. This sequence has been saved in the file ``ff_H11_F9_B5_C2.fq``.

````
@m64105_240622_213943/132618/ccs	1990_dir	      ------      	   0	(#0)				1990_rc 	ff_H11_F9_B5_C2.fq	  42	(#4)	4	4	4
````

The following example is from the ``Atrandi_4k_summary.txt`` file, specifically for read ``@m64105_240622_213943/132618/ccs``. It shows that this read did not contain any barcodes in its direct sequence `Direct	 1990	...	 -1	...	 -1	...	 -1	...	 -1	``, but CleanBar detected four barcodes in the reverse complement sequence ``Rev_com	 1990	H11	  6	F9	 18	B5	 30	C2	 42``. The barcode labels (H11, F9, B5, and C2) are shown along with the positions of their first nucleotides (6, 18, 30, and 42).

````
@m64105_240622_213943/132618/ccs	Direct	 1990	...	 -1	...	 -1	...	 -1	...	 -1	Rev_com	 1990	H11	  6	F9	 18	B5	 30	C2	 42
````


###  Default options
For information on the arguments CleanBar accepts, run:
````
./CB  --help
````

The program will then display the following output:


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

For example, if you want to see only 50 sequences in the screen output, the command would be the following:

````
bash prepare.sh
./CB -s 50 barcodes.txt <inputFASTQ>.fastq
````



## How to use CleanBar with other split-and-pool barcoding platforms?

CleanBar's default setting are specificaly set to detect Atrandi barcodes, which consist of 8 bp sequences generated through 4 rounds of ligation, with each ligation round using groups of 24 barcodes. But what if the barcodes used in your assay have a different length than the 8 bp used in the Atrandi platform? Or what if your split-and-pool barcoding setup uses only three rounds of ligation? Or what if each group contains 96 barcodes instead of 24? We provide some examples of how to proceed:

### Using four barcode groups with a custom number of barcodes per group:
To adjust CleanBar for barcode groups of a different size, modify the ``barcodes.txt`` file and update the command with the ``-bn`` flag. We provide a sample file called ``barcodes_bn20.txt``, where each group contains only 20 barcodes. When creating your own ``barcodes.txt`` file, make sure the formatting is correct. Each group must begin with a header line such as ``#Barcode A``, ``#Barcode B``, and so on. Groups should be separated by an empty line. Each barcode entry must be on its own line, starting with a barcode label followed by the sequence, separated by a space, for example, ``A1 GTAACCGA``.

````
bash prepare.sh
./CB   -bn 20 barcodes_bn20.txt Atrandi_4k.fastq
````

### Using four barcode groups, but two groups contains a lower number of barcodes:
In this example, the barcodes are arranged as follows: group A with 24 barcodes, group B with 24 barcodes, group C with 12 barcodes, and group D with 12 barcodes. CleanBar requires each group to have the same number of barcodes, and the ``-bn`` option accepts only a single value, which must be the highest number of barcodes used in any group. Therefore, the ``barcodes.txt`` file must include placeholder barcodes, such as ``H12 XXXXXXXX``, to fill out the shorter groups. This format is illustrated in the example file ``barcodes_halfxx.txt``. If the maximum number of barcodes per group remains at the default (24), and both the barcode length and expected linker length match the default settings (8 and 4), there is no need to modify the default CleanBar command, you only need to use the modified barcode list: 

````
bash prepare.sh
./CB  barcodes_halfxx.txt Atrandi_4k.fastq
````

### Using only three barcode groups with non-default barcode and linker lengths:
If your barcodes and linkers have lengths different from the default values, you must modify the ``barcodes.txt`` file and update the command with the ``-bs`` and ``-ls`` flags. If the resulting expected total length of the barcodes string is different from default, you should also modify the search length with the ``-l`` flag. 
In this example file ``barcodes_bs12_ls6_xx.txt``, there are three groups of barcodes, individual barcodes have size of 12 bp and the expected length of linkers between them is 6 bp, this means that the expected length of the complete string is 12 + 6 + 12 + 6 + 12 = 48 bp. The string in this example is longer that the string formed by the Atrandi platform (8 + 4 + 8 + 4 + 8 + 4 + 8 = 44 bp), thus it would be convenient to increase the seach length from the default 88 to 92 by the ``-l`` flag. Please note that the ``barcodes.txt`` file must always contain four groups of barcodes. Therefore, in our example file ``barcodes_bs12_ls6_xx.txt``, the barcodes listed under the ``#Barcode D`` group contain placeholder sequences, such as ``A10 XXXXXXXXXXXX``. 

````
bash prepare.sh
./CB   -l 92 -bs 12 -ls 6 barcodes_bs12_ls6_xx.txt Atrandi_4k.fastq
````

### Combining all options:
Here is an example command line for an extreme barcoding assay, where each barcode group contains up to 96 barcodes  ``-bn 96``, the barcode size is 20 bp ``-bs 20``, and the expected linker length is 6 bp ``-ls 6``.  The assay uses four rounds of ligation, so the expected total barcode string length is calculated as: 24 + 6 + 24 + 6 + 24 + 6 + 24 = 114 bp. Since the default search length is 88 (corresponding to the default expected string length of 44), it should be increased to at least 158 (114 + 44) using the ``-l 158`` flag. In this example, we also limit the screen output to 40 sequences with ``-s 40``. The flags should be used in the following order:

````
bash prepare.sh
./CB   -s 40 -l 158 -bn 96 -bs 20 -ls 6  <barcodesfilename>.txt <inputFASTQ>.fastq
````

## Citation:
If you use CleanBar program, please, cite our article:

Vicente Arnau, Alicia Ortiz-Maiques, Juan Valero-Tebar, Lucas Mora-Quilis, Vaida Kurmauskaite, Lorea Campos Dopazo, Pilar Domingo-Calap, Mária Džunková. CleanBar: A Versatile Demultiplexing Tool for Split-and-Pool Barcoding in Single-Cell Omics. Under review in ISME Communications
