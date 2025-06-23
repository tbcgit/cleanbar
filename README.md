# CleanBar

CleanBar is a flexible tool for demultiplexing reads tagged with sequentially ligated barcodes, accommodating variations 
in barcode positions and linker lengths while preventing misclassification of natural barcode-like sequences and handling diverse ligation errors. 
It also provides statistics useful for optimizing laboratory procedures. 
CleanBar is a program originally designed to process reads from the Atrandi platform, yet adaptable to a wide range of barcode configurations interspersed with any type of linker.
CleanBar searches for barcodes from both ends of the sequence.

To start using Cleanbar, you need to have a file in FASTQ format to be processed. 
This file will contain barcodes at the ends of its sequences, in the direct sequence, in the complementary inverse sequence or at both ends.
Create a working folder and copy your FASTQ file into it. 
Download the CleanBar program into this folder. You should have a file with the barcodes you have used. 
If you have used the Atrandi technology, you can use the file ``barcodes.txt`` that we provide you with. 
Download it and put it in the same working folder as CleanBar and your FASTQ file.

When using our ClenBar application, all sequences in your FASTQ file with the same set of 4 barcodes are written to a single file, whose name is generated from the labels of the 4 barcodes found and will be stored in the ``res_4barcodes`` folder.

On the other hand, all analysed sequences with 2 or 3 identical barcodes are stored in a single file which will be stored in the ``res_23barcodes`` folder.

If ``file_imput.fq`` is the name of your file to be analysed, a file ``file_input_1_0_bar.fq`` is generated with all readings where no barcodes have been detected at either end, or only one BARCODE. And the file ``file_input_3_2_bar.fq`` will be generated storing all the sequences where 3 or 2 barcodes have been detected.


## How instal CleanBar:

**Installation from source:**


Clone GitHub repository of CleanBar:
````
 git clone  https://github.com/tbcgit/cleanbar
````

We access the ‘cleanbar’ work folder.
````
 cd cleanbar
````

Display the files in the working folder:

````
 dir
````

``
Atrandi_100.fq  Atrandi_4k.fastq  LICENSE    barcodes.txt     barcodes_s6.txt  prepare.sh
Atrandi_1k.fq   CleanBar.c        README.md  barcodes_20.txt  barcodes_xx.txt
``

We can copy our FASTQ file that we are going to process into this folder. 

Or use the examples we have downloaded:  *Atrandi_1k.fq* or *Atrandi_4k.fastq*


**Installation from Windows:**

If you are a Windows or Mac user and a C or C++ programmer, you can download all the files from the repository 
to your work folder and compile the CleanBar.c program on your computer with your usual C compiler.

You will be able to work under the same conditions as if you had cloned this repository with GitHub.
However, you must follow the same instructions described in GitHub.








## How to use CleanBar:

The easiest way to use ClenBar is to download the file "preparer.sh" and run it in the working folder. 
In this folder we must have the CleanBar.c file, our FASTQ file to analyse and the file with the barcodes to detect.
To do this we will launch the following command:

````
sh prepare.sh
````
If we are going to use the Atrandi technology, we can download the file ‘barcodes.txt’ from this GirHub.
After launching the "prepare.sh" script, a CB executable file will be generated that will accept at least two arguments, the file with the barcodes and the FASTQ file to be processed.
If we don't have a FASTQ file, we can use the GitHub file "Atrandi_1k.fq" or "Atrandi_4k.fastq".

For example:

````
./CB barcodes.txt Arandi_1k.fq
````

Every time we run the CB command again, it is necessary to delete all the files in the res_4barcodes and res_32_barcodes folders, 
which we can do by re-running the prepare.sh script
 	

## The script ``prepare.sh``:

The "prepare.sh" file has 3 lines of code. They are as follows:

````
mkdir  res_4barcodes  ||  rm -r  res_4barcodes/*.fq
mkdir  res_23barcodes ||  rm -r  res_23barcodes/*.fq
./CB  --help  || ( gcc  -O2 -o CB  CleanBar.c &&  ./CB  --help )
````

The first line creates the "res_4barcodes" folder. If it already exists, it will delete all its files so that it can run CleanBar and save the new files with 4 matching barcodes in this folder. If we don't delete the files in the folder, the newly processed sequences will be added to the already created files.
The second line does the same, but with the "res_23barcodes" folder, saving the reads with 2 or 3 matching barcodes in the same files.
In both cases, the file names are generated from the shared barcode labels.

The second line does the same, but in the folder ‘res_23barcodes’, saving the readings with 2 or 3 matching barcodes in the same files.
In both cases, the file names are obtained from the shared barcode labels.

The third line of this command runs the CB program with the option
 "--help" to show us how to use cleanBar.
If the CB executable, which is generated by compiling CleanBar, has not been generated this line of code will compile the ‘CleanBar.c’ program and generate the CB executable.

You can run the prepare.sh command like this:

````
sh prepare.sh
````

or 

````
bash prepare.sh
````

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

