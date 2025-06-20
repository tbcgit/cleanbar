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

./CB barcodes.txt Arandi_1k.fq

Every time we run the CB command again, it is necessary to delete all the files in the res_4barcodes and res_32_barcodes folders, 
which we can do by re-running the prepare.sh script
 

For information on what arguments the ``CleanBar`` program accepts, write:

````
./CB  --help
````

## An example of how to view the CleanBar options:

````
./CB  --help
````

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
		
		

## The script ``prepare.sh``:

You can download the "prepare.sh" script. 
When you run it, if the ``res_4barcodes`` and ``res_23barcodes`` folders aren't created, it will create them.
And if they're already created, it will delete all the files inside them.
Next, it will show you how to run CleanBar using the CB executable.


You can run the prepare.sh command like this:

````
bash prepare.sh
````

## An example of how to run the program on the ``Atrandi_1k.fq`` file:

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

````
./CB   -s 50 -l 78 -bn 5 -bs 6  barcodes_s6.txt Atrandi_4k.fastq
````

Now we're going to use CleanBar with a file containing 4 groups of 5 barcodes, each containing barcodes of only 6 nucleotides. 
It will also only display the first 50 sequences analyzed. 
We will only analyze the first 78 nucleotides of each sequence.
