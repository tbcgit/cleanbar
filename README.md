# CleanBar

CleanBar is a flexible tool for demultiplexing reads tagged with sequentially ligated barcodes, accommodating variations in barcode positions and linker lengths while preventing misclassification of natural barcode-like sequences and handling diverse ligation errors. 
It also provides statistics useful for optimizing laboratory procedures. 
CleanBar is a program originally designed to process reads from the Atrandi platform, yet adaptable to a wide range of barcode configurations interspersed with any type of linker.
CleanBar searches for barcodes from both ends of the sequence.


All sequences in which the same set of 4 barcodes has been detected are written to a single file, whose name is generated from the labels of the 4 barcodes found. 
All these generated files are stored in the ``res_4barcodes`` folder.


On the other hand, all sequences analyzed with 2 or 3 identical barcodes are stored in different files in the ``res_23barcodes`` folder. 
The file names are generated using the 3 or 2 labels shared by the sequences.


Additionally, a file ``<file_input>_1_0_bar.fq`` is generated with all the reads where no barcodes have been detected at either end.
And ``<file_input>_3_2_bar.fq`` stores all the sequences in which 3 or 2 barcodes have been detected.

## How to used:
First, create the folders ``res_4barcodes`` and ``res_23barcodes`` in your linux sistem:

````
mkdir res_4barcodes
mkdir res_23barcodes
````

Second, if the folders are already created, delete their contents:

````
rm -r res_4barcodes
rm -r res_23barcodes
````

Third, Download the "CleanBar.c", "Atrandi_1k.fq" and "barcodes.txt" files and copy them to your working folder.

Finally, compile the CleanBar.c program and generate the CB executable using the program ggc, present in all the linux systems 

````
gcc  -O2 -o CB  CleanBar.c
```` 

For information on what arguments the program accepts, write:

````
./CB  --help
````

## An example of how to view the CleanBar options:

````
./CB  --help
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


