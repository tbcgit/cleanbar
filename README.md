# cleanbar
#CleanBar is a flexible tool for demultiplexing reads tagged with sequentially ligated barcodes,
#accommodating variations in barcode positions and linker lengths while preventing 
#misclassification of natural barcode-like sequences and handling diverse ligation errors. 
#It also provides statistics useful for optimizing laboratory procedures.
#----------------------------------------------------------------------------------------------

#First, create the folders "res_4barcodes" and "res_23barcodes" in your linux sistem:
# mkdir res_4barcodes
# mkdir res_23barcodes

#Second, if the folders are already created, delete their contents:
# rm -r res_4barcodes
# rm -r res_23barcodes

#Third, Download the "CleanBar.c", "Atrandi_1k.fq" and "barcodes.txt" files and copy them to your working folder.

#Finally, compile the CleanBar.c program and generate the CB executable. 
# gcc  -O2 -o CB  CleanBar.c

#For information on what arguments the program accepts, type:
# ./CB  --help

#An example of how to run the program on the "Atrandi_1k.fq" file:
# ./CB  barcodes.txt Atrandi_1k.fq
