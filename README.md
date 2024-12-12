# cleanbar
 CleanBar as open-source software for demultiplexing microbial single-cell genomics data
#First, create the folders "res_4barcodes" and "res_23barcodes"
mkdir res_4barcodes
mkdir res_23barcodes

#Second, If the folders are already created, delete their contents:
rm -r res_4barcodes
rm -r res_23barcodes

#Third, compile the program CleanBar.c:
gcc  -O2 -o CB  CleanBar.c

#Fourth, run the program
./CB  -help
