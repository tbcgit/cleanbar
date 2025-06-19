mkdir  res_4barcodes  ||  rm -r  res_4barcodes/*.fq
mkdir  res_23barcodes ||  rm -r  res_23barcodes/*.fq

./CB  --help  || ( gcc  -O2 -o CB  CleanBar.c &&  ./CB  --help )
