#!/bin/sh
#Code to create local database and look for a query of sequences. 
#fasta files you want to query need to end in .fa
inputdb= #Set to file you want as database with path
outputdb= #Set name of file you want as database 
outfile= #Set path to file of choice. 

makeblastdb -in $inputdb -dbtype nucl -out $outputdb;
for x in *fa; do blastn -db $outputdb -query $x -out $outfile -outfmt 6; done


