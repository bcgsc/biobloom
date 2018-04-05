#!/bin/bash
#Placeholder script for generating .rv (reversed files) needed for indexing a present time
#Requires gnu awk, tac, tr, and seqtk to be in your path
#usage genRV.sh filename.fa
for filename; do
	cat $filename | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | tac | tr "\t" "\n" |seqtk seq -r - > $filename.rv
	echo $filename".rv written"
done
