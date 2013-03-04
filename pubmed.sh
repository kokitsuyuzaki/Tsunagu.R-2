#! /bin/sh
ls Put_Data_Here/ > journal.txt

while read line
do
txt=${line}
echo $txt
R --vanilla --slave --args <pubmed.R> log.txt $txt $line
done<journal.txt