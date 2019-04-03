# read.aligment.review
read.aligment.review

Test data can be found in read.alignment.review/toy.vx.example/toy.test.data

toy.reference.fa is a reference fasta file 10,000 bp long.

toy.short_reads.fq is a fastq file with 10 reads each 100 bp long. The reads are substrings of toy.reference.fa

toy.long_reads.fq is a fastq file with 10 reads each 2000 bp long. The reads are substrings of toy.reference.fa

The script move_test_data.sh simply copies the above toy files in toy.test.data to a specified folder.
i.e. ./move_test_data.sh cushaw3 //copies toy.reference.fa toy.short_reads.fq toy.long_reads.fq to the cushaw3 folder.
