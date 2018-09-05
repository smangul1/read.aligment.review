bowtie2 -x reference.fasta -U reads.toy.example.fastq | samtools view -bS - >reads.toy.example.bowtie2.bam
