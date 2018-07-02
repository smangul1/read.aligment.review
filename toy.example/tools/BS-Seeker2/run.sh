build genome index:
python bs_seeker2-build.py -f …/ref.new.fas --aligner=bowtie -p …/bowtie

alignment:
bs_seeker2-align.py -I …/toy.short.reads.fastq --aligner=bowtie -p …/bowtie -o …./output.sam -f sam -g …/ref.new.fas
