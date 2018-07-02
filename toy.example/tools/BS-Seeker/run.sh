build genome index:
python Prepocessing_genome.py -f .../ref.new.fas -p .../bowtie > .../log_Preprocessing_genome.txt

alignment:
python BS_Seeker.py -i …/toy.short.reads.fastq -t N -p …/bowtie -d …/reference_genome/ -o …/output.txt

t: tag

convert output to SAM format:
python BSSout2SAM.py -d ./reference_genome/ -f .../output_short.txt -r .../ref.new.fas
