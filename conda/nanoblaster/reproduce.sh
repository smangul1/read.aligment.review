#!/bin/bash
nanoblaster -C10 -r reference.fasta -i reads.toy.example.fastq -o reads.toy.example.nanoblaster

#Creates a sam file. The output sam is missing the appropriate headers(?), so I couldn't convert to bam. 
# Error message: [samopen] no @SQ lines in the header.[sam_read1] missing header? Abort!
