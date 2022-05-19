## Simple Nextflow script to truncate FASTQ files ##

Take the first N reads of a pair of FASTQ files and output truncated FASTQ files.
This provides a simple way to get a small testset of reads

based on:
https://stackoverflow.com/questions/69010242/nextflow-how-should-i-truncate-a-fastq-file-at-line-x-process-fails-with-error

usgage example:
```
nextflow run readpickn/main.nf -with-docker [nextflow-image] --datadir [path or s3 folder uri] --outdir [path or s3 folder uri]
```

