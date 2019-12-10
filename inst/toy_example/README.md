This directory contains some artificial data. Its purpose is to serve
as a toy example for the HMVAR package.

The directory **merged.snps** represents the hypothetical results
of calling bi-allelic SNPs on one genome, across 10 samples, from two
groups. The format is the same as the one produced by MIDAS after
calling `midas_merge.py snps`.

The file **map.txt** associates each hypothetical sample with
one of two groups.

The directory **expected** contains the expected results from
different functions and analysis.

At some point, the data here might also be included as R data objects,
however, it is included as flat files in order to fully test the I/O
capabilities of HMVAR.
