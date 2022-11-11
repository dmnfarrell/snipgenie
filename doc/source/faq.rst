FAQ
===

**The run was stopped during execution, can it be resumed?**

Yes, by default the program won't overwrite intermediate files when re-run. So just run it again.
Make sure there are no old tmp.****.bam files in the mapped folder if an alignment got interrupted.

**My sample files are not being parsed properly.**

This may be because your sample names are unusual. The program extracts the unique sample names from
the files by using the '_' symbol as delimeter. If your names differ you can supply a different
delimeter with the labelsep option.

**I added new files and tried to re-run but it failed.**

This is because the samples don't match the previous variant call output. You might see a different
number of samples warning. By default the results of mpileup are not overwritten as this is the slowest step. You should first delete the file raw.bcf in the output folder and run again.

**I have more than 1000 samples and the bcftools mpileup step fails.**

This is likely due to the limit on the number of files that can be opened at the same time. You can
increase this limit on Linux using ulimit -n 2000 or whatever value you need up to 9999. Note that for
many samples this step could take several days to run.
