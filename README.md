## Merge VCFs

**Note!** Please use a real tool for this such as https://github.com/mkirsche/Jasmine.  If you must use something along the lines of this tool, an improved and better maintained fork of this code can be found at https://github.com/papaemmelab/mergeSVvcf.

This set of routines merges VCF files by calls, labelling calls by the callers that
made them in an `INFO` field, `Callers=`.   Most of the work goes into normalizing
SV calls, which are treated as paired of breakpoints that are equal if they fall within
a user-specified window.

Here’s a very quick example of merging three VCFs; 

```bash
mergevcf -l broad,dkfz,sanger -n -m 2 svs_broad.vcf svs_dkfz.vcf svs_sanger.vcf > merged_svs.vcf
```

where each SV is labeled by the caller that saw it, with the labels given (`-l broad,dkfz,sanger`) and the number of callers that saw it (`-n`), at least two have to see the breakpoint for it to be PASS (`-m `2), and the vcf files are given.

An overview of how it works can be found on the [Simpsonlab blog](http://simpsonlab.github.io/2015/06/15/merging-sv-calls/).
