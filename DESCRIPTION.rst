Merge VCFs
==========

This set of routines merges VCF files by calls, labelling calls by the callers that
made them in an INFO field, Callers=.   Most of the work goes into normalizing
SV calls, which are treated as paired of breakpoints that are equal if they fall within
a user-specified window.

