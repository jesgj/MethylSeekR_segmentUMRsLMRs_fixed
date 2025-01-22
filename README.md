# MethylSeekR_segmentUMRsLMRs_fixed
segmentUMRsLMRs function (MethylSeekR) fixed

A member of my lab was using the R package **MethylSeekR**, but encountered an error when running the function.**segmentUMRsLMRs()**:

```
Error in h(simpleError(msg, call)) :
  error in evaluating the argument 'subject' in selecting a method for function 'vcountPattern': unused argument (as.character = FALSE)
```
The problem is in the function vcountPattern() that belongs to the package Biostrings so I fixed it.

Versions of the packages causing this issue:

- Biostrings_2.74.1
- MethylSeekR_1.46.0
