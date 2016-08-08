.. _sec-output:

==========================
Making sense of the output
==========================


************
Output files
************

``DEploid`` outputs text files with user-specified prefix with flag **-o**.

*prefix*.log
    Log file records _DEploid_ version, input file paths, parameter used and proportion estimates at the final iteration.

*prefix*.llk
    Log likelihood of the MCMC chain.

*prefix*.prop
    MCMC updates of the proportion estimates.

*prefix*.hap
    Haplotypes at the final iteration in plain text file.

*prefix*.vcf
    Haplotypes at the final iteration in VCF format.


***********************
Interpreting the output
***********************

.. todo::
    Include r scripts calling for interpreting the data
