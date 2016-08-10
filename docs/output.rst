.. _sec-output:

==========================
Making sense of the output
==========================


************
Output files
************

``dEploid`` outputs text files with user-specified prefix with flag **-o**.

*prefix*.log
    Log file records ``dEploid`` version, input file paths, parameter used and proportion estimates at the final iteration.

*prefix*.llk
    Log likelihood of the MCMC chain.

*prefix*.prop
    MCMC updates of the proportion estimates.

*prefix*.hap
    Haplotypes at the final iteration in plain text file.

*prefix*.vcf
    When flag ``-vcfOut`` is turned on, haplotypes are saved at the final iteration in VCF format.

*prefix*.single*i*
    When flag ``-exportPostProb`` is turned on, posterior probabilities of the final iteration of strain *i*.

***********************
Interpreting the output
***********************

.. todo::
    Include r scripts calling for interpreting the data
