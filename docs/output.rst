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

*prefix*.single[i]
    When flag ``-exportPostProb`` is turned on, posterior probabilities of the final iteration of strain [i].


******************************
Example of output interpretion
******************************


Example 1
*********

::

    $ ./dEploid -vcf data/exampleData/PG0390-C.eg.vcf.gz \
    -plaf data/exampleData/labStrains.eg.PLAF.txt \
    -noPanel -o PG0390-CNopanel
    $ R --slave "--args -vcf data/exampleData/PG0390-C.eg.vcf.gz
    -plaf data/exampleData/labStrains.eg.PLAF.txt
    -dEprefix PG0390-CNopanel
    -o PG0390-CNopanel " < utilities/interpretDEploid.r

.. image:: _static/PG0390-CNopanel.interpretDEploidFigure.1.png
   :width: 1024px
   :alt: interpretDEploidFigure.1

The top three figures are the same as figures show in :ref:`data example <sec-eg>`, with a small addition of inferred WSAF marked in blue, in the top right figure.

- The bottom left figure show the relative proportion change history of the MCMC chain.
- The middle figure show the correlation between the expected and observed allele frequency in sample.
- The right figure shows changes in MCMC likelihood .

.. image:: _static/PG0390-CNopanel.interpretDEploidFigure.2.png
   :width: 1024px
   :alt: interpretDEploidFigure.2

This panel figure shows all allele frequencies within sample across all 14 chromosomes. Expected and observed WSAF are marked in blue and red respectively.


Example 2
*********

::

    $ ./dEploid -ref data/exampleData/PD0577-C_ref.trim.txt \
    -alt data/exampleData/PD0577-C_alt.trim.txt \
    -plaf data/exampleData/asia-1_PLAF.trim.txt \
    -panel data/exampleData/asia-1_panel.trim.txt \
    -o PD0577-CPanel \
    -k 3

    $ R --slave "--args -ref data/exampleData/PD0577-C_ref.trim.txt
    -alt data/exampleData/PD0577-C_alt.trim.txt
    -plaf data/exampleData/asia-1_PLAF.trim.txt
    -o PD0577-CPanel \
    -dEprefix PD0577-CPanel
    -inbreeding" < utilities/interpretDEploid.r

    $ ./dEploid -ref data/exampleData/PD0577-C_ref.trim.txt \
    -alt data/exampleData/PD0577-C_alt.trim.txt \
    -plaf data/exampleData/asia-1_PLAF.trim.txt \
    -panel data/exampleData/asia-1_panel.trim.txt \
    -o PD0577-CPanel \
    -painting PD0577-CPanel.hap \
    -k 3 -initialP 0.5 0.25 0.25

    $ R --slave "--args -ref data/exampleData/PD0577-C_ref.trim.txt
    -alt data/exampleData/PD0577-C_alt.trim.txt
    -plaf data/exampleData/asia-1_PLAF.trim.txt
    -o PD0577-CPanel \
    -dEprefix PD0577-CPanel
    -inbreeding" < utilities/interpretDEploid.r
