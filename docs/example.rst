.. _sec-example:

=======
Example
=======

.. todo::
    Full example of working pipeline will be provided with the *Pf3k* pilot paper.


.. .. note::
..     Caveat: need to run the program multiple times, because some models are harder than the others.


.. ::

..     $ ./dEploid -ref labStrains/PG0390_first100ref.txt \
..     -alt labStrains/PG0390_first100alt.txt \
..     -plaf labStrains/labStrains_first100_PLAF.txt \
..     -panel labStrains/lab_first100_Panel.txt \
..     -o tmp1

..     $ ./dEploid -vcf tests/testData/PG0389-C.vcf -plaf tests/testData/labStrains_samples_PLAF.txt -noPanel -o PG0389-CNopanel
..     $ ./dEploid -vcf tests/testData/PG0389-C.vcf -plaf tests/testData/labStrains_samples_PLAF.txt -panel tests/testData/clonalPanel.csv -o PG0389-Cpanel -exportPostProb -seed 2

.. ###Input examples :
.. ```bash
.. ./dEploid -ref labStrains/PG0390_first100ref.txt -alt labStrains/PG0390_first100alt.txt -plaf labStrains/labStrains_first100_PLAF.txt -panel labStrains/lab_first100_Panel.txt -o tmp1
.. ./dEploid -ref labStrains/PG0390_first100ref.txt -alt labStrains/PG0390_first100alt.txt -plaf labStrains/labStrains_first100_PLAF.txt -panel labStrains/lab_first100_Panel.txt -nSample 100 -rate 3
.. ```

