DEploid workflow
================

DEploid-Bestpractice workflow
-----------------------------

Enabled by the flag `-best`, Similar to the first version, the main workflow consist with two steps:

  0. Use `dEploid` on clonal samples, and build a reference panel.
  1. For the mixed samples, infer both proportions and haplotypes at the same time with the reference panel provided.

The following key actions are now hidden internally in step 2:

  0. Use high-quality sites and reference panel to learn the number of strains (k) in malaria mixed infection
  1. Fix k, infer strain proportion and within sample relatedness
  2. Fix k and proportion, use a lasso algorithm to select the effective reference panel.

Most of the parameters (set to default) have now been tuned (see [DEploid-Bestpractices-github-actions](https://github.com/DEploid-dev/DEploid-Bestpractices-github-actions) for details).


DEploid-LASSO workflow (third version)
--------------------------------------

This version aims to improve the haplotype inference with reference panel of any size. It was quickly converted into the Best Practice workflow.


DEploid-IBD workflow (second version)
-------------------------------------

The previous workflow struggle for samples with high within-host relatedness. Hence we introduced new procedures and workflow:

  0. Use `dEploid` on clonal samples, and build a reference panel.
  1. Use the IBD method to infer the proportions without a reference panel.
  2. Tune the haplotype with the given reference panel with fixed strain proportions

![workflow](_static/scheme.png "Pf3k workflow")

Black boxes indicate the key deconvolution steps when our program DEploid is used. Boxes in blue and purple represent the input and output respectively at each step. Steps **Deconvolution with IBD** and **Deconvolution with a reference panel** can be combined by using the flag `-ibd`.

*Caveat: need to run the program multiple times, because some models are harder than the others.*


DEploid workflow (first version)
--------------------------------

Our main workflow consist with two steps:

  0. Use `dEploid` on clonal samples, and build a reference panel.
  1. For the mixed samples, infer both proportions and haplotypes at the same time with the reference panel provided.

*Caveat: need to run the program multiple times, because some models are harder than the others.*
