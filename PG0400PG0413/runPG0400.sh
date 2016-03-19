#!/bin/bash
#$ -cwd
#$ -V
##$ -P mcvean.prja -q short.qa
#$ -e ErrFiles
#$ -o OutFiles
#$ -N PG0400
#$ -t 1-6

sample="PG0400.C"
#root="/well/mcvean/joezhu/pf3k/pf3k_5_1/labStrains/"
root="../labStrains/"


pfdeconv=../pfDeconv
#pfdeconv=../pfDeconv_dbg

#dbgsuffix="dbg"
dbgsuffix=""

ref=${root}${sample}_ref.txt
alt=${root}${sample}_alt.txt
plaf=${root}Hb3_7g8_plaf.txt
panel=${root}Hb3_7g8.csv
seed=${SGE_TASK_ID}

exclude=${sample}.exclude.csv

prefix=${root}${sample}/${sample}_exclude_panelOff_seed${seed}${dbgsuffix}

time ${pfdeconv} -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-noPanel \
-seed ${seed} \
-exclude ${exclude} \
-o ${prefix} \
-burn 0.5 -nSample 1000 > ${prefix}.dbg

R --slave "--args ${prefix} ${ref} ${alt} ${plaf} ${exclude}" < ../utilities/makePlots.r

prefix=${root}${sample}/${sample}_exclude_panel_seed${seed}miss0recom1${dbgsuffix}

time ${pfdeconv} -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-panel ${panel} \
-seed ${seed} \
-exclude ${exclude} \
-o ${prefix} \
-burn 0.5 -nSample 1000 \
-miss 0 \
-recomb 1 > ${prefix}.dbg

R --slave "--args ${prefix} ${ref} ${alt} ${plaf} ${exclude}" < ../utilities/makePlots.r

prefix=${root}${sample}/${sample}_exclude_panel_seed${seed}${dbgsuffix}

time ${pfdeconv} -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-panel ${panel} \
-seed ${seed} \
-exclude ${exclude} \
-o ${prefix} \
-burn 0.5 -nSample 1000 > ${prefix}.dbg

R --slave "--args ${prefix} ${ref} ${alt} ${plaf} ${exclude}" < ../utilities/makePlots.r


prefix=${root}${sample}/${sample}_exclude_panel_seed${seed}bothSmall${dbgsuffix}

time ${pfdeconv} -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-panel ${panel} \
-seed ${seed} \
-exclude ${exclude} \
-o ${prefix} \
-burn 0.5 -nSample 1000 \
-miss 0.01 \
-recomb 0.01 > ${prefix}.dbg

R --slave "--args ${prefix} ${ref} ${alt} ${plaf} ${exclude}" < ../utilities/makePlots.r

prefix=${root}${sample}/${sample}_exclude_panel_seed${seed}bothMedian${dbgsuffix}

time ${pfdeconv} -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-panel ${panel} \
-seed ${seed} \
-exclude ${exclude} \
-o ${prefix} \
-burn 0.5 -nSample 1000 \
-miss 0.1 \
-recomb 0.1 > ${prefix}.dbg

R --slave "--args ${prefix} ${ref} ${alt} ${plaf} ${exclude}" < ../utilities/makePlots.r



