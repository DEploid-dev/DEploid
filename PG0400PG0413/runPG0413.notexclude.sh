#!/bin/bash
#$ -cwd
#$ -V
#$ -P mcvean.prja -q short.qa
#$ -e ErrFiles
#$ -o OutFiles
#$ -N PG0413
#$ -t 1-30

sample="PG0413.C"
root="/well/mcvean/joezhu/pf3k/pf3k_5_1/labStrains/"
add=""

ref=${root}${sample}_ref.txt
alt=${root}${sample}_alt.txt
plaf=${root}Hb3_7g8_plaf.txt
panel=${root}Hb3_7g8.csv
seed=${SGE_TASK_ID}

prefix=${root}${sample}/${sample}_panel_seed${seed}${add}

time pfDeconv -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-panel ${panel} \
-seed ${seed} \
-o ${prefix} \
-burn 0

R --slave "--args ${prefix} ${ref} ${alt} ${plaf}" < ../utilities/makePlots.r

prefix=${root}${sample}/${sample}_panel_seed${seed}miss0${add}

time pfDeconv -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-panel ${panel} \
-seed ${seed} \
-o ${prefix} \
-burn 0 \
-miss 0

R --slave "--args ${prefix} ${ref} ${alt} ${plaf}" < ../utilities/makePlots.r

prefix=${root}${sample}/${sample}_panel_seed${seed}missSmall${add}

time pfDeconv -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-panel ${panel} \
-seed ${seed} \
-o ${prefix} \
-burn 0 \
-miss 0.00000000000000000000000000001

R --slave "--args ${prefix} ${ref} ${alt} ${plaf}" < ../utilities/makePlots.r


prefix=${root}${sample}/${sample}_panel_seed${seed}miss0recom1${add}

time pfDeconv -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-panel ${panel} \
-seed ${seed} \
-o ${prefix} \
-burn 0 \
-miss 0 \
-recomb 1

R --slave "--args ${prefix} ${ref} ${alt} ${plaf}" < ../utilities/makePlots.r

prefix=${root}${sample}/${sample}_noPanel_seed${seed}${add}

time pfDeconv -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-noPanel \
-seed ${seed} \
-o ${prefix} \
-burn 0 \

R --slave "--args ${prefix} ${ref} ${alt} ${plaf}" < ../utilities/makePlots.r


prefix=${root}${sample}/${sample}_panel_seed${seed}bothSmall

time pfDeconv -ref ${ref} \
-alt ${alt} \
-plaf  ${plaf} \
-panel ${panel} \
-seed ${seed} \
-o ${prefix} \
-burn 0 \
-miss 0.001 \
-recomb 0.001

R --slave "--args ${prefix} ${ref} ${alt} ${plaf} " < ../utilities/makePlots.r
