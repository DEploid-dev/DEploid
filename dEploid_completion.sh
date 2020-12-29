#/usr/bin/env bash

_dEploid_completions()
{
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts="-h -help \
    -ref -alt -plaf -panel -exclude\
    -vcf -sample -plafFromVcf \
    -o -seed -noPanel\
    -ibd -lasso -best"


    if [[ ${cur} == -* || ${prev} == "dEploid" ]] ; then
        COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
        return 0
    else
        if [[ ${prev} == "-ref" || ${prev} == "-alt" || ${prev} == "-plaf" || ${prev} == "-panel" || ${prev} == "-exclude" ]] ; then
            COMPREPLY=( $(compgen -o plusdirs -f -X '!*.@(txt|gz)' -- ${cur}) )
            for ((i=0; i < ${#COMPREPLY[@]}; i++)); do
                [ -d "${COMPREPLY[$i]}" ] && COMPREPLY[$i]=${COMPREPLY[$i]}/
            done
            return 0
        else
            if [[ ${prev} == "-vcf" ]] ; then
                COMPREPLY=( $(compgen -o plusdirs -f -X '!*.@(vcf|gz)' -- ${cur}) )
                for ((i=0; i < ${#COMPREPLY[@]}; i++)); do
                    [ -d "${COMPREPLY[$i]}" ] && COMPREPLY[$i]=${COMPREPLY[$i]}/
                done
                return 0
            else
                return 0
            fi
        fi
    fi
    return 0
}

complete -o filenames -o nospace -o bashdefault -F _dEploid_completions dEploid
