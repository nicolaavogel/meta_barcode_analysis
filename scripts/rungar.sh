module load samtools

samtools faidx ${1}.fa

/projects/wintherpedersen/apps/gargammel/src/fragSim -n $3 --loc 3.7344 --scale 0.35 ${1}.fa > reads_${2}.fa
/projects/wintherpedersen/apps/gargammel/src/deamSim -matfile dhigh reads_${2}.fa > reads_${2}High.fa
/projects/wintherpedersen/apps/gargammel/src/adptSim -l 140 -artp reads_${2}HighAdapt.fa reads_${2}High.fa
/projects/wintherpedersen/apps/gargammel/art_src_MountRainier_Linux/art_illumina -ss HS25 -amp -na -p -l 140 -c 1 -i reads_${2}HighAdapt.fa -o ${2}High
/projects/wintherpedersen/apps/gargammel/leeHom --ancientdna -fq1 ${2}High1.fq -fq2 ${2}High2.fq -fqo ${2}_${3}

