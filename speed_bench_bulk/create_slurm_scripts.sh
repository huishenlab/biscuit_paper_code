# Script to create SLURM scripts for running the BISCUIT benchmarking
mkdir -p slurm

SUBS=( 001M 005M 010M 025M 050M 100M 250M )
SMLL=( 001M 005M 010M 025M 050M )
{
    while read sname; do
        for sub in "${SUBS[@]}"; do
            sed "s/QQQ/${sname}_${sub}/g" template_hg38_biscuit > slurm/bisc_${sname}_${sub}.slurm
            sed "s/QQQ/${sname}_${sub}/g" template_hg38_bwameth > slurm/bwam_${sname}_${sub}.slurm
            sed "s/QQQ/${sname}_${sub}/g" template_hg38_gembees > slurm/gemb_${sname}_${sub}.slurm
            sed "s/QQQ/${sname}_${sub}/g" template_hg38_bsbohlt > slurm/bolt_${sname}_${sub}.slurm

            if [[ ${sname} == "TCGA_COAD_A00R" ]] || [[ ${sname} == "TCGA_LUSC_2600" ]]; then
                sed "s/QQQ/${sname}_${sub}/g; s/IS_SRA=true/IS_SRA=false/g" template_hg38_bismark > slurm/bism_${sname}_${sub}.slurm
            else
                sed "s/QQQ/${sname}_${sub}/g" template_hg38_bismark > slurm/bism_${sname}_${sub}.slurm
            fi
        done
    done
} < human_samples

{
    while read sname; do
        for sub in "${SUBS[@]}"; do
            sed "s/QQQ/${sname}_${sub}/g" template_mm10_biscuit > slurm/bisc_${sname}_${sub}.slurm
            sed "s/QQQ/${sname}_${sub}/g" template_mm10_bismark > slurm/bism_${sname}_${sub}.slurm
            sed "s/QQQ/${sname}_${sub}/g" template_mm10_bwameth > slurm/bwam_${sname}_${sub}.slurm
            sed "s/QQQ/${sname}_${sub}/g" template_mm10_gembees > slurm/gemb_${sname}_${sub}.slurm
            sed "s/QQQ/${sname}_${sub}/g" template_mm10_bsbohlt > slurm/bolt_${sname}_${sub}.slurm
        done
    done
} < mouse_samples

{
    while read sname; do
        if [[ ${sname} == "SRR11614917_SRR11614920" ]]; then
            for sml in "${SMLL[@]}"; do
                sed "s/QQQ/${sname}_${sml}/g" template_z11_biscuit > slurm/bisc_${sname}_${sml}.slurm
                sed "s/QQQ/${sname}_${sml}/g" template_z11_bismark > slurm/bism_${sname}_${sml}.slurm
                sed "s/QQQ/${sname}_${sml}/g" template_z11_bwameth > slurm/bwam_${sname}_${sml}.slurm
                sed "s/QQQ/${sname}_${sml}/g" template_z11_gembees > slurm/gemb_${sname}_${sml}.slurm
                sed "s/QQQ/${sname}_${sml}/g" template_z11_bsbohlt > slurm/bolt_${sname}_${sml}.slurm
            done
        else
            for sub in "${SUBS[@]}"; do
                sed "s/QQQ/${sname}_${sub}/g" template_z11_biscuit > slurm/bisc_${sname}_${sub}.slurm
                sed "s/QQQ/${sname}_${sub}/g" template_z11_bismark > slurm/bism_${sname}_${sub}.slurm
                sed "s/QQQ/${sname}_${sub}/g" template_z11_bwameth > slurm/bwam_${sname}_${sub}.slurm
                sed "s/QQQ/${sname}_${sub}/g" template_z11_gembees > slurm/gemb_${sname}_${sub}.slurm
                sed "s/QQQ/${sname}_${sub}/g" template_z11_bsbohlt > slurm/bolt_${sname}_${sub}.slurm
            done
        fi
    done
} < zebrafish_samples

{
    while read sname; do
        sed "s/QQQ/${sname}_XXXM/g; s/IS_TRUSEQ=false/IS_TRUSEQ=true/g" template_hg38_biscuit > slurm/bisc_${sname}_XXXM.slurm
        sed "s/QQQ/${sname}_XXXM/g; s/IS_TRUSEQ=false/IS_TRUSEQ=true/g" template_hg38_bismark > slurm/bism_${sname}_XXXM.slurm
        sed "s/QQQ/${sname}_XXXM/g; s/IS_TRUSEQ=false/IS_TRUSEQ=true/g" template_hg38_bwameth > slurm/bwam_${sname}_XXXM.slurm
        sed "s/QQQ/${sname}_XXXM/g; s/IS_TRUSEQ=false/IS_TRUSEQ=true/g" template_hg38_gembees > slurm/gemb_${sname}_XXXM.slurm
        sed "s/QQQ/${sname}_XXXM/g; s/IS_TRUSEQ=false/IS_TRUSEQ=true/g" template_hg38_bsbohlt > slurm/bolt_${sname}_XXXM.slurm
    done
} < truseq_samples

# Create a bash script to loop through slurm scripts and submit them to the queue
cat > slurm/submit_slurm_scripts.sh <<EOF
for FILE in \`ls *.slurm\`; do
    sbatch \${FILE}
    sleep 1.0s
done
EOF
