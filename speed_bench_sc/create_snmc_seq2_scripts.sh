# Script to create PBS scripts for running the BISCUIT benchmarking
mkdir -p slurm

{
    while read sname; do
        sed "s/QQQ/${sname}/g" snmc_seq2_hg38_biscuit > slurm/bisc_${sname}.slurm
        sed "s/QQQ/${sname}/g" snmc_seq2_hg38_bismark > slurm/bism_${sname}.slurm
        sed "s/QQQ/${sname}/g" snmc_seq2_hg38_bwameth > slurm/bwam_${sname}.slurm
        sed "s/QQQ/${sname}/g" snmc_seq2_hg38_gembees > slurm/gemb_${sname}.slurm
        sed "s/QQQ/${sname}/g" snmc_seq2_hg38_bsbohlt > slurm/bolt_${sname}.slurm
    done
} < snmc_seq2_human_samples

{
    while read sname; do
        sed "s/QQQ/${sname}/g" snmc_seq2_mm10_biscuit > slurm/bisc_${sname}.slurm
        sed "s/QQQ/${sname}/g" snmc_seq2_mm10_bismark > slurm/bism_${sname}.slurm
        sed "s/QQQ/${sname}/g" snmc_seq2_mm10_bwameth > slurm/bwam_${sname}.slurm
        sed "s/QQQ/${sname}/g" snmc_seq2_mm10_gembees > slurm/gemb_${sname}.slurm
        sed "s/QQQ/${sname}/g" snmc_seq2_mm10_bsbohlt > slurm/bolt_${sname}.slurm
    done
} < snmc_seq2_mouse_samples
