# Script to create PBS scripts for running the BISCUIT benchmarking
mkdir -p slurm

{
    while read sname; do
        sed "s/QQQ/${sname}/g" smallwood_mm10_biscuit > slurm/bisc_${sname}.slurm
        sed "s/QQQ/${sname}/g" smallwood_mm10_bismark > slurm/bism_${sname}.slurm
        sed "s/QQQ/${sname}/g" smallwood_mm10_bwameth > slurm/bwam_${sname}.slurm
        sed "s/QQQ/${sname}/g" smallwood_mm10_gembees > slurm/gemb_${sname}.slurm
        sed "s/QQQ/${sname}/g" smallwood_mm10_bsbohlt > slurm/bolt_${sname}.slurm
    done
} < smallwood_samples
