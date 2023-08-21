# Script to create SLURM scripts for running the BISCUIT benchmarking
{
    while read sname; do
        sed "s/QQQ/${sname}/g" slurm_sub_trim_all > ${sname}.slurm
    done
} < sub_trim_all_samples

{
    while read sname; do
        sed "s/QQQ/${sname}/g" slurm_sub_trim_some > ${sname}.slurm
    done
} < sub_trim_some_samples

{
    while read sname; do
        sed "s/QQQ/${sname}/g" slurm_sub_all > ${sname}.slurm
    done
} < sub_all_samples

{
    while read sname; do
        sed "s/QQQ/${sname}/g" slurm_trim_all > ${sname}.slurm
    done
} < trim_all_samples

# Create a bash script to loop through slurm scripts and submit them to the queue
cat > submit_slurm_scripts.sh <<EOF
for FILE in \`ls *.slurm\`; do
    sbatch \${FILE}
    sleep 1.0s
done
EOF
