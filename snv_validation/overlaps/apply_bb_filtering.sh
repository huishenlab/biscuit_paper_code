WORKDIR=/path/to/working/directory
DATA=${WORKDIR}/snv_validation/overlaps/wgbs_overlaps/filtered_vcfs
ASSETS=${WORKDIR}/snv_validation/assets

# Annotate file
bcftools annotate \
    -R ${ASSETS}/hg38.whitelist.sorted.bed.gz \
    -O z \
    -a ${ASSETS}/dbSnp153Common.snv.reformatted.txt.gz \
    -h ${ASSETS}/dbsnp_common.hdr \
    -c CHROM,FROM,TO,TYPE,COMMON_SOME,COMMON_ALL,REF_MIN,ALT_MIN,REF_DBSNP,ALT_DBSNP,REF_ALL,ALT_ALL,RSID,MAX_MAF \
    ${DATA}/removed_lowqual_genotype_chromosome.vcf.gz | \
bcftools view \
    -O z \
    -i'ALT!="N" & ALT!="." & ( (COUNT(GQ>=15)>=1 & COMMON_ALL==1 & MAX_MAF>=0.05) | (COUNT(GQ>=60)>=1) )' \
    - \
> ${DATA}/bb_filtered.vcf.gz
tabix -p vcf ${DATA}/bb_filtered.vcf.gz

