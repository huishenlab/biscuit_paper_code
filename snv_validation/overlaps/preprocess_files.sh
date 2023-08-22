WORKDIR=/path/to/working/directory
DIRLOC=${WORKDIR}/snv_validation
ASSETS=${DIRLOC}/assets

FILE=${DIRLOC}/overlaps/wgbs_overlaps/filtered_vcfs/removed_base_only.vcf.gz
TRUTH=${DIRLOC}/giab/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_indels_removed.vcf.gz

# Merge GIAB and bisucit
bcftools merge --threads 1 --regions chr11:1-22000000 -O z \
    ${FILE} \
    ${TRUTH} \
> merged.hg38.vcf.gz
tabix -p vcf merged.hg38.vcf.gz

# Remove phased haplotypes - logic requires unphased
zcat merged.hg38.vcf.gz | \
perl -ne \
    '$repcol=10; $line=$_; chomp $line; if ($line=~/^\s*\#/) {print $line."\n";} else {@f=split(/\t/,$line); $gt=$f[$repcol]; @ff=split(/:/,$gt); if ($ff[0]=~/^(\d+)\|(\d+)$/) {$ff[0]=($2>$1)?"$1/$2":"$2/$1"}; $f[$repcol]=join(":",@ff); print join("\t",@f)."\n";}' | \
bcftools sort | \
bgzip --threads 1 > merged.rmhap.hg38.vcf.gz
tabix -p vcf merged.rmhap.hg38.vcf.gz

# Intersect with whitelist and dbSNP
bcftools annotate --threads 1 -O z \
    -R ${ASSETS}/hg38.whitelist.sorted.bed.gz \
    -a ${ASSETS}/dbSnp153Common.snv.reformatted.txt.gz \
    -h ${ASSETS}/dbsnp_common.hdr \
    -c CHROM,FROM,TO,TYPE,COMMON_SOME,COMMON_ALL,REF_MIN,ALT_MIN,REF_DBSNP,ALT_DBSNP,REF_ALL,ALT_ALL,RSID,MAX_MAF \
    merged.rmhap.hg38.vcf.gz \
> merged.whitelist.dbsnp.hg38.vcf.gz

# Format for processing and plotting
bcftools query -u \
    --format '%CHROM\t%POS\t%RSID\t%COMMON_ALL\t%MAX_MAF[\t%GT][\t%GQ][\t%AF1][\t%DP]\n' \
    merged.whitelist.dbsnp.hg38.vcf.gz \
> merged.whitelist.dbsnp.hg38.processed.tsv
