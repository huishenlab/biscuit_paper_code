# dbSNP-related
#-------------------------------------------------------------------------------
# Download file from UCSC
wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb

# Convert to BED file and extract only SNVs
bigBedToBed dbSnp153Common.bb stdout | \
grep "snv" | \
bgzip > dbSnp153Common.snv.bed.gz

# Reformat file for input to bcftools annotate
python parse_dbSNP_bed.py \
    dbSnp153Common.snv.bed.gz | \
grep "snv" | \
bgzip > dbSnp153Common.snv.reformatted.txt.gz
tabix -s 1 -b 2 -e 3 dbSnp153Common.snv.reformatted.txt.gz

# Clean up
rm dbSnp153Common.bb
rm dbSnp153Common.snv.bed.gz
#-------------------------------------------------------------------------------

# Imprinting regions for hg38
#-------------------------------------------------------------------------------
# Retrieve needed files
wget https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg19.imprinting.tsv.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Remove header line
zcat EPIC.hg19.imprinting.tsv.gz | \
tail -n+2 | \
sort -k1,1 -k2,2n > EPIC.hg19.imprinting.tsv

# Lift over from hg19 to hg38
liftOver EPIC.hg19.imprinting.tsv hg19ToHg38.over.chain.gz EPIC.hg38.imprinting.unsorted.tsv EPIC.hg38.imprinting.unmapped

# Sort lifted output
sort -k1,1 -k2,2n EPIC.hg38.imprinting.unsorted.tsv > EPIC.hg38.imprinting.tsv

cat EPIC.hg38.imprinting.tsv | \
cut -f5 | \
sort | \
uniq -c | \
awk '{ print $2 }' | \
sed 's/[:-]/ /g' | \
awk '{ print $1"\t"$2"\t"$3 }' | \
sort -k1,1 -k2,2n | \
bedtools merge -i - -d 100 | \
awk '{ print $1"\t"$2"\t"$3 }' > EPIC.hg19.imprinting.unique_regions.bed

liftOver \
    EPIC.hg19.imprinting.unique_regions.bed \
    hg19ToHg38.over.chain.gz \
    EPIC.hg38.imprinting.unique_regions.unsorted.bed  \
    EPIC.hg38.imprinting.unique_regions.unmapped

sort -k1,1 -k2,2n EPIC.hg38.imprinting.unique_regions.unsorted.bed | \
awk '{
    print $1":"$2"-"$3 > "EPIC.hg38.imprinting.unique_regions.txt"
    print
}' > EPIC.hg38.imprinting.unique_regions.bed

# Clean up
rm EPIC.hg19.imprinting.tsv.gz
rm EPIC.hg38.imprinting.unsorted.tsv
rm EPIC.hg38.imprinting.unmapped
rm EPIC.hg38.imprinting.unique_regions.unsorted.bed
rm EPIC.hg38.imprinting.unique_regions.unmapped
rm EPIC.hg19.imprinting.tsv
rm hg19ToHg38.over.chain.gz
#-------------------------------------------------------------------------------

# Whitelist region
#-------------------------------------------------------------------------------
# Download file
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz

# Create whitelist
WORKDIR=/path/to/working/directory
zcat hg38.blacklist.bed.gz | \
sort -k1,1 -k2,2n | \
bedtools complement -i - -g <( sort -k1,1 ${WORKDIR}/biscuit/hg38/hg38_noContig.fa.fai) | \
gzip > hg38.whitelist.sorted.bed.gz

# Clean up
rm hg38.blacklist.bed.gz
#-------------------------------------------------------------------------------
