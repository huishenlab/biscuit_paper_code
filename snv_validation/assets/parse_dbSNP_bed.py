# Parse a BED file generated from the UCSC dbSNP bigBed file
#
# Code is based on:
#     https://github.com/ekushele/methylseq/blob/dev/bin/processUcscDbsnp.pl
#
# Example of preparing inputs:
#
# wget http://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153Common.bb
#
# bigBedToBed dbSnp153Common.bb stdout | \
# grep "snv" | \
# bgzip > dbSnp153Common.bed.gz
# tabix -p bed dbSnp153Common.bed.gz
#
import gzip
import sys

def uniq_vals(vals):
    """Pick off unique values from comma-separated list.

    Inputs -
        vals - comma-separated list of values
    Returns -
        list
    """
    vals_list = vals.split(',')

    out = {}
    for v in vals_list:
        # Oftentimes there is an extra comma at the end of the list, so ignore blank entries
        if v == '':
            continue

        try:
            out[v] += 1
        except KeyError:
            out[v] = 1

    return list(out.keys())

def max_val(vals):
    """Find the maximum value from a list of strings.

    Inputs -
        vals - list of strings that hold float values
    Returns -
        float
    """
    return max( [float(v) for v in vals] )

def create_header():
    # chr11 87268   87268   snv 1   1   T   C   T   C   T   C   rs574324672 0.227436
    out = [
        '##INFO=<ID=CHROM,Number=".",Type=String,Description="chromosome">',
        '##INFO=<ID=FROM,Number=".",Type=Integer,Description="1-based inclusive start position">',
        '##INFO=<ID=TO,Number=".",Type=Integer,Description="1-based inclusive end position">',
        '##INFO=<ID=TYPE,Number=".",Type=String,Description="Type of variant">',
        '##INFO=<ID=COMMON_SOME,Number=".",Type=Integer,Description="commonSome dbSNP allele">',
        '##INFO=<ID=COMMON_ALL,Number=".",Type=Integer,Description="commonAll dbSNP allele">',
        '##INFO=<ID=REF_MIN,Number=".",Type=String,Description="lower number of REFs between REF_DBSNP and REF_ALL">',
        '##INFO=<ID=ALT_MIN,Number=".",Type=String,Description="lower number of ALTs between ALT_DBSNP and ALT_ALL">',
        '##INFO=<ID=REF_DBSNP,Number=".",Type=String,Description="dbSNP REFs">',
        '##INFO=<ID=ALT_DBSNP,Number=".",Type=String,Description="dbSNP ALTs">',
        '##INFO=<ID=REF_ALL,Number=".",Type=String,Description="All REFs">',
        '##INFO=<ID=ALT_ALL,Number=".",Type=String,Description="All ALTs">',
        '##INFO=<ID=RSID,Number=".",Type=String,Description="dbSNP rs ID">',
        '##INFO=<ID=MAX_MAF,Number=".",Type=Float,Description="Largest MAF value">',
    ]

    with open('dbsnp_common.hdr', 'w') as f:
        f.write('\n'.join(out))

    return None

fin = sys.argv[1]

with gzip.open(fin, 'rt') as f:
    for line in f.readlines():
        l = line.strip().split('\t')

        # Maximum minor allele frequency (MAF)
        all_mafs = uniq_vals(l[9])
        max_maf  = max_val(all_mafs)

        # Get common alleles
        all_refs = uniq_vals(l[10])
        all_alts = uniq_vals(l[11])

        # dbSNP variants
        dbsnp_ref = uniq_vals(l[4])
        dbsnp_alt = uniq_vals(l[6])

        used_dbsnp = (len(dbsnp_alt) < len(all_alts) or len(dbsnp_ref) < len(all_refs))
        min_ref = dbsnp_ref if used_dbsnp else all_refs
        min_alt = dbsnp_alt if used_dbsnp else all_alts

        # Notes
        common_some = '1' if "commonSome" in l[14] else '0'
        common_all  = '1' if "commonAll"  in l[14] else '0'

        # Output
        out = [
            l[0]               , # chromosome
            int(l[1])+1             , # 1-based start
            l[2]               , # 1-based end
            l[13]              , # type
            common_some        , # 1 if commonSome in notes column, else 0
            common_all         , # 1 if commonAll in notes column, else 0
            ','.join(min_ref)  , # lower number of REFs between dbsnp_ref and all_refs
            ','.join(min_alt)  , # lower number of ALTs between dbsnp_alt and all_alts
            ','.join(dbsnp_ref), # list of dbsnp REFs
            ','.join(dbsnp_alt), # list of dbsnp ALTs
            ','.join(all_refs) , # list of all REFs
            ','.join(all_alts) , # list of all ALTs
            l[3]               , # dbSNP reference ID
            max_maf              # Largest MAF value
        ]
        out = [ str(i) for i in out ]

        print('\t'.join(out))

create_header()
