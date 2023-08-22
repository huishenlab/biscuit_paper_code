################################################################################
#####
##### Code to apply filters to BISCUIT VCF file to try and wittle variants down
##### to those that are likely germline variants
#####
##### Creator: Jacob Morrison
#####
##### Created: Jan 2020
#####
################################################################################

# Usage
function usage
{
    echo "bash apply_qual_genotype_chrom_filters.sh [optional arguments] -v biscuit_pileup.vcf.gz -o output"
    echo -e "\tRequired:"
    echo -e "\t\t-v, --vcf - gzipped output file from biscuit pileup"
    echo -e "\t\t-o, --output - name of output file, don't include file extension"
    echo -e "\tOptional:"
    echo -e "\t\t-f, --filter - Filter LOW QUAL quality filter variants"
    echo -e "\t\t-g, --genotype - Filter genotype="0/0" variants"
    echo -e "\t\t-c, --chromosome - Filter non-canonical chromosome variants"
    echo -e "\t\t-b, --biscuitgatk - Apply BISCUIT and GATK filters"
}

# Base filter
function base_filter
{
    zcat ${FILE} |
    awk '{
        if(substr($1,1,1) == "#") {
            print
        }
        if(substr($1,1,1) != "#" && $5 != "." && $5 != "N") {
            print
        }
    }' - > ${OUTY}
    bgzip ${OUTY}
    tabix -p vcf ${OUTY}.gz
}

# Low quality filter filter
function low_qual
{
    zcat ${FILE} |
    awk '{
        if(substr($1,1,1) == "#") {
            print
        }
        if(substr($1,1,1) != "#" && $5 != "." && $5 != "N" && $7 != "LowQual") {
            print
        }
    }' - > ${OUTY}
    bgzip ${OUTY}
    tabix -p vcf ${OUTY}.gz
}

# Genotype filter
function genotype
{
    zcat ${FILE} |
    awk '{
        if(substr($1,1,1) == "#") {
            print
        }
        if(substr($1,1,1) != "#" && $5 != "." && $5 != "N") {
            n_info = split($9,info,":")
            n_form = split($10,form,":")
            for (i=1; i<n_info; ++i) {
                if (info[i] == "GT") { gt_idx = i }
            }
            if (form[gt_idx] != "0/0") {
                print
            }
        }
    }' - > ${OUTY}
    bgzip ${OUTY}
    tabix -p vcf ${OUTY}.gz
}

# Chromosome filter
function chromosome
{
    zcat ${FILE} |
    awk '{
        if(substr($1,1,1) == "#") {
            print
        }
        if(substr($1,1,1) != "#" && $5 != "." && $5 != "N") {
            if ($1 ~ /^chr[1234567890X]{1,2}$/) {
                print
            }
        }
    }' - > ${OUTY}
    bgzip ${OUTY}
    tabix -p vcf ${OUTY}.gz
}

# Low quality filter + Genotype filter
function low_qual_genotype
{
    zcat ${FILE} |
    awk '{
        if(substr($1,1,1) == "#") {
            print
        }
        if(substr($1,1,1) != "#" && $5 != "." && $5 != "N" && $7 != "LowQual") {
            n_info = split($9,info,":")
            n_form = split($10,form,":")
            for (i=1; i<n_info; ++i) {
                if (info[i] == "GT") { gt_idx = i }
            }
            if (form[gt_idx] != "0/0") {
                print
            }
        }
    }' - > ${OUTY}
    bgzip ${OUTY}
    tabix -p vcf ${OUTY}.gz
}

# Low quality filter + Chromosome filter
function low_qual_chromosome
{
    zcat ${FILE} |
    awk '{
        if(substr($1,1,1) == "#") {
            print
        }
        if(substr($1,1,1) != "#" && $5 != "." && $5 != "N" && $7 != "LowQual") {
            if ($1 ~ /^chr[1234567890X]{1,2}$/) {
                print
            }
        }
    }' - > ${OUTY}
    bgzip ${OUTY}
    tabix -p vcf ${OUTY}.gz
}

# Low quality filter + Genotype + Chromosome filter
function genotype_chromosome
{
    echo ${FILE}
    echo ${OUTY}
    zcat ${FILE} |
    awk '{
        if(substr($1,1,1) == "#") {
            print
        }
        if(substr($1,1,1) != "#" && $5 != "." && $5 != "N") {
            if ($1 ~ /^chr[1234567890X]{1,2}$/) {
                n_info = split($9,info,":")
                n_form = split($10,form,":")
                for (i=1; i<n_info; ++i) {
                    if (info[i] == "GT") { gt_idx = i }
                }
                if (form[gt_idx] != "0/0") {
                    print
                }
            }
        }
    }' - > ${OUTY}
    bgzip ${OUTY}
    tabix -p vcf ${OUTY}.gz
}

# Low quality filter + Genotype + Chromosome filter
function low_qual_genotype_chromosome
{
    zcat ${FILE} |
    awk '{
        if(substr($1,1,1) == "#") {
            print
        }
        if(substr($1,1,1) != "#" && $5 != "." && $5 != "N" && $7 != "LowQual") {
            if ($1 ~ /^chr[1234567890X]{1,2}$/) {
                n_info = split($9,info,":")
                n_form = split($10,form,":")
                for (i=1; i<n_info; ++i) {
                    if (info[i] == "GT") { gt_idx = i }
                }
                if (form[gt_idx] != "0/0") {
                    print
                }
            }
        }
    }' - > ${OUTY}
    bgzip ${OUTY}
    tabix -p vcf ${OUTY}.gz
}

# Low quality filter + Genotype + Chromosome filter
function biscuit_gatk
{
    zcat ${FILE} |
    awk '{
        if(substr($1,1,1) == "#") {
            print
        }
        if(substr($1,1,1) != "#" && $5 != "." && $5 != "N" && $7 != "LowQual") {
            if ($1 ~ /^chr[1234567890X]{1,2}$/ && $6 > 30) {
                n_info = split($9,info,":")
                n_form = split($10,form,":")
                for (i=1; i<n_info; ++i) {
                    if (info[i] == "GT") { gt_idx = i }
                    if (info[i] == "DP") { dp_idx = i }
                }
                if (form[gt_idx] != "0/0" && $6/form[dp_idx] < 2) {
                    print
                }
            }
        }
    }' - > ${OUTY}
    bgzip ${OUTY}
    tabix -p vcf ${OUTY}.gz
}

# Init variables
FILT="false"
GENO="false"
CHRM="false"
GATK="false"
FILE=
OUTY=

# Get command line arguments
if [ "${1}" = "" ]; then
    echo "No command line arguments provided"
    usage
    exit 1
fi

while [ "${1}" != "" ]; do
    case ${1} in
        -b | --biscuitgatk )
            GATK="true"
            ;;
        -f | --filter )
            FILT="true"
            ;;
        -g | --genotype )
            GENO="true"
            ;;
        -c | --chromosome )
            CHRM="true"
            ;;
        -v | --vcf )
            shift
            FILE=${1}
            ;;
        -o | --output )
            shift
            OUTY=${1}.vcf
            ;;
        -h | --help )
            usage
            exit
            ;;
        * )
            echo "Unknown command line argument"
            usage
            exit 1
    esac
    shift
done

# Check VCF file was supplied
if [[ "${FILE}" = ""  || "${OUTY}" = "" ]]; then
    echo "VCF or output file missing"
    usage
    exit 1
fi

# Output information about filtering process
echo "The following variants will be filterd from the VCF file: ${FILE}"
if [ "${FILT}" = "true" ]; then
    echo "Filtering LOW QUAL quality filter variants"
fi
if [ "${GENO}" = "true" ]; then
    echo "Filtering genotype="0/0" variants"
fi
if [ "${CHRM}" = "true" ]; then
    echo "Filtering non-canonical chromosome variants"
fi
if [ "${GATK}" = "true" ]; then
    echo "Applying filters to help match GATK and BISCUIT filters"
fi

if [[ "${FILT}" != "true" && "${GENO}" != "true" && "${CHRM}" != "true" && "${GATK}" != "true" ]]; then
    echo "Only variants with no ALT allele (ALT = '.') or ambiguous ALT alleles (ALT = 'N') will be filtered"
else
    echo "Additionally, variants with no ALT allele (ALT = '.') or ambiguous ALT alleles (ALT = 'N') will be filtered"
fi
echo "Variants passing the filter will be written to: ${OUTY}"

# Apply filters
# Base filter only
if [[ "${FILT}" != "true" && "${GENO}" != "true" && "${CHRM}" != "true" && "${GATK}" != "true" ]]; then
    base_filter
fi

# Single filter
if [[ "${FILT}" = "true" && "${GENO}" != "true" && "${CHRM}" != "true" ]]; then
    low_qual
fi
if [[ "${FILT}" != "true" && "${GENO}" = "true" && "${CHRM}" != "true" ]]; then
    genotype
fi
if [[ "${FILT}" != "true" && "${GENO}" != "true" && "${CHRM}" = "true" ]]; then
    chromosome
fi

# Two filters
if [[ "${FILT}" = "true" && "${GENO}" = "true" && "${CHRM}" != "true" ]]; then
    low_qual_genotype
fi
if [[ "${FILT}" = "true" && "${GENO}" != "true" && "${CHRM}" = "true" ]]; then
    low_qual_chromosome
fi
if [[ "${FILT}" != "true" && "${GENO}" = "true" && "${CHRM}" = "true" ]]; then
    genotype_chromosome
fi

# All filters
if [[ "${FILT}" = "true" && "${GENO}" = "true" && "${CHRM}" = "true" ]]; then
    low_qual_genotype_chromosome
fi

if [[ "${FILT}" != "true" && "${GENO}" != "true" && "${CHRM}" != "true" && "${GATK}" = "true" ]]; then
    biscuit_gatk
fi
