{
    read
    while IFS=$'\t' read -r FQ1 FQ1_MD5 FQ2 FQ2_MD5 NAME; do
        wget --quiet ${FQ1}
        fq1_md5_download=$(md5sum "$(basename -- "${FQ1}")" | cut -d " " -f1)
        fq1_md5_results="Input: ${FQ1_MD5}\nFile: $fq1_md5_download"
        if [[ ${FQ1_MD5} == $fq1_md5_download ]]; then
            echo -e "\n\e[92mSUCCESS\e[39m\n$fq1_md5_results"
        else
            echo -e "\n\e[91mFAILURE\e[39m\n$fq1_md5_results"
        fi

        wget --quiet ${FQ2}
        fq2_md5_download=$(md5sum "$(basename -- "${FQ2}")" | cut -d " " -f1)
        fq2_md5_results="Input: ${FQ2_MD5}\nFile: $fq2_md5_download"
        if [[ ${FQ2_MD5} == $fq2_md5_download ]]; then
            echo -e "\n\e[92mSUCCESS\e[39m\n$fq2_md5_results"
        else
            echo -e "\n\e[91mFAILURE\e[39m\n$fq2_md5_results"
        fi
    done
} < sequence.index.NA12878_Illumina300X_wgs_09252015
