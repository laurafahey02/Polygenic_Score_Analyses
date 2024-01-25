#!/bin/bash
# PRSet

set -e -x -o pipefail

main() {
    # print input data objects / links
    echo "Value of base_assoc: '$base_assoc'"
    echo "Value of plink_bed: '$plink_bed'"
    echo "Value of plink_bim: '$plink_bim'"
    echo "Value of plink_fam: '$plink_fam'"

    docker load -i /docker.tar.gz

    # initialize variables
    prsice2_input_options=""

    # download and prep inputs
    base_assoc_prefix="${base_assoc_name/.assoc/}"
    dx download "$base_assoc" -o "${base_assoc_prefix}.assoc"

    # download pheno file if present
    if [[ -n "${pheno_txt}" ]]; then
        pheno_txt_prefix="${pheno_txt_name/.pheno/}"
        dx download "$pheno_txt" -o "${pheno_txt_prefix}.pheno"
        prsice2_input_options="--pheno /home/dnanexus/$pheno_txt_prefix.pheno"
    fi
    # Download samples2remove file if present
    if [ -n "$remove_samples" ]
    then
        dx download "$remove_samples" -o remove_samples
    fi
   # Download gtf file
    if [ -n "$gtf" ]
    then
        dx download "$gtf" -o gtf
    fi
    # Download pathway file
    if [ -n "$msigdb" ]
    then
        dx download "$msigdb" -o msigdb
    # Download filt with snps to extract
    fi
    if [ -n "$extract_snps" ]
    then
        dx download "$extract_snps" -o extract_snps
    fi

    # renaming BED, BIM, FAM files since PRSice2 expects a common prefix
    dx download "$plink_bed" -o target.bed
    dx download "$plink_bim" -o target.bim
    dx download "$plink_fam" -o target.fam

    # adding additional command line arguments
    if [ "$extra_options" != "" ]; then
        prsice2_input_options="$prsice2_input_options $extra_options"
    fi

    # running PRSice2
    docker run -v /home/dnanexus/:/home/dnanexus prsice2:latest \
       Rscript /opt/PRSice.R \
       --prsice ./opt/PRSice_linux \
       --base /home/dnanexus/"${base_assoc_prefix}.assoc" \
       --target /home/dnanexus/target  --binary-target F --thread $(nproc) \
       ${prsice2_input_options} \
       --gtf /home/dnanexus/gtf \
       --msigdb /home/dnanexus/msigdb \
       --extract /home/dnanexus/extract_snps \
       --remove /home/dnanexus/remove_samples --ignore-fid\
       --out "/home/dnanexus/${base_assoc_prefix}"

    # upload outputs
    summary_txt=$(dx upload ${base_assoc_prefix}.summary --brief)
    best_txt=$(dx upload ${base_assoc_prefix}.best  --brief)
    #log_txt=$(dx upload ${base_assoc_prefix}.log  --brief)

    dx-jobutil-add-output summary_txt "${summary_txt}" --class=file
    dx-jobutil-add-output best_txt "${best_txt}" --class=file
   # dx-jobutil-add-output best_txt "${log_txt}" --class=file

}
