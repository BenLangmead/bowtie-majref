# $1 is the prefix, the file should look like ${1}.vcf.gz
bgzip -cd ${1}.vcf.gz > ${1}_v42.vcf
sed -i 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/' ${1}_v42.vcf
bgzip ${1}_v42.vcf
