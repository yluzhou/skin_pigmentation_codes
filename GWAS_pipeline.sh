# GWAS pipeline
plink --vcf input.vcf.gz --maf 0.01 --geno 0.2 --hwe 1e-5 --recode 12 --output-missing-genotype 0 --transpose --out ${tped_prefix} --allow-extra-chr
emmax-kin-intel64 -v -s -d 10 ${tped_prefix}
emmax-intel64 -v -d 10 -t ${tped_prefix} -p ${pheno_file} -k ${tped_prefix}.aIBS.kinf -c ${cov_file} -o ${GWAS_output}

# Conditional analysis
gcta64 --bfile ${bfile_prefix} --maf 0.01 --cojo-p ${p_value} --cojo-file ${summary_statistics_of_phenotype}.ma --cojo-slct --cojo-actual-geno --out ${cojo_output}


# LocusZoom
bcftools view -r chrx:start-end input.vcf.gz -Oz -o target_region.vcf.gz
tabix -p vcf target_region.vcf.gz  
plink --vcf target_region.vcf.gz --make-bed --out target_data

# calculate LD
plink --bfile target_data \
      --r2 \
      --ld-snp rs123 \       # lead variant or interested variant
      --ld-window-kb 1000 \  
      --ld-window 999999 \    
      --ld-window-r2 0 \     
      --out ld_results
      
locuszoom --metal input.file --ld input.ld --refsnp "chrx:pos" --build hg38 --flank 1000kb    # input.file was made in metal format including two columns (MarkerName and P-value)
