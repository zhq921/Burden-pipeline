#!/bin/bash
#SBATCH --mem 80g
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH -p global

# activate conda environment
source ~/.bashrc                   # this cmd will activate conda base env.
eval "$(conda shell.bash hook)"    # this will allow runing conda inside shell scripts.
conda activate burden

# configuration
project_name="CVM/MRKH"
echo "${project_name} project starting..."
basename=./${project_name}/vcf/${project_name}
input_vcf=CVM/MRKH.vcf.gz

# make working directory
cd /mnt/work/research/zhaohq/burden
mkdir -p ./${project_name}
mkdir -p ./${project_name}/vcf
mkdir -p ./${project_name}/samInfo
mkdir -p ./${project_name}/burdenTest/qualified_variant
mkdir -p ./${project_name}/burdenTest/ACAT
mkdir -p ./${project_name}/log
mkdir -p ./${project_name}/gwas

################## Step 1 Preprocessing ##################
in_1_1=${input_vcf}
out_1_1=${basename}.split.gz
# Step 1.1 split multiallelic sites into biallelic records
bcftools norm -m -both -O z \
--threads 10 \
-o ${out_1_1} \
${in_1_1} \
> ./${project_name}/log/1.1-Split_multiallelic.log 2>&1

##################### Step 2 variant/genotype filtration (hard filter 1) #####################
in_2_1=${basename}.split.gz
out_2_1=${basename}.split.hard
in_2_2=${basename}.split.hard.recode.vcf
out_2_2=${basename}.split.hard.AB

# Step 2.1 hard filter
vcftools --gzvcf ${in_2_1} \
--max-missing 0.9 \
--minDP 10 --minGQ 20 --hwe 0.000001 \
--recode --recode-INFO-all \
--out ${out_2_1} \
> ./${project_name}/log/2.1-hard_filter.log 2>&1

# Step 2.2 AB filter
Rscript ./code/AB_filter.R \
-v ${in_2_2} \
-o ${out_2_2} -i 0.25 \
> ./${project_name}/log/2.2-AB_filter.log 2>&1


################# Step 3 GATK VQSR + hard filter2: QD/SOR ##################
in_3_1=${basename}.split.hard.AB
out_3_1=${basename}.split.hard.AB.VQSR.snp
in_3_2=${basename}.split.hard.AB
out_3_2=${basename}.split.hard.AB.VQSR.snp.recal
in_3_3=${basename}.split.hard.AB.VQSR.snp.recal
out_3_3=${basename}.split.hard.AB.VQSR.snp.indel
in_3_4=${basename}.split.hard.AB.VQSR.snp.recal
out_3_4=${basename}.split.hard.AB.VQSR.snp.indel.recal
in_3_5=${basename}.split.hard.AB.VQSR.snp.indel.recal
out_3_5=${basename}.split.hard.AB.VQSR
in_3_6=${basename}.split.hard.AB.VQSR
out_3_6=${basename}.split.hard.AB.VQSR.QDSOR

# Step 3.1 VQSR 1: snp recalibraion
gatk VariantRecalibrator \
--reference /mnt/share/database/hg19/hg19.fa \
--variant ${in_3_1} \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 /mnt/share/database/hg19/hapmap_3.3.hg19.sites.vcf \
--resource:omni,known=false,training=true,truth=true,prior=12.0 /mnt/share/database/hg19/1000G_omni2.5.hg19.sites.vcf \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 /mnt/share/database/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/share/database/hg19/dbsnp_138.hg19.vcf \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum \
-an FS -an SOR -mode SNP -tranche 100.0 -tranche 99.9 \
-tranche 99.0 -tranche 95.0 -tranche 90.0 \
--output ${out_3_1}.recalFile \
--tranches-file ${out_3_1}.tranches \
--rscript-file ${out_3_1}.plots.R \
> ./${project_name}/log/3.1-VQSR_snp_recal.log 2>&1

# Step 3.2 VQSR 2: snp apply
gatk ApplyVQSR \
--reference /mnt/share/database/hg19/hg19.fa \
--variant ${in_3_2} \
--truth-sensitivity-filter-level 99 \
--tranches-file ${out_3_1}.tranches \
--recal-file ${out_3_1}.recalFile \
-mode SNP -O ${out_3_2} \
> ./${project_name}/log/3.2-VQSR_snp_apply.log 2>&1

# Step 3.3 VQSR 3: indel recalibration
gatk VariantRecalibrator \
--reference /mnt/share/database/hg19/hg19.fa \
--variant ${in_3_3} \
--resource:mills,known=true,training=true,truth=true,prior=12.0 /mnt/share/database/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/share/database/hg19/dbsnp_138.hg19.vcf \
-an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
-mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 \
-tranche 90.0 \
--output ${out_3_3}.recalFile \
--max-gaussians 4 \
--tranches-file ${out_3_3}.tranches \
--rscript-file ${out_3_3}.plots.R \
> ./${project_name}/log/3.3-VQSR_indel_recal.log 2>&1

# Step 3.4 VQSR 4: indel apply
gatk ApplyVQSR \
--reference /mnt/share/database/hg19/hg19.fa \
--variant ${in_3_4} \
--truth-sensitivity-filter-level 99 \
--tranches-file ${out_3_3}.tranches \
--recal-file ${out_3_3}.recalFile \
-mode INDEL -O ${out_3_4} \
> ./${project_name}/log/3.4-VQSR_indel_apply.log 2>&1

# Step 3.5 VQSR 5: filter PASS
gatk SelectVariants \
-R /mnt/share/database/hg19/hg19.fa \
--variant ${in_3_5} \
--exclude-filtered --preserve-alleles --output ${out_3_5} \
> ./${project_name}/log/3.5-VQSR_filter_PASS.log 2>&1

# Step 3.6 hard filter 2: QD/SOR
Rscript ./code/QDSOR_filter.R \
-v ${in_3_6} \
-o ${out_3_6} \
> ./${project_name}/log/3.6-QDSOR_filter.log 2>&1

mv ${basename}.split.hard.AB.VQSR.QDSOR ${basename}.split.hard.AB.VQSR.QDSOR.runid

# # get runid2samid
# cd /Users/hqzhao/Desktop/burden/cs_bosheng_1006/saminfo
# cut -f 1,3 sam_info.tsv \
# | sed -e 's/.hc.g.vcf.gz$//g' \
# | awk 'BEGIN{OFS=" "}{print $2,$1}' \
# > runid2samid.txt

bcftools reheader --samples ./${project_name}/samInfo/runid2samid.txt -o \
${basename}.split.hard.AB.VQSR.QDSOR \
${basename}.split.hard.AB.VQSR.QDSOR.runid \
> ./${project_name}/log/3.7-reheader.log 2>&1


##################### Step 4 sample filtration #####################
in_4_1=${basename}.split.hard.AB.VQSR.QDSOR
in_4_2=${basename}.split.hard.AB.VQSR.QDSOR
in_4_3=${basename}.split.hard.AB.VQSR.QDSOR
in_4_4=${basename}.split.hard.AB.VQSR.QDSOR
in_4_5_1=${basename}.split.hard.AB.VQSR.QDSOR
out_4_5_1=${basename}.split.hard.AB.VQSR.QDSOR.nochr
in_4_5_2=${basename}.split.hard.AB.VQSR.QDSOR.nochr
out_4_5_2=${basename}.split.hard.AB.VQSR.QDSOR.nochr.setid
in_4_5_2_1=${basename}.split.hard.AB.VQSR.QDSOR.nochr.setid
out_4_5_2_1=${basename}.split.hard.AB.VQSR.QDSOR.nochr.setid.rmdup
in_4_5_3=${basename}.split.hard.AB.VQSR.QDSOR.nochr.setid.rmdup
out_4_5_3=${basename}
in_4_5_4=${basename}
out_4_5_4=${basename}.binary
out_4_5_4_2=${basename}.binary.clean
in_4_5_5=${basename}.binary
out_4_5_5=./${project_name}/samInfo/${project_name}_pca
out_4_5_5_2=./${project_name}/samInfo/${project_name}_clean_pca
in_4_5_6=${basename}.binary
out_4_5_6=${basename}.1000g.merge
in_4_5_7=${basename}.1000g.merge
out_4_5_7=${basename}.1000g.merge.filter
in_4_5_8=${basename}.1000g.merge.filter
out_4_5_8=./${project_name}/samInfo/${project_name}_1000g_pca
in_4_6=${basename}.split.hard.AB.VQSR.QDSOR
out_4_6=${basename}.split.hard.AB.VQSR.QDSOR.samFlt
out_4_7_1=./${project_name}/samInfo/${project_name}_pca_after_samFlt
out_4_7_2=./${project_name}/samInfo/${project_name}_1000g_pca_after_samFlt
out_4_7_3=./${project_name}/samInfo/${project_name}_clean_pca_after_samFlt
in_4_8=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.recode.vcf
out_4_8=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew

# Step 4.1 missing rate
vcftools --vcf \
${in_4_1} \
--missing-indv \
--out ./${project_name}/samInfo/${project_name} \
> ./${project_name}/log/4.1-sam_miss.log 2>&1

# Step 4.2 depth
vcftools --vcf \
${in_4_2} \
--depth \
--out ./${project_name}/samInfo/${project_name} \
> ./${project_name}/log/4.2-sam_depth.log 2>&1

# Step 4.3 het
vcftools --vcf \
${in_4_3} \
--het \
--out ./${project_name}/samInfo/${project_name} \
> ./${project_name}/log/4.3-sam_het.log 2>&1

# Step 4.4 relatedness
vcftools --vcf ${in_4_4} \
--maf 0.01 \
--out ${in_4_4}.maf1 \
--recode --recode-INFO-all \
> ./${project_name}/log/4.4.1-extract_maf1.log 2>&1
## relatedness
vcftools --vcf \
${in_4_4}.maf1.recode.vcf \
--relatedness2 \
--out ./${project_name}/samInfo/${project_name} \
> ./${project_name}/log/4.4.2-sam_relatedness.log 2>&1

# Step 4.5 PCA filter
# 4.5.1 remove chr (PCA-specific)
awk '{ 
        if($0 !~ /^#/) 
            {gsub(/^chr/,""); print}
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            {gsub(/chr/,""); print}
        else print $0 
      }' ${in_4_5_1} \
       > ${out_4_5_1} \
       2> ./${project_name}/log/4.5.1-rmchr.log

## add chr
## awk '{ 
##         if($0 !~ /^#/) 
##             print "chr"$0;
##         else if(match($0,/(##contig=<ID=)(.*)/,m))
##             print m[1]"chr"m[2];
##         else print $0 
##       }' no_chr.vcf > with_chr.vcf

# 4.5.2 annotate set ID
bcftools annotate --remove ID \
--set-id +'%CHROM:%POS:%REF:%ALT' -O v \
-o ${out_4_5_2} \
${in_4_5_2} \
> ./${project_name}/log/4.5.2-SetId.log 2>&1
## 4.5.2.1 rm dup 
bcftools norm -O v --rm-dup all -o ${out_4_5_2_1} ${in_4_5_2_1} \
> ./${project_name}/log/4.5.2.1-rm_dup_alleles.log 2>&1
## 4.5.3 convert VCF into Plink readable format (tfam,tped)
vcftools --vcf ${in_4_5_3} \
--plink-tped \
--out ${out_4_5_3} \
> ./${project_name}/log/4.5.3-vcf2plink.log 2>&1

## 4.5.4 Plink binary format (fam,bed,bim)
plink --tfile ${in_4_5_4} \
--allow-no-sex --make-bed --noweb \
--out ${out_4_5_4} \
> ./${project_name}/log/4.5.4-plink_binary.log 2>&1
## 4.5.4.2 clean data for PCA
plink --bfile ${in_4_5_5} \
--maf 0.01 --indep-pairwise 50 5 0.2 \
--out ${out_4_5_4_2} \
> ./${project_name}/log/4.5.4.2-plink_binary_clean.log 2>&1

# 4.5.5 PCA
plink --bfile ${in_4_5_5} \
--pca 10 --out ${out_4_5_5} \
> ./${project_name}/log/4.5.5-pca.log 2>&1

# 4.5.5.2 PCA
plink --bfile ${in_4_5_5} \
--pca 10 \
--extract ${out_4_5_4_2}.prune.in \
--out ${out_4_5_5_2} \
> ./${project_name}/log/4.5.5.2-clean_pca.log 2>&1


##4.5.6 Merge 1000 genome and cohort study
plink \
--bfile ${in_4_5_6} \
--bmerge ./source/data1000G/Merge.bed ./source/data1000G/Merge.bim ./source/data1000G/Merge.fam \
--make-bed --out ${out_4_5_6} \
> ./${project_name}/log/4.5.6-merge_1000g.log 2>&1

##4.5.7 filter
plink \
--bfile ${in_4_5_7} \
--hwe 0.000001 --geno 0.05 --maf 0.01 \
--make-bed --out ${out_4_5_7} \
> ./${project_name}/log/4.5.7-merge_1000g.filter.log 2>&1

##4.5.8 perform PCA with 1000G
plink --bfile ${in_4_5_8} \
--pca 10 --out ${out_4_5_8} \
> ./${project_name}/log/4.5.8-merge_1000g.pca.log 2>&1

## Step 4.6 missing rate>0.2
awk 'BEGIN{OFS="\t"}{if(NR>1){if($5>0.2) {print $1,"high missing"}}}' \
./${project_name}/samInfo/${project_name}.imiss \
> ./${project_name}/samInfo/sam2rm_info
## mean depth
depth=20
awk 'BEGIN{OFS="\t"}{if(NR>1){if($3<'"$depth"') {print $1,"low depth"}}}' \
./${project_name}/samInfo/${project_name}.idepth \
>> ./${project_name}/samInfo/sam2rm_info
## relatedness
awk '{if(NR>1){if($7>0.354) print}}' \
./${project_name}/samInfo/${project_name}.relatedness2 \
| awk '{if($7!=0.5) print}' > ./${project_name}/samInfo/${project_name}.related.sam
Rscript ./code/get_related_sam.R -i ./${project_name}/samInfo/${project_name}.related.sam \
-o ./${project_name}/samInfo/sam2rm_info \
> ./${project_name}/log/4.6.1-get_related_sam.log 2>&1
## hetertozygosity
Rscript ./code/get_het_sam_ver2.R -i ./${project_name}/samInfo/${project_name}.het \
-o ./${project_name}/samInfo/sam2rm_info \
> ./${project_name}/log/4.6.2-get_het_sam.log 2>&1
## PCA
Rscript ./code/get_PCA_outlier.R -i ${out_4_5_5}.eigenvec \
-o ./${project_name}/samInfo/sam2rm_info \
> ./${project_name}/log/4.6.3-get_PCA_outlier.log 2>&1

cut -f 1 ./${project_name}/samInfo/sam2rm_info | sort | uniq > ./${project_name}/samInfo/sam2rm
vcftools --vcf ${in_4_6} \
--remove ./${project_name}/samInfo/sam2rm --recode --recode-INFO-all \
--out ${out_4_6} \
> ./${project_name}/log/4.6.4-sample_filter.log 2>&1

# Step 4.7 rerun PCA after sample removal
paste ./${project_name}/samInfo/sam2rm ./${project_name}/samInfo/sam2rm > ./${project_name}/samInfo/sam2rm_plink
# 4.7.1 PCA
plink --bfile ${in_4_5_5} \
--pca 10 --remove ./${project_name}/samInfo/sam2rm_plink \
--out ${out_4_7_1} \
> ./${project_name}/log/4.7.1-pca_after_samFlt.log 2>&1
# 4.7.2 perform PCA with 1000G
plink --bfile ${in_4_5_8} \
--pca 10 --remove ./${project_name}/samInfo/sam2rm_plink \
--out ${out_4_7_2} \
> ./${project_name}/log/4.7.2-merge_1000g.pca_after_samFlt.log 2>&1
# 4.7.3 clean data PCA
plink --bfile ${in_4_5_5} \
--pca 10 --remove ./${project_name}/samInfo/sam2rm_plink \
--extract ${out_4_5_4_2}.prune.in \
--out ${out_4_7_3} \
> ./${project_name}/log/4.7.3-clean_pca_after_samFlt.log 2>&1

# Step 4.8 update AC AF AN
vcffixup ${in_4_8} \
1> ${out_4_8} \
2> ./${project_name}/log/4.8-vcf_renew.log

############################# Step 5 VEP ##################
in_5=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew
out_5=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.vcf

# vep annotation
vep \
-i ${in_5} \
--offline --cache \
--dir_cache /mnt/work/database/vep/vep-104/caches \
--fasta /mnt/share/database/vep/vep-ref/hg19/hg19_vep.fa \
--custom /mnt/share/database/gnomad/2.1.1/hg19/gnomad.exomes.r2.1.1.sites.vcf.bgz,gnomADg,vcf,exact,0,AC,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
--plugin SpliceAI,snv=/mnt/work/database/spliceAI/spliceai_scores.raw.snv.hg19.vcf.gz,indel=/mnt/work/database/spliceAI/spliceai_scores.raw.indel.hg19.vcf.gz \
--plugin dbNSFP,/mnt/work/database/dbNSFP/dbNSFP4.3a/dbNSFP4.3a_grch37.gz,ALL \
--plugin CADD,/mnt/work/database/vep/plugins/CADD/v1.6_GRCh37/whole_genome_SNVs.tsv.gz,/mnt/work/database/vep/plugins/CADD/v1.6_GRCh37/InDels.tsv.gz \
--plugin LoF,loftee_path:/mnt/work/database/vep/plugins/loftee/loftee/,human_ancestor_fa:/mnt/work/database/vep/plugins/loftee/human_ancestor.fa.gz,conservation_file:/mnt/work/database/vep/plugins/loftee/phylocsf_gerp.sql,gerp_file:/mnt/work/database/vep/plugins/loftee/GERP_scores.final.sorted.txt.gz \
--force_overwrite \
--sift b \
--polyphen b \
--regulatory \
--numbers \
--ccds \
--hgvs \
--symbol \
--xref_refseq \
--buffer_size 10000 \
--canonical \
--protein \
--biotype \
--check_existing \
--exclude_null_alleles \
--af \
--af_1kg \
--af_esp \
--af_gnomad \
--max_af \
--vcf \
--fork 32 \
--output_file ${out_5} \
> ./${project_name}/log/5-vep.log 2>&1

####### repeat region removal #######
in_6_1=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.vcf
out_6_1=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.snpOnly
vcftools --vcf ${in_6_1} \
--mac 1 \
--remove-indels --recode \
--recode-INFO-all \
--out ${out_6_1} \
> ./${project_name}/log/6.1-SNP_only.log 2>&1

#### indel only
in_6_1=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.vcf
out_6_3=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.indelOnly
vcftools --vcf ${in_6_1} \
--mac 1 \
--keep-only-indels --recode \
--recode-INFO-all \
--out ${out_6_3} \
> ./${project_name}/log/6.3-indel_only.log 2>&1

#### indel repeat region removal
in_6_2=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.indelOnly.recode.vcf
vcftools --vcf ${in_6_2} \
--exclude-bed /mnt/work/research/zhaohq/burden/source/repeat/rptmsk_merged_0430.bed \
--recode --recode-INFO-all \
--out ${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.indelOnly.rmRepeat \
> ./${project_name}/log/6.4-indel_repeat_removal.log 2>&1

### compress into bgzip
bgzip ${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.indelOnly.rmRepeat.recode.vcf \
2> ./${project_name}/log/6.5.1-compress_indelOnly.log
bgzip ${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.snpOnly.recode.vcf \
2> ./${project_name}/log/6.5.2-compress_snpOnly.log
tabix -p vcf \
${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.indelOnly.rmRepeat.recode.vcf.gz
tabix -p vcf \
${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.snpOnly.recode.vcf.gz

# ### combine a SNP VCF and an indel VCF into one
bcftools concat \
-a -Ov \
-o ${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.indelRmrptSNPraw.vcf \
${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.indelOnly.rmRepeat.recode.vcf.gz \
${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.snpOnly.recode.vcf.gz \
> ./${project_name}/log/6.5.3-merge_indelRmrpt_SNPraw.log 2>&1

Rscript /mnt/work/research/zhaohq/burden/cs_bosheng_1006/code/get_sam_class.R > ./${project_name}/log/6.0.1-get_samClass_afFlt.log 2>&1
######################### Step 6 qualifying variant selection ########################
## update internal AF (3e-04) (AC<=3): (> 3 / (sample_num * 2)) & (< 4 / (sample_num * 2))
## 0.01 0.001 AF
## indelRmrptSNPraw
in_6=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.indelRmrptSNPraw.vcf
## qualifying variant selection
parallel --xapply --header : Rscript ./code/qualifying_variant.R \
-v ${in_6} -o ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw -e {external} -i {internal} -t refseq \
> ./${project_name}/log/6.5-qualify_variant_indelRmrptSNPraw.log 2>&1 ::: external 0.0001 0.001 0.01 ::: internal 3e-04 0.001 0.01
###subgroup
parallel --xapply --header : Rscript ./code/qualifying_subgroup.R \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw -e {external} -i {internal} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} -l ./${project_name}/samInfo/sam_class_afFlt.txt \
> ./${project_name}/log/6.6-qualify_variant_indelRmrptSNPraw_subgroup.log \
2>&1 ::: external 1e-04 0.001 0.01 ::: internal 3e-04 0.001 0.01 ::: subgroupCol 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7

########################## Step 7 burden test ########################
## 7.1 fisher + ACAT + permutation 
## subgroup indelRmrptSNPraw
parallel --xapply -j1 --header : Rscript ./code/ACAT_subgroup.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.1-burden_subgroup_1.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 1 1 1

parallel --xapply -j1 --header : Rscript ./code/ACAT_subgroup.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.1-burden_subgroup_2.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 2 2 2

parallel --xapply -j1 --header : Rscript ./code/ACAT_subgroup.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.1-burden_subgroup_3.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 3 3 3

parallel --xapply -j1 --header : Rscript ./code/ACAT_subgroup.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.1-burden_subgroup_4.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 4 4 4

parallel --xapply -j1 --header : Rscript ./code/ACAT_subgroup.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.1-burden_subgroup_5.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 5 5 5

########################## Step 8 GWAS #######################
##### GWAS for cs sample
##### scoliosis_2021Sep/code/get_sam_class.R to get cs sample information
in_8_1_0=${basename}.split.hard.AB.VQSR.QDSOR.samFlt.renew.vep.indelRmrptSNPraw.vcf
out_8_1_0=./${project_name}/gwas/${project_name}.cs.vcf
in_8_1_1=./${project_name}/gwas/${project_name}.cs.vcf
out_8_1_1=./${project_name}/gwas/${project_name}.gwas.vcf
in_8_1_2=./${project_name}/gwas/${project_name}.gwas.vcf
out_8_1_2=./${project_name}/gwas/${project_name}.gwas.gnomadAFfilter.vcf
in_8_2=./${project_name}/gwas/${project_name}.gwas.gnomadAFfilter.vcf
out_8_2=./${project_name}/gwas/${project_name}.gwas
in_8_3=./${project_name}/gwas/${project_name}.gwas
out_8_3=./${project_name}/gwas/${project_name}.gwas.binary
in_8_4_1=./${project_name}/gwas/${project_name}.gwas.binary
out_8_4_1=./${project_name}/gwas/${project_name}.assoc
out_8_4_2=./${project_name}/gwas/${project_name}.logistic
out_8_4_3=./${project_name}/gwas/${project_name}.selfPCA.logistic
out_8_4_4=./${project_name}/gwas/${project_name}.1000PCA.logistic

# # 8.1 data preparation(setid & getalt)
# # bcftools view -S ./${project_name}/gwas/gwasSam.txt ${in_8_1_0} \
# 1>${out_8_1_0} 2>./${project_name}/log/8.1.0-subset_sample.log
# ## 8.1.1 rm dup
# bcftools norm -O v --rm-dup all -o ${out_8_1_1} ${in_8_1_1}\
# > ./${project_name}/log/8.1.1-rm_dup_alleles.log 2>&1

## 8.1.2 filter gnomad AF 0.01
Rscript ./code/gnomadAF_filter4GWAS.R \
-v ${in_8_1_2} \
-o ${out_8_1_2} -e 0.01 \
> ./${project_name}/log/8.1.2-gnomadAFfilter4GWAS.log 2>&1
## get alternative allele
awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}' ${out_8_1_2} \
1> ./${project_name}/gwas/alt_alleles \
2> ./${project_name}/log/8.1.3-get_alt_alleles.log
## 8.2 convert VCF into Plink readable format (tfam,tped)
vcftools --vcf ${in_8_2} \
--plink-tped \
--out ${out_8_2} \
> ./${project_name}/log/8.2-gwas_vcf2plink.log 2>&1
## 8.3 Plink binary format (fam,bed,bim)
plink --tfile ${in_8_3} \
--allow-no-sex --make-bed --noweb --maf 0.01 \
--snps-only --out ${out_8_3} \
> ./${project_name}/log/8.3-plink_binary.log 2>&1
# 8.4 gwas runing
## sample information
awk 'BEGIN{FS="\t";OFS="\t";}{print $1,$1,$2}' \
./${project_name}/gwas/gwasSam_class.txt \
> ./${project_name}/gwas/sam.pheno
## 8.4.1 Run a simple association analysis
plink --bfile ${in_8_4_1} \
--make-pheno ./${project_name}/gwas/sam.pheno "case" \
--assoc perm --reference-allele ./${project_name}/gwas/alt_alleles \
--allow-no-sex --adjust --noweb \
--out ${out_8_4_1} \
> ./${project_name}/log/8.4.1-assoc_test_chisq.log 2>&1
## 8.4.2 Run a logistic/lasso regression
## logistic without covariate
plink --bfile ${in_8_4_1} \
--make-pheno ./${project_name}/gwas/sam.pheno "case" \
--ci 0.95 \
--logistic perm \
--reference-allele ./${project_name}/gwas/alt_alleles \
--allow-no-sex --adjust --noweb \
--out ${out_8_4_2} \
> ./${project_name}/log/8.4.2-assoc_test_logistic.log 2>&1

### sample information # self-PCA
echo -e "FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10" > ./${project_name}/gwas/sam_selfPCA.cov
awk 'BEGIN{FS=" ";OFS="\t";}{print $2,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
./${project_name}/samInfo/${project_name}_clean_pca_after_samFlt.eigenvec \
>> ./${project_name}/gwas/sam_selfPCA.cov
### logistic
plink --bfile ${in_8_4_1} \
--make-pheno ./${project_name}/gwas/sam.pheno "case" \
--ci 0.95 \
--logistic perm hide-covar \
--reference-allele ./${project_name}/gwas/alt_alleles \
--allow-no-sex --adjust --noweb \
--covar ./${project_name}/gwas/sam_selfPCA.cov \
--covar-name PC1-PC3 \
--out ${out_8_4_3} \
> ./${project_name}/log/8.4.3-assoc_test_logistic_with_selfPCA.log 2>&1

### sample information # 1000G PCA
echo -e "FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10" > ./${project_name}/gwas/sam_1000PCA.cov
awk 'BEGIN{FS=" ";OFS="\t";}{print $2,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
./${project_name}/samInfo/${project_name}_1000g_pca_after_samFlt.eigenvec \
>> ./${project_name}/gwas/sam_1000PCA.cov
### logistic
plink --bfile ${in_8_4_1} \
--make-pheno ./${project_name}/gwas/sam.pheno "case" \
--ci 0.95 \
--logistic perm hide-covar \
--reference-allele ./${project_name}/gwas/alt_alleles \
--allow-no-sex --adjust --noweb \
--covar ./${project_name}/gwas/sam_1000PCA.cov \
--covar-name PC1-PC3 \
--out ${out_8_4_4} \
> ./${project_name}/log/8.4.4-assoc_test_logistic_with_1000PCA.log 2>&1
# ### lasso
# plink --bfile ${in_8_4_1} \
# --make-pheno ./${project_name}/gwas/sam.pheno "case" \
# --ci 0.95 \
# --lasso hide-covar \
# --reference-allele ./${project_name}/gwas/alt_alleles \
# --allow-no-sex --adjust --noweb \
# --covar ./${project_name}/gwas/sam.cov \
# --covar-name PC1-PC5 \
# --out ${out_8_4_2} \
# > ./${project_name}/log/8.4.2-assoc_test_logistic.log 2>&1

## subgroup pathway-level burden
## indelRmrptSNPraw
srun="srun -n1 -N1 --exclusive"

parallel --xapply --header : $srun Rscript ./code/burden_pathway_ebio_gene.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-g /mnt/work/research/zhaohq/burden/source/HPO_MPO_MsigDB.RData \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.4-burden_pathway_subgroup_2.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 2 2 2 &

parallel --xapply --header : $srun Rscript ./code/burden_pathway_ebio_gene.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-g /mnt/work/research/zhaohq/burden/source/HPO_MPO_MsigDB.RData \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.4-burden_pathway_subgroup_1.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 1 1 1 &

srun="srun -n1 -N1 --exclusive"
parallel --xapply --header : $srun Rscript ./code/burden_pathway_ebio_gene.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-g /mnt/work/research/zhaohq/burden/source/HPO_MPO_MsigDB.RData \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.4-burden_pathway_subgroup_1.log \
2>&1 ::: lm syn \
::: externalLof 1e-04 ::: internalLof 3e-04 \
::: subgroupCol 1 &

parallel --xapply --header : $srun Rscript ./code/burden_pathway_ebio_gene.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-g /mnt/work/research/zhaohq/burden/source/HPO_MPO_MsigDB.RData \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.4-burden_pathway_subgroup_2.log \
2>&1 ::: lm syn \
::: externalLof 1e-04 ::: internalLof 3e-04 \
::: subgroupCol 2 &

parallel --xapply --header : $srun Rscript ./code/burden_pathway_ebio_gene.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-g /mnt/work/research/zhaohq/burden/source/HPO_MPO_MsigDB.RData \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.4-burden_pathway_subgroup_3.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 3 3 3 &

parallel --xapply --header : $srun Rscript ./code/burden_pathway_ebio_gene.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-g /mnt/work/research/zhaohq/burden/source/HPO_MPO_MsigDB.RData \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.4-burden_pathway_subgroup_4.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 4 4 4 &

parallel --xapply --header : $srun Rscript ./code/burden_pathway_ebio_gene.R -v {lm} -l \
./${project_name}/samInfo/sam_class_afFlt.txt \
-p ./${project_name}/burdenTest/qualified_variant/${project_name}.indelRmrptSNPraw \
-g /mnt/work/research/zhaohq/burden/source/HPO_MPO_MsigDB.RData \
-e {externalLof} -i {internalLof} -t refseq \
-s ./${project_name}/samInfo/subgroup.RData -n {subgroupCol} > \
./${project_name}/log/7.4-burden_pathway_subgroup_5.log \
2>&1 ::: lm lof lof lof lof_mis lof_mis lof_mis \
::: externalLof 1e-04 0.001 0.01 ::: internalLof 3e-04 0.001 0.01 \
::: subgroupCol 5 5 5 &

wait

### patient level gene distribution in a pathway
srun="srun -n1 -N1 --exclusive"
parallel --header : \
$srun Rscript ./code/patientLevel_pathwayGene_stat.R \
-j ./${project_name} \
-p {pathway} -s {subgroup} -v {variantType} \
-l "./samInfo/sam_class_afFlt.txt" \
-e {externalLof} -i {internalLof} \
-f "./samInfo/subgroup.RData" \
-g "/mnt/work/research/zhaohq/burden/cs_0615/samInfo/HPO4CS.RData" \
> ./${project_name}/log/7.6-patientLevel_pathwayGene_stat_HPO_1.log 2>&1 \
::: pathway "HPO_Abnormal_oral_morphology" \
"HPO_Abnormality_of_forebrain_morphology" "HPO_Hypotonia" \
"HPO_Global_developmental_delay" "HPO_Abnormality_of_the_vertebral_column" \
::: subgroup 1 2 3 \
::: externalLof 0.001 ::: internalLof 0.001 \
::: variantType lof lof_mis &

parallel --header : \
$srun Rscript ./code/patientLevel_pathwayGene_stat.R \
-j ./${project_name} \
-p {pathway} -s {subgroup} -v {variantType} \
-l "./samInfo/sam_class_afFlt.txt" \
-e {externalLof} -i {internalLof} \
-f "./samInfo/subgroup.RData" \
-g "/mnt/work/research/zhaohq/burden/cs_0615/samInfo/HPO4CS.RData" \
> ./${project_name}/log/7.6-patientLevel_pathwayGene_stat_HPO_2.log 2>&1 \
::: pathway "HPO_Abnormal_oral_morphology" \
"HPO_Abnormality_of_forebrain_morphology" "HPO_Hypotonia" \
"HPO_Global_developmental_delay" "HPO_Abnormality_of_the_vertebral_column" \
::: subgroup 1 2 3 \
::: externalLof 0.01 ::: internalLof 0.01 \
::: variantType lof lof_mis &

parallel --header : \
$srun Rscript ./code/patientLevel_pathwayGene_stat.R \
-j ./${project_name} \
-p {pathway} -s {subgroup} -v {variantType} \
-l "./samInfo/sam_class_afFlt.txt" \
-e {externalLof} -i {internalLof} \
-f "./samInfo/subgroup.RData" \
-g "/mnt/work/research/zhaohq/burden/cs_0615/samInfo/HPO4CS.RData" \
> ./${project_name}/log/7.6-patientLevel_pathwayGene_stat_HPO_3.log 2>&1 \
::: pathway "HPO_Abnormal_oral_morphology" \
"HPO_Abnormality_of_forebrain_morphology" "HPO_Hypotonia" \
"HPO_Global_developmental_delay" "HPO_Abnormality_of_the_vertebral_column" \
::: subgroup 1 2 3 \
::: externalLof 1e-04 ::: internalLof 3e-04 \
::: variantType lof lof_mis &