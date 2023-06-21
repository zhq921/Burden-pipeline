# Produce PCA bi-plot for 1000 Genomes Phase III
# ref: https://www.biostars.org/p/335605/
plink=/share/pub/huangyk/conda/envs/cfDNA/bin/plink
bcftools=/share/pub/huangyk/conda/envs/cfDNA/bin/bcftools

# 1, Download the files as VCF.gz (and tab-indices)
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" ;
suffix=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" ;
for chr in {1..22}; do
    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done
# 2, Download 1000 Genomes PED file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;
# 3, Download the GRCh37 / hg19 reference genome
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai ;
gunzip human_g1k_v37.fasta.gz ;
# 4, Convert the 1000 Genomes files to BCF
# Ensure that multi-allelic calls are split and that indels are left-aligned compared to reference genome (1st pipe)
# Sets the ID field to a unique value: CHROM:POS:REF:ALT (2nd pipe)
# Removes duplicates (3rd pipe)
# -I +'%CHROM:%POS:%REF:%ALT' means that unset IDs will be set to CHROM:POS:REF:ALT
# 
# -x ID -I +'%CHROM:%POS:%REF:%ALT' first erases the current ID and then sets it to CHROM:POS:REF:ALT
for chr in {1..22}; do
    $bcftools norm -m-any --check-ref w -f human_g1k_v37.fasta \
      ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
      $bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        $bcftools norm -Ob --rm-dup both \
          > ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf ;

    $bcftools index ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf ;
done
# 5, Convert the BCF files to PLINK format
for chr in {1..22}; do
    $plink --noweb \
      --bcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;
done

# 7, Prune variants from each chromosome
# --maf 0.10, only retain SNPs with MAF greater than 10%
# --indep [window size] [step size/variant count)] [Variance inflation factor (VIF) threshold]
# e.g. indep 50 5 1.5, Generates a list of markers in approx. linkage equilibrium - takes 50 SNPs at a time and then shifts by 5 for the window. VIF (1/(1-r^2)) is the cut-off for linkage disequilibrium
# 
mkdir Pruned ;

for chr in {1..22}; do
    $plink --noweb \
      --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
      --maf 0.10 --indep 50 5 1.5 \
      --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;

    $plink --noweb \
      --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
      --extract Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.prune.in \
      --make-bed \
      --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;
done
# 8, Get a list of all PLINK files
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;
sed -i 's/.bim//g' ForMerge.list ;
# 9, Merge all projects into a single PLINK file
$plink --merge-list ForMerge.list --out Merge ;