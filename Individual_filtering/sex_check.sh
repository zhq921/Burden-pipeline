
$plink --bfile ${out_4_5_4} --split-x 'b37' --make-bed --out ${out_4_5_4}.split-x \
> ./${project_name}/log/4.5.4.1-split_x.log 2>&1

$plink --bfile ${out_4_5_4} --check-sex --out ${out_4_5_4}.chech_sex

## XYdepth，indexedbam
## reference： https://www.biostars.org/p/192157/
### guess sex from two tests. 
### First test is weather the average coverage of chrX is greater 
### than minimum of average coverages of autosomal chromosomes 
### then the sample is probably female. Second test is weather 
### the average coverage of chrY is 10 times smaller than average 
### coverage of autosomal chromosomes then the sample is probably female.

coveragestats - Script to print filename and coverage statistics from samtools indexed files

echo $*
samtools idxstats $* 2>/dev/null | awk '{ tmp=($3)/($2) ; printf"%s %0.5f\n", $1, tmp }' 2>/dev/null | head -21
autosomalmin - Script to print minimum average coverage of autosomal chromosomes

coveragestats $* | cut -d' ' -f2 | head -20 | tail -19 | awk 'NR == 1 || $1 < min {min = $1}END{printf "%f\n", min}'
xcoverage - Script to print average coverage of chrX

coveragestats $* | grep ^X | cut -d" " -f2
ycoverage - Script to print average coverage of chrY

coveragestats $* | grep ^Y | cut -d" " -f2
maleorfemale - Script to perform first two described sex tests

min=$(autosomalmin $*)
xcov=$(xcoverage $*)
ycov=$(ycoverage $*)
test=$(echo "$xcov"'>'"$min" | bc -l)
if [ $test -eq 1 ]
then
    echo "$*" is female
else
    echo "$*" is male
fi
test=$(echo "$ycov"'<('"$min"/10\) | bc -l)
if [ $test -eq 1 ]
then
        echo "$*" is female
else
        echo "$*" is male
fi
