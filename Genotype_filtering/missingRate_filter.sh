#!/usr/bin/env bash

start=$(date +%s) 

if [[ -z "$3" ]]; then
echo -e "Description: Missing rate filter for two groups in a gvcf file by missing rate threshold\nUsage is missing_filter.sh gvcffile missing_rate_threshold_per_group name_for_output cont_num case_num"
exit 1
fi

echo "Input vcf: ${1}"
echo "Output vcf: ${3}"
echo "Max missing rate for each group: ${2}"
echo "Sample size in control group: ${4}"
echo "Sample size in case group: ${5}"

# th: threshold
g1_miss_th=$2
g2_miss_th=$2
cont_num=$4
case_num=$5
miss_th=1

awk 'BEGIN{
	FS="\t";OFS="\t"; # 
	g1_stt=10; g1_end=g1_stt+'"$cont_num"'-1; # group 1，9
	g2_stt=g1_stt+'"$cont_num"'; g2_end=g1_stt+'"$cont_num"'+'"$case_num"'-1;
	g1_n=g1_end-g1_stt+1; # group 1
	g2_n=g2_end-g2_stt+1;
	all_n=g1_n+g2_n
	}
	{
		if($1~/^#/) {print $0;} else{ # #
			g1_miss_n=0; g2_miss_n=0;
			for(i=1;i<=NF;i++){ # 
				if($i~/^\.\/\./){ # ./.
					if(i<=g1_end){g1_miss_n=g1_miss_n+1}
					if(i>=g2_stt){g2_miss_n=g2_miss_n+1}
				}
			}
			g1_miss_rt=g1_miss_n/g1_n; # group 1
			g2_miss_rt=g2_miss_n/g2_n;
			# miss_n=g1_miss_n+g2_miss_n;
			# miss_rt=miss_n/all_n;
			if(g1_miss_rt<='"$g1_miss_th"' && g2_miss_rt<='"$g2_miss_th"'){
				print $0
			}
		}
}' $1 > $3

end=$(date +%s)
take=$(( end - start ))
echo Time taken to execute commands is ${take} seconds.

