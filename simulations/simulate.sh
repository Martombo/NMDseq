rm simus
k=0
for k in {0,0.04,0.08,0.12,0.16,0.2,0.24,0.28}
do
	Rscript simulate.R $k
	join <(cut -d' ' -f1,4 SMG6_results) <(cut -d' ' -f-2 results_simu|sort) >tmp
	join tmp genes_exp >compare
	echo $k >>simus 
	cut -d' ' -f2,3 compare | sed 's/NA/0/g' | plotit_noRM pearson | grep -A1 -B2 "     cor" >>simus
done
