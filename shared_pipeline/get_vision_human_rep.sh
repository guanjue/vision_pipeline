#script_folder=/storage/home/gzx103/scratch/vision/human/data2run/shared_pipeline/
#input_folder=/storage/home/gzx103/scratch/vision/human/data2run/
###### Dependence
### R
# metap
### python
# numpy, matplotlib, scipy

script_folder=$1
input_folder=$2
ref_ref=$3
### H3K4me3_H1_R1
mark_ref_cell=$4
### H1
mark_ref_rep_id=$5

cd $input_folder

###### (1) Reads-count --> (2) NB-p-value 
#sleep 20000
for ct in $(cat cell_list.txt)
do
	echo $ct
	for mk in $(cat mark_list.txt)
	do
		echo $mk
		###### 2_nbp
		for r_id in $(cat replicate_list.txt)
		do
			time Rscript $script_folder'negative_binomial_p_2r_bgadj.R' $mk'_'$ct'_'$r_id'.ip.bed' $mk'_'$ct'_'$r_id'.ct.bed' $mk'_'$ct'_'$r_id
		done
	done
done

###### (3) Prepare PKnorm list (between reference sample)
rm 'pknorm_list_reference.rep.txt'
for mk in $(cat mark_list.txt)
do
	echo $mk
	### select the reference sample for different mark reference pknorm normalization
	echo $ref_ref'.nbp_2r_bgadj.txt' >> 'pknorm_list.1.txt'
	echo $mk'_'$mark_ref_cell'_'$mark_ref_rep_id'.nbp_2r_bgadj.txt' >> 'pknorm_list.2.txt'
done
paste 'pknorm_list.1.txt' 'pknorm_list.2.txt' > 'pknorm_list_reference.rep.txt'
rm 'pknorm_list.1.txt'
rm 'pknorm_list.2.txt'

###### (4) Prepare PKnorm list (between reference & target sample)
for mk in $(cat mark_list.txt)
do
	echo $mk
	rm $mk'.pknorm_list.rep.txt'
	for ct in $(cat cell_list.txt)
	do
		echo $ct
		echo $mk'_'$mark_ref_cell'_'$mark_ref_rep_id'.pknorm.ref.txt' >> $mk'.pknorm_list.1.txt'
		echo $mk'_'$ct'_R1.nbp_2r_bgadj.txt' >> $mk'.pknorm_list.2.txt'
		echo $mk'_'$mark_ref_cell'_'$mark_ref_rep_id'.pknorm.ref.txt' >> $mk'.pknorm_list.1.txt'
		echo $mk'_'$ct'_R2.nbp_2r_bgadj.txt' >> $mk'.pknorm_list.2.txt'
	done
	paste $mk'.pknorm_list.1.txt' $mk'.pknorm_list.2.txt' > $mk'.pknorm_list.rep.txt'
	rm $mk'.pknorm_list.1.txt' 
	rm $mk'.pknorm_list.2.txt'
done

###### (5) PKnorm normalization (between reference sample)
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1}')
	upperlim=500
	lowerlim=0
	echo $sig1 
	echo $sig2
	echo $sig2_celltype
	### set upper limit
	cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
	cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt'
	### peak norm
	time python $script_folder'peaknorm_rotate_log_ref_mean.py' -w 1 -p 1 -n 500000 -a 1 -b $sig1'.upperlim.txt' -c 1 -d $sig2'.upperlim.txt' -u $upperlim -l $lowerlim
	### rm tmp files
	rm $sig1'.upperlim.txt' $sig2'.upperlim.txt'
	mv $sig2_celltype'_nbp_2r_bgadj.pknorm.txt' $sig2_celltype'.pknorm.ref.txt'
	mv $sig2_celltype'_nbp_2r_bgadj.info.txt' $sig2_celltype'.info.ref.txt'
	mv $sig2_celltype'_nbp_2r_bgadj.pknorm.scatterplot.png' $sig2_celltype'.pknorm.scatterplot.ref.png'
	mv $sig2_celltype'_nbp_2r_bgadj.scatterplot.png' $sig2_celltype'.scatterplot.ref.png'
done < pknorm_list_reference.rep.txt

###### (6) PKnorm normalization (between reference & target sample)
for mk in $(cat mark_list.txt)
do
	echo $mk
	while read LINE
	do
		sig1=$(echo "$LINE" | awk '{print $1}')
		sig2=$(echo "$LINE" | awk '{print $2}')
		sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1}')
		upperlim=500
		lowerlim=0
		echo $sig1 
		echo $sig2
		echo $sig2_celltype
		### set upper limit
		cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
		cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt'
		### peak norm
		time python $script_folder'peaknorm_rotate_log_z_mean.py' -w 1 -p 1 -n 500000 -a 1 -b $sig1'.upperlim.txt' -c 1 -d $sig2'.upperlim.txt' -u $upperlim -l $lowerlim
		### rm tmp files
		rm $sig1'.upperlim.txt' $sig2'.upperlim.txt'
		cat $sig2_celltype'.pknorm.txt' | awk -F '\t' -v OFS='\t' '{if ($1>=16) print 16; else print $1 }' > $sig2_celltype'.pknorm.16lim.txt'
		cat $sig2_celltype'.pknorm.txt' | awk -F '\t' -v OFS='\t' '{if ($1>=16) print 16; else if ($1<=2) print 2; else print $1}' > $sig2_celltype'.pknorm.2_16lim.txt'
	done < $mk'.pknorm_list.rep.txt'
done

mkdir vision_data_output/
mv $sig2_celltype'.pknorm.16lim.txt' vision_data_output/



 
