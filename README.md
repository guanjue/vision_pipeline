# vision_pipeline

## Dependence:
#### Python/2.7
###### numpy, matplotlib, scipy
#### R



## Input data
###### cell_list.txt in input folder
```
cell_list.txt
>>> head cell_list.txt 
H1
H2
CD34

```

###### mark_list.txt in input folder
```
mark_list.txt
>>> head mark_list.txt 
H3K4me1
H3K4me3
H3K27ac

```

## run pipeline
###### The 1st parameter is the script folder of the pipeline; 2nd parameter is the input folder of the data
```
time sh /storage/home/gzx103/scratch/vision/human/data2run/shared_pipeline/get_vision_human_rep.sh /storage/home/gzx103/scratch/vision/human/data2run/shared_pipeline/ /storage/home/gzx103/scratch/vision/human/data2run/test_data/ H3K4me3_H1_R1 H1 R1
```

