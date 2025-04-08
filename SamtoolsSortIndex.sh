#!/bin/bash
#Script to resort the .bam files by coordinate


#Parameters to configure:
#{project_root}: folder where all the files for the projects are stored
project_root="/mnt/data2/Clara/NGS_274_825/pipeline/align_hisat/sorted_bams/"
#{samtools_path}
samtools_path="/mnt/software/samtools/samtools-1.17"
#{output_dir}
output_dir="/mnt/data2/Clara/NGS_274_825/pipeline/align_hisat/sorted_bams/sorted_c"

######Project files
#Paht to the folder where the BAM files are
BAM_path=${project_root}

#Define file lists to loop over
BAM_File_List=$(ls ${BAM_path}/*.bam)



#########Run Picard tools over the files


	#Reorder cleanBAM files by sorting by queryname
    for BAM in ${BAM_File_List}; do
		
		
		#Define the Sample_Name
        Sample_Name=$(basename ${BAM} _R1_250115.sortn.bam)


        echo "Samtools sort Started working on ${BAM}"
	
	${samtools_path}/samtools sort ${BAM} -@ 6 -o ${output_dir}/${Sample_Name}.sortc.bam

	#${samtools_path}/samtools index ${output_dir}/${Sample_Name}.sortc.bam

    done
