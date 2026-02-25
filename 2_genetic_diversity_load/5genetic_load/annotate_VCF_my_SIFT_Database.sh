#step1 make database
#cd /home/zhangwp/Carya/mappingout_Cca202201/mutation_road/sift4g/SIFT4G_Create_Genomic_DB-master
#mkdir -p ./test_files/ChuNE/combine_all
#cd ./test_files/ChuNE/combine_all
#find ../ChuNE_* -name "v20241218" -type d -exec echo "ln -s {}/* . " \; |sh
#step1 make database


#step2 annotate for every ind vcf, filter depth, Q20, Bi-allele, 20% missing
java="/home/ZhangWP/software/SIFT4G_Create_Genomic_DB-master"
vcffile="/home/ZhangWP/water_pine/mutation_road/new_4Jun2025/Gpen147_out2_20_missing/sp_ind_vcf"
database="/home/ZhangWP/software/SIFT4G_Create_Genomic_DB-master/test_files/Gpen/combine_all"
results="/home/ZhangWP/water_pine/mutation_road/new_4Jun2025/Gpen147_out2_20_missing/sift4g_result"

for id in `cat $1`
do
java -jar $java/SIFT4G_Annotator.jar -c -i $vcffile/${id}.recode.vcf -d $database -r $results -t
#java -Xms2G -Xmx8G -jar $java/SIFT4G_Annotator.jar -c -i $vcffile/${id}.recode.vcf -d $database -r $results -t
#java -Xms2G -Xmx8G -XX:ActiveProcessorCount=10 -jar $java/SIFT4G_Annotator.jar -c -i $vcffile/${id}.recode.vcf -d $database -r $results -t

done

#`-c` 表示使用坐标模式（Coordinate mode），即根据VCF文件中的染色体和位置信息来注释变异（必须的参数）
#`-i` 表示输入vcf文件的路径
#`-d` 表示SIFT4G数据库目录的路径
#`-r` 表示结果文件夹的路径
#`-t` 表示提取多个转录本的注释（Optional），可选参数
