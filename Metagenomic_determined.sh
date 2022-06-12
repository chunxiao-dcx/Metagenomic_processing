###########################################################
# File Name: Metagenomic_determined.sh
# Publication: Metagenomic analysis reveals indole signaling effect on bacteria community in activated sludge: quorum sensing inhibition and antibiotic resistance enrichment
# Author: Dai Chunxiao
# mail: 2446378845@qq.com
# Group: A-G1, B-G2, C-G3
###########################################################

for i in A1 A2 A3 B1 B2 B3 C1 C2 C3; do 
{
trimmomatic PE \
-threads 12 \
-phred33 \
./data/${i}_1.fastq.gz ./data/${i}_1.fastq.gz \
./trim/${i}_1.fq.gz  ./data/output_forward_unpaired_${i}_1.fq.gz ./trim/${i}_2.fq.gz ./data/output_reverse_unpaired_${i}_2.fq.gz \
ILLUMINACLIP:./soft/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:true \
SLIDINGWINDOW:4:15 \
LEADING:3 \
TRAILING:3 \
MINLEN:36 
}
done

for i in A1 A2 A3 B1 B2 B3 C1 C2 C3; do 
{
./megahit -1 ./trim/${i}_1.fq.gz -2 ./trim/${i}_2.fq.gz \
-o ./trim/MY_OUTPUT_${i}
}
done

for i in A1 A2 A3 B1 B2 B3 C1 C2 C3; do 
{
./gmhmmp -m MetaGeneMark_v1.mod -a -d \
-o ${i}.gff -f G \
-A ${i}.faa \
-D ${i}.gene.ffn -r -s . \
-g 11 -p 1 ./trim/MY_OUTPUT_${i}/final.contigs.fasta
}
done
# Gene abundance by salmon 
for i in A1 A2 A3 B1 B2 B3 C1 C2 C3; do 
{
salmon index -t ${i}.gene.ffn -p 9 -k 31 -i ./index${i}
salmon quant --validateMappings -i ./index${i} -l A -p 4 -1 ./Trim/${i}_1.fq.gz -2 ./Trim/${i}_2.fq.gz -o ./salmon/${i}.quant
}
done

# Annotation using Megan
for i in A1 A2 A3 B1 B2 B3 C1 C2 C3; do
{
diamond blastx --db ./database/nr -q ./${i}.gene.ffn -t /tmp --daa ./megan/${i}_NR.daa
daa2rma -i ./megan/${i}_NR.daa -ms 50 -me 0.01 -top 50  -mdb ./db/megan-map-Oct2019.db  -o ./megan/${i}.rma
rma2info -i ./megan/${i}.rma -r2c EGGNOG -n true --paths true --ranks true --list true --listMore true -v >  ./megan/${i}eggnog.txt
rma2info -i ./megan/${i}.rma -r2c Taxonomy -n true --paths true --ranks true --list true --listMore true -v >  ./megan/${i}Taxonomy.txt
}
done

# Quorum sensing genes prediction
for i in A1 A2 A3 B1 B2 B3 C1 C2 C3; do
{
diamond blastx --db ./db/QSG -q ./${i}.gene.faa -o ./QSG/${i}_QS_DB_matches.txt --more-sensitive --top 1 -e 0.00001 --id 50
}
done

# 16s DNA extraction
for i in A1 A2 A3 B1 B2 B3 C1 C2 C3; do
{
gzip -dc ./trim/${i}_1.fq.gz > ./trim/${i}_clean.R1.fq
gzip -dc ./trim/${i}_2.fq.gz > ./trim/${i}_clean.R2.fq

/usr/bin/python2 emirge.py RESULT/16S_MEC -1 ./trim/${i}_clean.R1.fq -2 ./trim/${i}_clean.R2.fq \
-f ../emirge_16S/silva_bt2/ge1200bp.le2000bp.0.97.fixed.SILVA_138_SSURef_NR99_tax_silva_trunc \
-b ../emirge_16S/silva_bt2/ge1200bp.le2000bp.0.97.fixed -a 64 --phred33
}
done

# Fungi and protozoa communities’ composition determined using Kraken 2
kraken2-build  --threads 8 --db ./database/kraken --download-taxonomy &&
kraken2-build  --threads 8 --db ./database/kraken  --download-library protozoa &&
kraken2-build  --threads 8 --db ./database/kraken  --download-library fungi &&

for i in A1 A2 A3 B1 B2 B3 C1 C2 C3; do
{
kraken2 --db ./database/kraken --threads 8  --report ./kraken/${i}.report --output ./kraken/${i}.output --paired ./trim/${i}_clean.R1.fq ./trim/${i}_clean.R1.fq
}
done

# ARGs prediction using ARG-OAP
perl argoap_pipeline_stageone_version2.3 -i ./trim -o ./testoutdir -m ./dcx-meta-data.txt -n 8
perl argoap_pipeline_stagetwo_version2 -i ./testoutdir/extracted.fa -m ./testoutdir/meta_data_online.txt -o ./ARG-resultout -l 25 -d 80 -e 1e-5

