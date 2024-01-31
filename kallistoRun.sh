wget http://ftp.ensembl.org/pub/release-105/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-105/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.105.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/ena/assembly/GCA_000/GCA_000001/GCA_000001215.4_sequence_report.txt

head -n 8 GCA_000001215.4_sequence_report.txt | cut -f 2-3 |tail -n+2 > Dm_BDGP6.32.chromosomes.txt

kallisto index -i Dm_index.idx Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz

OE_DIR=/home/wheelerlab3/Data/Mierisch_Lab/Novogene_2022-01-19/raw_data

mkdir oe_results

echo "sample condition path" > oe_kallisto_results.txt

for sample in Rib_1 Rib_2 Rib_3 w1118_1 w1118_2 w1118_3
do
    kallisto quant -i Dm_index.idx --genomebam --gtf Drosophila_melanogaster.BDGP6.32.105.gtf.gz --chromosomes Dm_BDGP6.32.chromosomes.txt \
     -o oe_results/$sample -b 100 -t 4 ${OE_DIR}/${sample}/${sample}_1.fq.gz ${OE_DIR}/${sample}/${sample}_2.fq.gz

    [[ $sample = Rib* ]] && echo "$sample Rib oe_results/$sample" >> oe_kallisto_results.txt
    [[ $sample = w1118* ]] && echo "$sample control oe_results/$sample" >> oe_kallisto_results.txt
done

zcat /home/madeline/Research/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz | grep ">" |cut -d' ' -f1,4,7,8- |cut -d'>' -f2 > tmp_annotation.txt
python3 /home/wheelerlab3/Data/Mierisch_Lab/heather_analysis_2022-03-07/format_Dm_annotation_file.py
rm tmp_annotation.txt
