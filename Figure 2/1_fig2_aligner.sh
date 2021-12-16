for ((i=1;i<=149;i++)); do
    bwa mem fasta/IL10RA.fasta del_$i.fastq.gz | grep -vE "^@" | cut -f6 > del_$i.bwa_mem
    bwa mem fasta/IL10RA.fasta ins_$i.fastq.gz | grep -vE "^@" | cut -f6 > ins_$i.bwa_mem
    bowtie2 -x fasta/IL10RA -U del_$i.fastq.gz | grep -vE "^@" | cut -f6 > del_$i.bt2
    bowtie2 -x fasta/IL10RA -U ins_$i.fastq.gz | grep -vE "^@" | cut -f6 > ins_$i.bt2
done


