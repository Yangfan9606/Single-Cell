# INPUT with read1.fastq.gz, read2.fastq.gz
INPUT="input_sample_name"

perl Hydrogel_seq_Read1_Barcode.check.pl $INPUT.read1.fastq.gz -b1 barcode1_TSO.txt -b2 barcode2.txt -b3 barcode3.txt -o $INPUT.R1
perl trim.fastq.pl $INPUT.R1.BC_check.gz $INPUT.read1.fastq.gz $INPUT.read2.fastq.gz -o $INPUT
perl fix.R1R2.clipper.pl $INPUT.1.fastq.clipper.gz $INPUT.2.fastq.clipper.gz -o $INPUT
perl To.10x.fastq.pl 3M-february-2018.txt.gz $INPUT.2.fix.fastq.gz $INPUT.R1.BC_check.gz $INPUT.1.fix.fastq.gz -o $INPUT
ln -s $INPUT.For10x.1.fastq.gz ${INPUT}_S1_L001_R1_001.fastq.gz
ln -s $INPUT.For10x.2.fastq.gz ${INPUT}_S1_L001_R2_001.fastq.gz

# Cteate mix genome in GRCh38_mm10_mix_ref use cellranger mkref
# ---
# Run with cellranger 
INPUT_ID="sample_ID"
cellranger count --id=$INPUT_ID --transcriptome=./GRCh38_mm10_mix_ref/ --fastqs=./ --sample=$INPUT --chemistry=SC3Pv3 --localcores 12
