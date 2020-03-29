docker run -i -t -v /home/blanca/tutorial:/tutorial -w /tutorial ceciliaklein/teaching:uvic

export PATH=$PATH:/tutorial/teaching-utils/;

1. Differential expression analysis edgeR

cd analysis

edgeR.analysis.R --input_matrix ../quantifications/encode.mouse.gene.expected_count.idr_NA.tsv \
                 --metadata /tutorial/data/gene.quantifications.index.tsv \
                 --fields tissue \
                 --coefficient 3 \
                 --output brain_X_liver

awk '$NF<0.01 && $2<-10{print $1"\tover_brain_X_liver"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.over_brain_X_liver.txt

awk '$NF<0.01 && $2>10 {print $1"\tover_liver_X_brain"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.over_liver_X_brain.txt

wc -l edgeR.over*.txt

Output: 

 118 edgeR.over_brain_X_liver.txt
 107 edgeR.over_liver_X_brain.txt
 225 total

- Heatmap

awk '$3=="gene"{ match($0, /gene_id "([^"]+).+gene_type "([^"]+)/, var); print var[1],var[2] }' OFS="\t" /tutorial/refs/gencode.vM4.gtf \
| join.py --file1 stdin \
          --file2 <(cat edgeR.over*.txt) \
| sed '1igene\tedgeR\tgene_type' > gene.edgeR.2.tsv


cut -f1 gene.edgeR.2.tsv \
| tail -n+2 \
| selectMatrixRows.sh - ../quantifications/encode.mouse.gene.TPM.idr_NA.tsv \
| ggheatmap.R --width 5 \
              --height 8 \
              --col_metadata /tutorial/data/gene.quantifications.index.tsv \
              --colSide_by tissue \
              --col_labels labExpId \
              --row_metadata gene.edgeR.2.tsv \
              --merge_row_mdata_on gene \
              --rowSide_by edgeR,gene_type \
              --row_labels none \
              --log \
              --pseudocount 0.1 \
              --col_dendro \
              --row_dendro \
              --matrix_palette /tutorial/palettes/palDiverging.txt \
              --colSide_palette /tutorial/palettes/palTissue.txt \
              --output heatmap.brain_X_liver.pdf

2. GO enrichment BP 

awk '{split($10,a,/\"|\./); print a[2]}' /tutorial/refs/gencode.vM4.gtf | sort -u > universe.txt

# BP brain X liver
awk '{split($1,a,"."); print a[1]}' edgeR.over_brain_X_liver.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_brain_X_liver \
                  --species mouse

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_brain_X_liver.BP.tsv

% GOBPID Pvalue
GO:0060322 3e-28
GO:0007417 9e-16
GO:0006355 9e-14
GO:2001141 1e-13
GO:0021543 2e-12
GO:0050767 8e-12
GO:0034654 1e-11
GO:0019219 2e-11
GO:0045944 4e-11
GO:0032502 4e-11
GO:0022008 1e-10
GO:1903508 2e-10
GO:0001764 2e-10
GO:0021983 3e-10
GO:0051254 4e-10
GO:0032501 5e-10
GO:0007275 5e-10
GO:1901576 2e-09
GO:0010628 4e-09
GO:0009889 4e-09
GO:0090304 5e-09
GO:0010557 5e-09
GO:0034645 6e-09
GO:0021879 7e-09
GO:0060255 9e-09
GO:0021798 1e-08
GO:2000026 2e-08
GO:0031328 3e-08
GO:0060041 3e-08
GO:0051173 4e-08
GO:0021772 4e-08
GO:0045595 7e-08
GO:0090596 7e-08
GO:0001654 1e-07
GO:0048664 1e-07
GO:0051961 1e-07
GO:0048513 2e-07
GO:0051962 2e-07
GO:0034641 4e-07
GO:0022037 4e-07
GO:0048665 4e-07
GO:0046483 4e-07
GO:0000904 4e-07
GO:0021537 4e-07
GO:0045165 5e-07
GO:0048858 5e-07
GO:1901360 6e-07
GO:0006725 6e-07
GO:0031175 6e-07
GO:0045665 8e-07
GO:0021779 9e-07
GO:0021780 9e-07
GO:0021854 1e-06
GO:0010720 2e-06
GO:0007411 2e-06
GO:0045666 2e-06
GO:0021954 2e-06
GO:0021978 3e-06
GO:0021511 6e-06
GO:0021695 6e-06
GO:2000177 6e-06
GO:0021796 7e-06
GO:0044708 9e-06
GO:0048562 1e-05
GO:0030855 1e-05
GO:0021846 1e-05
GO:0035881 1e-05
GO:0032989 2e-05
GO:0010721 2e-05
GO:0021885 2e-05
GO:0048708 2e-05
GO:0048666 2e-05
GO:0050789 2e-05
GO:0003326 2e-05
GO:0003329 2e-05
GO:0070593 2e-05
GO:0010629 2e-05
GO:0045685 2e-05
GO:0048593 3e-05
GO:0021522 3e-05
GO:0021527 3e-05
GO:0072091 3e-05
GO:0010001 4e-05
GO:0030182 4e-05
GO:0060579 5e-05
GO:0048646 5e-05
GO:0048869 5e-05
GO:0021510 6e-05
GO:0021530 6e-05
GO:0043583 7e-05
GO:0048522 7e-05
GO:0021542 7e-05
GO:0009953 7e-05
GO:0019827 7e-05
GO:2000113 7e-05
GO:0048870 8e-05
GO:0009890 9e-05
GO:0021575 9e-05
GO:0021795 1e-04
GO:1903507 1e-04
GO:0048663 1e-04
GO:0042472 1e-04
GO:0007610 1e-04
GO:0021895 1e-04
GO:0000122 1e-04
GO:0021761 2e-04
GO:0048714 2e-04
GO:0051253 2e-04
GO:0021764 2e-04
GO:0021797 2e-04
GO:0002052 2e-04
GO:0048863 2e-04
GO:0007548 2e-04
GO:0002088 2e-04
GO:0007492 3e-04
GO:0007224 3e-04
GO:0021516 3e-04
GO:0021960 3e-04
GO:0048852 3e-04
GO:0051172 3e-04
GO:0030900 3e-04
GO:0014003 4e-04
GO:0035883 4e-04
GO:0021877 4e-04
GO:0009790 4e-04
GO:0021540 4e-04
GO:0021902 4e-04
GO:0060235 4e-04
GO:0098700 4e-04
GO:0003002 4e-04
GO:0007409 4e-04

# BP liver X brain
awk '{split($1,a,"."); print a[1]}' edgeR.over_liver_X_brain.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_liver_X_brain \
                  --species mouse

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_liver_X_brain.BP.tsv

% GOBPID Pvalue
GO:1900047 2e-18
GO:0050819 4e-18
GO:0061045 8e-17
GO:0030193 2e-16
GO:1903034 1e-12
GO:0050878 4e-12
GO:0055088 1e-11
GO:0010466 4e-11
GO:0031639 1e-10
GO:0046464 1e-10
GO:0070328 2e-10
GO:1900048 2e-09
GO:0042632 3e-09
GO:0050820 4e-09
GO:0006638 5e-09
GO:0051336 7e-09
GO:0019752 9e-09
GO:0097006 1e-08
GO:0051918 2e-08
GO:0034368 4e-08
GO:0072378 4e-08
GO:0033344 4e-08
GO:0006953 4e-08
GO:0015918 5e-08
GO:0010898 7e-08
GO:0042730 8e-08
GO:0043086 9e-08
GO:0072376 1e-07
GO:0006631 1e-07
GO:0051919 1e-07
GO:0034375 2e-07
GO:0090303 2e-07
GO:0006958 3e-07
GO:1905952 3e-07
GO:0045834 5e-07
GO:0071827 6e-07
GO:0008202 7e-07
GO:0009605 8e-07
GO:0043434 9e-07
GO:0045723 2e-06
GO:0030162 2e-06
GO:0071704 2e-06
GO:0048585 2e-06
GO:0090207 2e-06
GO:0006959 6e-06
GO:0034384 6e-06
GO:0030212 6e-06
GO:0051180 7e-06
GO:0010873 9e-06
GO:0034370 9e-06
GO:0034374 9e-06
GO:0050994 9e-06
GO:0046486 1e-05
GO:0035296 1e-05
GO:0044710 1e-05
GO:0034380 1e-05
GO:0051006 1e-05
GO:0044712 1e-05
GO:0032269 2e-05
GO:0050880 2e-05
GO:0033700 2e-05
GO:0034116 2e-05
GO:0044242 2e-05
GO:0010903 2e-05
GO:0010876 2e-05
GO:1902578 3e-05
GO:0016053 3e-05
GO:0043691 3e-05
GO:0006954 3e-05
GO:0050708 4e-05
GO:0034434 4e-05
GO:0006950 4e-05
GO:0044711 5e-05
GO:0051235 7e-05
GO:0003331 7e-05
GO:0043152 7e-05
GO:0016485 7e-05
GO:0046890 9e-05
GO:1901700 9e-05
GO:0051239 9e-05
GO:0044238 1e-04
GO:0010038 1e-04
GO:0048584 1e-04
GO:0022600 1e-04
GO:0006082 1e-04
GO:0034383 1e-04
GO:0060193 1e-04
GO:0097756 1e-04
GO:0010884 2e-04
GO:0048523 2e-04
GO:0090277 2e-04
GO:0016064 2e-04
GO:0072330 2e-04
GO:0019369 2e-04
GO:0031100 2e-04
GO:0065005 2e-04
GO:0008203 2e-04
GO:0002740 2e-04
GO:0006548 2e-04
GO:0019538 3e-04
GO:0042592 3e-04
GO:0070374 3e-04
GO:0065009 3e-04
GO:1902042 3e-04
GO:2000352 3e-04
GO:0002034 3e-04
GO:0010534 3e-04
GO:0010757 3e-04
GO:0034382 3e-04
GO:0002250 3e-04
GO:0002790 4e-04
GO:0046883 4e-04
GO:0007548 4e-04
GO:0040008 4e-04
GO:0006810 4e-04
GO:0070613 4e-04
GO:0009892 5e-04
GO:0008406 5e-04

3. Differential splicing analysis 

cd ../splicing

3.1. Prepare input files 

# List of transcript ids
awk '$3=="transcript" && $0~/gene_type "protein_coding"/{ match($0, /transcript_id "([^"]+)/, id); print id[1] }' /tutorial/refs/gencode.vM4.gtf |sort -u > protein_coding_transcript_IDs.txt

# Genome annotation restricted to exon features and filtered by transcript type
cat /tutorial/refs/gencode.vM4.gtf |awk '$3=="exon"' |grep -Ff protein_coding_transcript_IDs.txt > exon-annot.gtf

# Filter transcript TPM matrix
selectMatrixRows.sh protein_coding_transcript_IDs.txt /tutorial/quantifications/encode.mouse.transcript.TPM.idr_NA.tsv > pc-tx.tsv

# Individual transcript expression matrices
for tissue in Brain Heart Liver; do
    selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 pc-tx.tsv > expr.${tissue}.tsv
done

3.2. Generate local AS events

SUPPA info in github:  https://github.com/comprna/SUPPA/blob/master/README.md#overview 


suppa.py generateEvents -i exon-annot.gtf -e SE SS MX RI FL -o localEvents -f ioe

wc -l localEvents*ioe

 6200 localEvents_A3_strict.ioe
    5583 localEvents_A5_strict.ioe
   20128 localEvents_AF_strict.ioe
    3401 localEvents_AL_strict.ioe
    1054 localEvents_MX_strict.ioe
    2763 localEvents_RI_strict.ioe
   12254 localEvents_SE_strict.ioe
   51383 total

3.3. Compute percent spliced in index (PSI) values for local events

Compute PSI for all events and create individuals files per tissue:
    
-Skipping Exon (SE):

event=SE; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=SE; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

- Retained Intron (RI):

event=RI; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=RI; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

- Mutually Exclusive Exons (MX):

event=MX; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=MX; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

- Alternative First Exons (AF):

event=AF; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}

event=AF; for tissue in Brain Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done


3.4. Differential splicing analysis for local events

- Compute differential splicing analysis for all events with SUPPA:

event=SE; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=RI; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=MX; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

event=AF; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

- Look at top cases and play with thresholds of p-value and delta PSI:

Con estos thresholds no sale ninguno:

event=SE; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.7 || $2<-0.7) && $3<0.05{print}' DS.${event}.dpsi

event=RI; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.7 || $2<-0.7) && $3<0.05{print}' DS.${event}.dpsi

event=MX; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.7 || $2<-0.7) && $3<0.05{print}' DS.${event}.dpsi

event=AF; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.7 || $2<-0.7) && $3<0.05{print}' DS.${event}.dpsi

We try less restrictive ones: 

event=SE; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.55 || $2<-0.55) && $3<0.05{print}' DS.${event}.dpsi

#Output: 
#Brain-Liver_dPSI	Brain-Liver_p-val
#ENSMUSG00000019302.12;SE:chr11:101048525-101049459:101049476-101055026:+	-0.5880491647	0.0
#ENSMUSG00000026276.13;SE:chr1:93479068-93489124:93489154-93491157:+	0.5531274965999999	0.0
#ENSMUSG00000073468.7;SE:chr17:8311285-8317850:8317871-8318844:+	-0.5941274244	0.0

event=RI; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi

#Output:

#Brain-Liver_dPSI	Brain-Liver_p-val
#ENSMUSG00000002010.13;RI:chrX:73782006:73782103-73782675:73782716:-	-0.3174466202	0.0
#ENSMUSG00000034211.10;RI:chr5:129715684:129715910-129716727:129716866:+	-0.39104938899999997	0.0
#ENSMUSG00000037933.11;RI:chr13:49383043:49383253-49384971:49387025:+	-0.5230577351	0.0
#ENSMUSG00000038241.12;RI:chr2:155988312:155988432-155989499:155989620:+	0.4438048007	0.0
#ENSMUSG00000060950.7;RI:chr12:111678105:111678300-111678603:111678962:+	-0.35722050590000004	0.0


event=MX; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi

#Output:
#Brain-Liver_dPSI	Brain-Liver_p-val
#ENSMUSG00000019843.10;MX:chr10:39526930-39529437:39529601-39532016:39526930-39529919:39530074-39532016:+	-0.3312477193	0.0
#ENSMUSG00000020167.10;MX:chr10:80410274-80410398:80410624-80413205:80410274-80412827:80413062-80413205:-	0.3931814618	0.0
#ENSMUSG00000028863.9;MX:chr4:125103031-125107517:125107546-125108323:125103031-125107652:125107709-125108323:+	-0.4476539615	0.0
#ENSMUSG00000033335.13;MX:chr9:21478348-21478899:21479037-21481333:21478348-21480407:21480545-21481333:+	0.3322162528	0.0
#ENSMUSG00000040774.11;MX:chr3:106505852-106512803:106512887-106531137:106505852-106520009:106520093-106531137:-	-0.35759227200000004	0.0
#ENSMUSG00000053436.10;MX:chr17:28728975-28737011:28737090-28741788:28728975-28740680:28740759-28741788:+	0.3119375901	0.0
#ENSMUSG00000054808.9;MX:chr7:28907064-28909903:28909988-28912240:28907064-28911490:28911575-28912240:-	-0.30150284969999996	0.0


event=AF; awk 'BEGIN{FS=OFS="\t"}NR==1{print }NR>1 && $2!="nan" && ($2>0.6 || $2<-0.6) && $3<0.05{print}' DS.${event}.dpsi

#Output: 
#Brain-Liver_dPSI	Brain-Liver_p-val
#ENSMUSG00000002345.13;AF:chr8:70139707:70139772-70141864:70140277:70140356-70141864:+	-0.6110412516	0.000999001
#ENSMUSG00000029922.11;AF:chr6:39410216-39419840:39420067:39410216-39420099:39420358:-	-0.6698497643	0.0
#ENSMUSG00000032126.11;AF:chr9:44342071-44342221:44342381:44342071-44344010:44344228:-	0.7026755148	0.0
#ENSMUSG00000036002.8;AF:chr4:43037206-43039966:43040288:43037206-43045648:43045797:-	0.7832400683	0.0
#ENSMUSG00000042066.11;AF:chr1:132381181-132389902:132390318:132381181-132390977:132391281:-	0.7294985313	0.0
#ENSMUSG00000046598.10;AF:chr16:31422280:31422512-31437992:31436091:31436166-31437992:+	0.6137716352	0.0007492507


Heatmap:

SE: 

# prepare input for heatmap
event=SE; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.5 || $2<-0.5) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-SE.txt
selectMatrixRows.sh top-examples-SE.txt DS.SE.psivec > matrix.top-examples-SE.tsv

# heatmap SE
ggheatmap.R -i matrix.top-examples-SE.tsv -o heatmap_top-examples-SE.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8

RI: 

# prepare input for heatmap
event=RI; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-RI.txt
selectMatrixRows.sh top-examples-RI.txt DS.RI.psivec > matrix.top-examples-RI.tsv

# heatmap SE
ggheatmap.R -i matrix.top-examples-RI.tsv -o heatmap_top-examples-RI.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8


MX: 

# prepare input for heatmap
event=MX; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-MX.txt
selectMatrixRows.sh top-examples-MX.txt DS.MX.psivec > matrix.top-examples-MX.tsv

# heatmap SE
ggheatmap.R -i matrix.top-examples-MX.tsv -o heatmap_top-examples-MX.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8

AF: 

# prepare input for heatmap
event=AF; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.6 || $2<-0.6) && $3<0.05{print}' DS.${event}.dpsi|cut -f1 > top-examples-AF.txt
selectMatrixRows.sh top-examples-AF.txt DS.AF.psivec > matrix.top-examples-AF.tsv

# heatmap SE
ggheatmap.R -i matrix.top-examples-AF.tsv -o heatmap_top-examples-AF.pdf --matrix_palette /tutorial/palettes/palSequential.txt --row_dendro  --matrix_fill_limits "0,1" -B 8

##########
# ChIP-seq #
##########

mkdir ../chip-analysis

cd chip-analysis

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoHeart.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak","Heart_coordinates","Heart_peak","intersection"}$NF!=0{printf ("%20s\t%15s\t%15s\t%10s\t%20s\t%15s\t%15s\t%10s\t\n", $1,$2,$3,$4,$11,$12,$13,$14)}' > common-peaks-Brain-heart.tsv

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoHeart.narrowPeak -wao|awk 'BEGIN{OFS="\t"; print "Brain_coordinates","Brain_peak","Heart_coordinates","Heart_peak","intersection"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14}' > Common-peaks-Brain-liver.bed

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak"}$NF==0{print $1":"$2"-"$3,$4}' > Brain-specific-peaks.tsv

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoHeart.narrowPeak -wao|awk 'BEGIN{OFS="\t"; print "Brain_coordinates","Brain_peak"}$NF!=0{print $1,$2,$3,$4}' > Brain-specific-peaks.bed

bedtools intersect -a /tutorial/results/CHIPembryoLiver.narrowPeak -b /tutorial/results/CHIPembryoBrain.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Liver_coordinates","Liver_peak"}$NF==0{print $1":"$2"-"$3,$4}' > Liver-specific-peaks.tsv

bedtools intersect -a /tutorial/results/CHIPembryoLiver.narrowPeak -b /tutorial/results/CHIPembryoBrain.narrowPeak -wao|awk 'BEGIN{OFS="\t"; print "Liver_coordinates","Liver_peak"}$NF!=0{print $1,$2,$3,$4}' > Liver-specific-peaks.bed


BARPLOT
#######

# prepare input

wc -l Brain-specific-peaks.tsv Liver-specific-peaks.tsv Common-peaks-Brain-liver.tsv |awk '$2=="total"{print $2"\t"$1}$2!="total"{split($2,a,"-");print a[1]"\t"$1}' |sed '1iType\tNumber'> input.tsv

nano order.txt

ggbarplot.R -i input.tsv -o common_peaks_brain_liver.pdf --title "H3K4me3 peaks" --y_title "Number of peaks" --x_title "Tissue" --palette_fill /tutorial/palettes/palTissue.txt  --fill_by 1 --header --x_order=order.txt

##############
## BED FILE ##
##############

mkdir ../tss 


## We remove mithocondrial chromosomes 

awk 'BEGIN{FS=OFS="\t"}$1!="chrM" && $3=="transcript" && $0~/gene_type "protein_coding"/ && $7=="+"{ match($0, /transcript_id "([^"]+)/, id); printf ("%s\t20%s\t20%s\t%s\t%s\t\n",$1,$4-200,$4+200,$7,id[1]) }$3=="transcript" && $0~/gene_type "protein_coding"/ && $7=="-"{ match($0, /transcript_id "([^"]+)/, id); printf ("%s\t20%s\t20%s\t%s\t%s\t\n",$1,$5-200,$5+200,$7,id[1]) }' /tutorial/refs/gencode.vM4.gtf |head -n -1 > ./tss/protein-coding-transcript-200up_downTSS.bed

awk 'BEGIN{FS=OFS="\t"}$3=="gene" && $0~/gene_type "protein_coding"/ && $7=="+"{ match($0, /transcript_id "([^"]+)/, id); print $1,$5-200,$5+200,id[1],$7 }$3=="gene" && $0~/gene_type "protein_coding"/ && $7=="-"{ match($0, /transcript_id "([^"]+)/, id); print $1,$5-200,$5+200,id[1],$7 }' /tutorial/refs/gencode.vM4.gtf |sort -u > protein-coding-genes-200up_downTSS.bed

cut -f1 ./analysis/edgeR.over_liver_X_brain.txt |sort -u| grep -Ff - ./tss/protein-coding-genes-200up_downTSS.bed > ./tss/degs-over-liver-TSS.bed

cut -f1 ./analysis/edgeR.over_brain_X_liver.txt |sort -u| grep -Ff - ./tss/protein-coding-genes-200up_downTSS.bed > ./tss/degs-over-brain-TSS.bed

## Intersect
#############

cd tss

## subset

cat Brain-specific-peaks.bed | tail -n+2 | sort -u > brain-header-peaks.bed 

cat Liver-specific-peaks.bed | tail -n+2 | sort -u > liver-header-peaks.bed 

cat Common-peaks-Brain-liver.bed | tail -n+2 | sort -u > common-header-peaks.bed

bedtools intersect -a degs-over-brain-TSS.bed -b brain-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > brain-peaks-over-brain.bed

bedtools intersect -a degs-over-liver-TSS.bed -b liver-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > liver-peaks-over-liver.bed

bedtools intersect -a degs-over-brain-TSS.bed -b liver-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > liver-peaks-over-brain.bed

bedtools intersect -a degs-over-liver-TSS.bed -b brain-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > brain-peaks-over-liver.bed

bedtools intersect -a degs-over-liver-TSS.bed -b common-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > common-peaks-over-liver.bed

bedtools intersect -a degs-over-brain-TSS.bed -b common-header-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > common-peaks-over-brain.bed


















