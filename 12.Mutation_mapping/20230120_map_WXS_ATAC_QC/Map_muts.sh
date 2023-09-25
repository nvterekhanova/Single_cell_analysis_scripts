#2023-01-20: check UCEC/CRC pairs.
#2022-07-01: do mapping using WXS-mutations. Need to provide piece_ids for each cancer.
#2022-06-26: based on /diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/Map_germline_v/CPT0023690004/

dir='/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/13.Mut_mapping/Objects.v.7/20230120_map_WXS_ATAC_QC'

###CRC:
disease=CRC
for sample in CM618C1-T1Y2 CM618C2-S1Y2 CM663C1-S1Y1 CM663C1-T1Y1 CM1563C1-S1Y1 CM1563C1-T1Y1 CM268C1-S1 CM268C1-T1

do 
perl /diskmnt/Projects/HTAN_analysis/snATAC/Tools/10Xmapping/10Xmapping.pl --bam $dir/inputs/$disease/BAMs/$sample.possorted_bam.bam --maf $dir/inputs/$disease/$disease\_PanCanAllNonSilent.maf --out $dir/out/$disease/muts.mapped/$sample/$sample.step1.out

perl /diskmnt/Projects/HTAN_analysis/snATAC/Tools/10Xmapping/parse_scrna_bc.pl $dir/out/$disease/muts.mapped/$sample/$sample.step1.out $dir/out/$disease/muts.mapped/$sample/$sample.step2.out

done


###Save lists of barcodes for each sample:
for sample in CM618C1-T1Y2 CM618C2-S1Y2 CM663C1-S1Y1 CM663C1-T1Y1 CM1563C1-S1Y1 CM1563C1-T1Y1 CM268C1-S1 CM268C1-T1

do
perl create_matrix.20220702.pl $sample $disease

done


###UCEC:
disease=UCEC

for sample in CPT4096DU-T1 CPT4096DU-S1 CPT4427DU-T1 CPT4427DU-S1 CPT1541DU-T1 CPT1541DU-S1 CPT2373DU-S1 CPT2373DU-T1 CPT704DU-S1 CPT704DU-T1

do 
perl /diskmnt/Projects/HTAN_analysis/snATAC/Tools/10Xmapping/10Xmapping.pl --bam $dir/inputs/$disease/BAMs/$sample.possorted_bam.bam --maf $dir/inputs/$disease/$disease\_PanCanAllNonSilent.maf --out $dir/out/$disease/muts.mapped/$sample/$sample.step1.out

perl /diskmnt/Projects/HTAN_analysis/snATAC/Tools/10Xmapping/parse_scrna_bc.pl $dir/out/$disease/muts.mapped/$sample/$sample.step1.out $dir/out/$disease/muts.mapped/$sample/$sample.step2.out

done


###Save lists of barcodes for each sample:
for sample in CPT4096DU-T1 CPT4096DU-S1 CPT4427DU-T1 CPT4427DU-S1 CPT1541DU-T1 CPT1541DU-S1 CPT2373DU-S1 CPT2373DU-T1 CPT704DU-S1 CPT704DU-T1

do
perl create_matrix.20220702.pl $sample $disease

done

