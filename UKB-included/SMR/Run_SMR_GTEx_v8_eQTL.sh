
#cat -A /media/user/Lenovo/GSEM/SMR/tissues.txt
#dos2unix /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {1..2}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {3..4}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {5..6}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {7..8}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {9..10}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {11..12}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {13..14}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {15..16}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt


cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {17..18}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {19..20}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#bifactor
for i in {21..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt


#cat -A /media/user/Lenovo/GSEM/SMR/tissues.txt
#dos2unix /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {1..2}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {3..4}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {5..6}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {7..8}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {9..10}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {11..12}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {13..14}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {15..16}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt


cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {17..18}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {19..20}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F1
for i in {21..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt



#cat -A /media/user/Lenovo/GSEM/SMR/tissues.txt
#dos2unix /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {1..2}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {3..4}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {5..6}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {7..8}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {9..10}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {11..12}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {13..14}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {15..16}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt


cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {17..18}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {19..20}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F2
for i in {21..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt



#cat -A /media/user/Lenovo/GSEM/SMR/tissues.txt
#dos2unix /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {1..2}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {3..4}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {5..6}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {7..8}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {9..10}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {11..12}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {13..14}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {15..16}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt


cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {17..18}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {19..20}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F3
for i in {21..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt




#cat -A /media/user/Lenovo/GSEM/SMR/tissues.txt
#dos2unix /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {1..2}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {3..4}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {5..6}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {7..8}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {9..10}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {11..12}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {13..14}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {15..16}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt


cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {17..18}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {19..20}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt

cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/
while read tissue; do
#F4
for i in {21..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/GTEx_V8_cis_eqtl/eQTL/$tissue.lite \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_GTEx_eqtl_smr_multi_$tissue.chr$i \
--smr-multi
done
done < /media/user/Lenovo/GSEM/SMR/tissues.txt


