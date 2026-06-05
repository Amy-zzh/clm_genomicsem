cd /media/user/Lenovo/GSEM/SMR/smr-1.4.0-linux/smr/

#bifactor
for i in {1..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_BF.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/Brain-mMeta/Brain-mMeta \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_BF_BrainmMeta_smr_multi_chr$i \
--smr-multi
done

#factor 1
for i in {1..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F1.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/Brain-mMeta/Brain-mMeta \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F1_BrainmMeta_smr_multi_chr$i \
--smr-multi
done

#factor 2
for i in {1..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F2.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/Brain-mMeta/Brain-mMeta \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F2_BrainmMeta_smr_multi_chr$i \
--smr-multi
done


#factor 3
for i in {1..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F3.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/Brain-mMeta/Brain-mMeta \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F3_BrainmMeta_smr_multi_chr$i \
--smr-multi
done


#factor 4
for i in {1..22}; do
./smr --bfile /media/user/Lenovo/GSEM/SMR/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
--gwas-summary /media/user/Lenovo/GSEM/SMR/CLMGWAS_SMR_F4.ma \
--chr $i \
--beqtl-summary /media/user/Lenovo/GSEM/SMR/Brain-mMeta/Brain-mMeta \
--out /media/user/Lenovo/GSEM/SMR/Results/CLM_F4_BrainmMeta_smr_multi_chr$i \
--smr-multi
done

