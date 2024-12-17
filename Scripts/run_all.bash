rm -rf output_files/*

PREF_TO_DATA=../Data/6-Crataegus/crataegus-all
DATA_COALESCENT="${PREF_TO_DATA}/crataegus-all_CDS_concord.no-missing.cf.tree"
DATA_CONCAT="${PREF_TO_DATA}/crataegus-all.244-concat.fsa"

TREE_QMC=../Tools/TREE-QMC/tree-qmc
TREE_QMC_OUT="./output_files/tree_qmc_out"

APRO3=../Tools/ASTER/bin/astral-pro3
APRO3_OUT="./output_files/astral_pro_3_out"

WASTRAL=../Tools/ASTER/bin/wastral
WASTRAL_OUT="./output_files/wastral_out"


IQ_TREE=../Tools/iqtree-2.3.6-Linux-intel/bin/iqtree2
IQ_TREE_OUT="./output_files/iq_tree_2_out"



$TREE_QMC -i $DATA_COALESCENT  -o $TREE_QMC_OUT
$APRO3 -i $DATA_COALESCENT  -o $APRO3_OUT

$IQ_TREE -redo -nt 24 -s $DATA_CONCAT
cp ../Data/6-Crataegus/crataegus-all/crataegus-all.244-concat.fsa.treefile output_files/iq_tree_2_out
