samples=$(ls /data/NKTCL/NKTCL_Multi/NKTCL-scRNA)

for s in $samples; do
  echo "copying /data/NKTCL/NKTCL_Multi/NKTCL_EBV_${s}/velocyto/NKTCL_EBV_${s}.loom"
  cp "/data/NKTCL/NKTCL_Multi/NKTCL_EBV_${s}/velocyto/NKTCL_EBV_${s}.loom" "/data/NKTCL/NKTCL_Multi/NKTCL-velo/"
done