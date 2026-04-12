mkdir -p ./numbat_results
wget http://secure.bioxai.cn/NKTCL_allele.tar.gz -O ./numbat_results/NKTCL_allele.tar.gz
wget http://secure.bioxai.cn/NKTCL_numbat.tar.gz -O ./numbat_results/NKTCL_numbat.tar.gz
tar -xzvf ./data/NKTCL_numbat.tar.gz -C ./numbat_results/
tar -xzvf ./data/NKTCL_allele.tar.gz -C ./numbat_results/

