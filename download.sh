mkdir -p ./data
mkdir -p ./figures
mkdir -p ./tables
wget http://secure.bioxai.cn/ENKTCL/gene_pos.txt -O ./data/gene_pos.txt
wget http://secure.bioxai.cn/ENKTCL/cellphonedb.zip -O ./data/cellphonedb.ziptar.gz
wget http://secure.bioxai.cn/ENKTCL/NKTCL_data-processed.tar.gz -O ./data/NKTCL_data-processed.tar.gz
tar -xzvf ./data/NKTCL_data-processed.tar.gz -C ./data
wget http://secure.bioxai.cn/ENKTCL/OEP000498_bulk-processed.tar.gz -O ./data/OEP000498_bulk-processed.tar.gz
tar -xzvf ./data/OEP000498_bulk-processed.tar.gz -C ./data

