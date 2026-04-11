mkdir -p ./data
mkdir -p ./figures
mkdir -p ./tables
wget http://secure.bioxai.cn/gene_pos.txt -O ./data/gene_pos.txt
wget http://secure.bioxai.cn/cellphonedb.zip -O ./data/cellphonedb.zip
wget http://secure.bioxai.cn/NKTCL_data-processed.tar.gz -O ./data/NKTCL_data-processed.tar.gz
tar -xzvf ./data/NKTCL_data-processed.tar.gz -C ./data

