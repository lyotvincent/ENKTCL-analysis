mkdir -p ./data
wget http://secure.bioxai.cn/OEP000498_bulk-processed.tar.gz -O ./data/OEP000498_bulk-processed. tar.gz
tar -xzvf ./data/OEP000498_bulk-processed.tar.gz -C ./data
wget http://secure.bioxai.cn/GSE203663_scRNA-processed.tar.gz -O ./data/GSE203663_scRNA-processed.tar.gz
tar -xzvf ./data/GSE203663_scRNA-processed.tar.gz -C ./data