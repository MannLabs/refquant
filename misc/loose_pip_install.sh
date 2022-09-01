conda create -n refquant python=3.8 -y
conda activate refquant
pip install -e '../.[development]'
refquant
conda deactivate
