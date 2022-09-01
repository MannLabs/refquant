conda create -n refquant python=3.8 -y
conda activate refquant
pip install -e '../.[stable,development-stable]'
refquant
conda deactivate
