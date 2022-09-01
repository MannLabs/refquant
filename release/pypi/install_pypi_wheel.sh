conda create -n refquant_pip_test python=3.8 -y
conda activate refquant_pip_test
pip install "refquant[stable]"
refquant
conda deactivate
