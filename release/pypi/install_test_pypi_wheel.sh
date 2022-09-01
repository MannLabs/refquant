conda create -n refquant_pip_test python=3.8 -y
conda activate refquant_pip_test
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple "refquant[stable]"
refquant
conda deactivate
