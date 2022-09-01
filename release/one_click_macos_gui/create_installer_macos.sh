#!bash

# Initial cleanup
rm -rf dist
rm -rf build
FILE=refquant.pkg
if test -f "$FILE"; then
  rm refquant.pkg
fi
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n refquantinstaller python=3.8 -y
conda activate refquantinstaller

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_macos_gui
pip install "../../dist/refquant-0.0.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.2
pyinstaller ../pyinstaller/refquant.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../refquant/data/*.fasta dist/refquant/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/refquant/Contents/Resources
cp ../logos/alpha_logo.icns dist/refquant/Contents/Resources
mv dist/refquant_gui dist/refquant/Contents/MacOS
cp Info.plist dist/refquant/Contents
cp refquant_terminal dist/refquant/Contents/MacOS
cp ../../LICENSE.txt Resources/LICENSE.txt
cp ../logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/refquant --identifier de.mpg.biochem.refquant.app --version 0.0.1 --install-location /Applications/refquant.app --scripts scripts refquant.pkg
productbuild --distribution distribution.xml --resources Resources --package-path refquant.pkg dist/refquant_gui_installer_macos.pkg
