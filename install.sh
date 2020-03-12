installationPath=$(pwd)
echo "#!"$installationPath"/src/conda/bin/python" >GRACy.py
echo "installationDirectory = \""$installationPath"/\"" >> GRACy.py
echo " " >> GRACy.py
cat GRACy.py ./src/.GRACy_main.py >temp ; mv temp GRACy.py
chmod +x GRACy.py

if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 already installed"
else


	cd src
	bash Miniconda3-latest-Linux-x86_64.sh -b -p ./conda
	cd ../
	if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 successfully installed"
	fi
fi

./src/conda/bin/conda install -c anaconda -y pillow
./src/conda/bin/conda install -c anaconda -y numpy


