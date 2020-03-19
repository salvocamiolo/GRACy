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



if  test -f "./src/conda2/bin/conda"; then
	echo "Miniconda2 already installed"
else


	cd src
	bash Miniconda2-latest-Linux-x86_64.sh -b -p ./conda2
	cd ../
	if  test -f "./src/conda2/bin/conda"; then
	echo "Miniconda2 successfully installed"
	fi
fi

#Check if a c compiler is installed. If not:
# sudo apt-get install build-essential

./src/conda/bin/conda install -c anaconda -y pillow
./src/conda/bin/conda install -c anaconda -y numpy
./src/conda/bin/conda install -c conda-forge -y matplotlib 
./src/conda/bin/conda config --add channels bioconda
./src/conda/bin/conda install -y  trim-galore
./src/conda/bin/conda install -c bioconda -y  bowtie2=2.3.5.1
./src/conda/bin/conda install -c bioconda -y  samtools=1.3.1
./src/conda/bin/conda install -c yuxiang -y  bam2fastq=1.1.0
./src/conda/bin/conda install -c bioconda -y  fastuniq=1.1
./src/conda/bin/conda install -c bioconda -y  cutadapt=2.6
./src/conda/bin/conda install -c conda-forge -y pypdf2=1.26.0
./src/conda/bin/conda install -c anaconda -y reportlab=3.5.9
./src/conda/bin/conda install -c bioconda -y  biopython=1.76
./src/conda/bin/conda install -c bioconda -y  jellyfish=2.2.10
./src/conda/bin/conda install -c bioconda -y  bwa=0.7.17
./src/conda/bin/conda install -c bioconda -y  prinseq=0.20.4
./src/conda/bin/conda install -c bioconda -y  khmer=3.0.0
./src/conda/bin/conda install -c bioconda -y  seqtk=1.3
./src/conda/bin/conda install -c bioconda -y  spades=3.12
./src/conda/bin/conda install -c bioconda -y  picard=2.21
./src/conda/bin/conda install -c bioconda -y  lastz=1.0.4
./src/conda2/bin/conda install -c bioconda -y  ragout=2.2
./src/conda/bin/conda install -c bioconda -y  perl-perl4-corelibs
./src/conda/bin/conda install -c bioconda -y  blast=2.9.0
./src/conda/bin/conda install -c bioconda -y  lofreq=2.1.4
./src/conda/bin/conda install -c bioconda -y  cd-hit=4.8.1
./src/conda/bin/conda install -c bioconda -y  cap3
./src/conda/bin/conda install -c bioconda -y  bedtools=2.29.2
./src/conda/bin/conda install -c bioconda -y  fastx_toolkit=0.0.14
./src/conda/bin/conda install -c bioconda -y  blat=36
./src/conda2/bin/conda install -c bioconda -y  mummer
./src/conda/bin/conda install -c bioconda -y  bcftools=1.8
./src/conda/bin/conda install -c bioconda -y  exonerate=2.4
./src/conda/bin/conda install -c anaconda -y pyqt=5.9.2




