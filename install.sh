installationPath=$(pwd)
echo "#!"$installationPath"/src/conda/bin/python" >GRACy.py
echo "installationDirectory = \""$installationPath"/\"" >> GRACy.py
echo " " >> GRACy.py
cat GRACy.py ./src/.GRACy_main.py >temp ; mv temp GRACy.py
chmod +x GRACy.py

if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 already installed"
else

	echo "Installing miniconda 3"
	cd src
	if  test -f "installation.log"; then
	rm installation.log
	fi
	touch installation.log
	bash Miniconda3-latest-Linux-x86_64.sh -b -p ./conda >> installation.log
	cd ../
	if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 successfully installed"
	./src/conda/bin/conda config --set notify_outdated_conda false

	echo "Installing pillow. Please wait...."
	./src/conda/bin/conda install -c anaconda -y pillow >> installation.log
	echo "Checking pillow installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq pillow condaList; then
		echo "Pillow was successfully installed"
	else
		echo "Pillow was not installed. Please check file installation.log for details"
	fi

	echo "Installing numpy. Please wait...."
	./src/conda/bin/conda install -c anaconda -y numpy >> installation.log
	echo "Checking numpy installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq numpy condaList; then
		echo "Numpy was successfully installed"
	else
		echo "Numpy was not installed. Please check file installation.log for details"
	fi

	echo "Installing matplotlib. Please wait...."
	./src/conda/bin/conda install -c conda-forge -y matplotlib >> installation.log
	echo "Checking matplotlib installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq matplotlib condaList; then
		echo "Matplotlib was successfully installed"
	else
		echo "Matplotlib was not installed. Please check file installation.log for details"
	fi

	./src/conda/bin/conda config --add channels bioconda >> installation.log


	echo "Installing Trimgalore. Please wait...."
	./src/conda/bin/conda install -y  trim-galore >> installation.log
	echo "Checking Trimgalore installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq trim-galore condaList; then
		echo "Trimgalore was successfully installed"
	else
		echo "Trimgalore was not installed. Please check file installation.log for details"
	fi

	echo "Installing Bowtie2. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  bowtie2=2.3.5.1 >> installation.log
	./src/conda/bin/conda install -y tbb=2020.2 >> installation.log
	echo "Checking bowtie2 installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq bowtie2 condaList; then
		echo "Bowtie2 was successfully installed"
	else
		echo "Bowtie2 was not installed. Please check file installation.log for details"
	fi

	echo "Installing bam2fastq. Please wait...."
	./src/conda/bin/conda install -c yuxiang -y  bam2fastq=1.1.0  >> installation.log
	echo "Checking bam2fastq installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq bam2fastq condaList; then
		echo "bam2fastq was successfully installed"
	else
		echo "bam2fastq was not installed. Please check file installation.log for details"
	fi

	echo "Installing fastuniq. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  fastuniq=1.1 >> installation.log
	echo "Checking fastuniq installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq fastuniq condaList; then
		echo "fastuniq was successfully installed"
	else
		echo "fastuniq was not installed. Please check file installation.log for details"
	fi

	echo "Installing cutadapt. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  cutadapt=2.6  >> installation.log
	echo "Checking cutadapt installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq cutadapt condaList; then
		echo "cutadapt was successfully installed"
	else
		echo "cutadapt was not installed. Please check file installation.log for details"
	fi

	echo "Installing pypdf2. Please wait...."
	./src/conda/bin/conda install -c conda-forge -y pypdf2=1.26.0 >> installation.log
	echo "Checking pypdf2 installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq pypdf2 condaList; then
		echo "pypdf2 was successfully installed"
	else
		echo "pypdf2 was not installed. Please check file installation.log for details"
	fi

	echo "Installing reportlab. Please wait...."
	./src/conda/bin/conda install -c anaconda -y reportlab=3.5.9 >> installation.log
	echo "Checking reportlab installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq reportlab condaList; then
		echo "reportlab was successfully installed"
	else
		echo "reportlab was not installed. Please check file installation.log for details"
	fi

	echo "Installing biopython. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  biopython=1.76 >> installation.log
	echo "Checking biopython installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq biopython condaList; then
		echo "biopython was successfully installed"
	else
		echo "biopython was not installed. Please check file installation.log for details"
	fi

	echo "Installing jellyfish. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  jellyfish=2.2.10 >> installation.log
	echo "Checking jellyfish installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq jellyfish condaList; then
		echo "jellyfish was successfully installed"
	else
		echo "jellyfish was not installed. Please check file installation.log for details"
	fi

	echo "Installing bwa. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  bwa=0.7.17 >> installation.log
	echo "Checking bwa installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq bwa condaList; then
		echo "bwa was successfully installed"
	else
		echo "bwa was not installed. Please check file installation.log for details"
	fi

	echo "Installing prinseq. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  prinseq=0.20.4 >> installation.log
	echo "Checking prinseq installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq prinseq condaList; then
		echo "prinseq was successfully installed"
	else
		echo "prinseq was not installed. Please check file installation.log for details"
	fi

	echo "Installing khmer. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  khmer=3.0.0  >> installation.log
	echo "Checking khmer installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq khmer condaList; then
		echo "khmer was successfully installed"
	else
		echo "khmer was not installed. Please check file installation.log for details"
	fi

	echo "Installing seqtk. Please wait...." 
	./src/conda/bin/conda install -c bioconda -y  seqtk=1.3 >> installation.log
	echo "Checking seqtk installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq seqtk condaList; then
		echo "seqtk was successfully installed"
	else
		echo "seqtk was not installed. Please check file installation.log for details"
	fi

	echo "Installing spades. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  spades=3.12 >> installation.log
	echo "Checking spades installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq spades condaList; then
		echo "spades was successfully installed"
	else
		echo "spades was not installed. Please check file installation.log for details"
	fi

	echo "Installing picard. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  picard=2.21 >> installation.log
	echo "Checking picard installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq picard condaList; then
		echo "picard was successfully installed"
	else
		echo "picard was not installed. Please check file installation.log for details"
	fi

	echo "Installing lastz. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  lastz=1.0.4  >> installation.log
	echo "Checking lastz installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq lastz condaList; then
		echo "lastz was successfully installed"
	else
		echo "lastz was not installed. Please check file installation.log for details"
	fi


	echo "Installing perl-perl4-corelibs. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  perl-perl4-corelibs  >> installation.log
	echo "Checking perl-perl4-corelibs installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq perl-perl4-corelibs condaList; then
		echo "perl-perl4-corelibs was successfully installed"
	else
		echo "perl-perl4-corelibs was not installed. Please check file installation.log for details"
	fi

	echo "Installing blast. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  blast=2.9.0 >> installation.log
	echo "Checking blast installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq blast condaList; then
		echo "blast was successfully installed"
	else
		echo "blast was not installed. Please check file installation.log for details"
	fi

	echo "Installing cd-hit. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  cd-hit=4.8.1 >> installation.log
	echo "Checking cd-hit installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq cd-hit condaList; then
		echo "cd-hit was successfully installed"
	else
		echo "cd-hit was not installed. Please check file installation.log for details"
	fi

	echo "Installing cap3. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  cap3 >> installation.log
	echo "Checking cap3 installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq cap3 condaList; then
		echo "cap3 was successfully installed"
	else
		echo "cap3 was not installed. Please check file installation.log for details"
	fi

	echo "Installing bedtools. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  bedtools=2.29.2  >> installation.log
	echo "Checking bedtools installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq bedtools condaList; then
		echo "bedtools was successfully installed"
	else
		echo "bedtools was not installed. Please check file installation.log for details"
	fi

	echo "Installing fastx_toolkit. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  fastx_toolkit=0.0.14 >> installation.log
	echo "Checking fastx_toolkit installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq fastx_toolkit condaList; then
		echo "fastx_toolkit was successfully installed"
	else
		echo "fastx_toolkit was not installed. Please check file installation.log for details"
	fi

	echo "Installing blat. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  blat=36 >> installation.log
	echo "Checking blat installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq blat condaList; then
		echo "blat was successfully installed"
	else
		echo "blat was not installed. Please check file installation.log for details"
	fi

	echo "Installing exonerate. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  exonerate=2.4 >> installation.log
	echo "Checking exonerate installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq exonerate condaList; then
		echo "exonerate was successfully installed"
	else
		echo "exonerate was not installed. Please check file installation.log for details"
	fi

	echo "Installing pyqt. Please wait...."
	./src/conda/bin/conda install -c anaconda -y pyqt=5.9.2 >> installation.log
	echo "Checking pyqt installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq pyqt condaList; then
		echo "pyqt was successfully installed"
	else
		echo "pyqt was not installed. Please check file installation.log for details"
	fi

	echo "Installing varscan. Please wait...."
	./src/conda/bin/conda install -c bioconda -y varscan=2.4.4 >> installation.log
	echo "Checking varscan installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq varscan condaList; then
		echo "varscan was successfully installed"
	else
		echo "varscan was not installed. Please check file installation.log for details"
	fi

	echo "Installing tabix. Please wait...."
	./src/conda/bin/conda install -c bioconda -y tabix >> installation.log
	echo "Checking tabix installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq tabix condaList; then
		echo "tabix was successfully installed"
	else
		echo "tabix was not installed. Please check file installation.log for details"
	fi

	echo "Installing samtools. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  samtools=1.3.1 >> installation.log
	echo "Checking samtools installation...."
	./src/conda/bin/conda list > condaList
	if grep -Fq samtools condaList; then
		echo "samtools was successfully installed"
	else
		echo "samtools was not installed. Please check file installation.log for details"
	fi

	else
	 echo "miniconda3 and all its packages were not installed. Please check file installation.log for details"
	fi
fi



if  test -f "./src/conda2/bin/conda"; then
	echo "Miniconda2 already installed"
else


	cd src
	bash Miniconda2-latest-Linux-x86_64.sh -b -p ./conda2 >> installation.log
	cd ../
	if  test -f "./src/conda2/bin/conda"; then
	echo "Miniconda2 successfully installed"

	./src/conda2/bin/conda config --set notify_outdated_conda false

	echo "Installing ragout. Please wait...."
	./src/conda2/bin/conda install -c bioconda -y  ragout=2.2 >> installation.log
	echo "Checking ragout installation...."
	./src/conda2/bin/conda list > condaList
	if grep -Fq ragout condaList; then
		echo "ragout was successfully installed"
	else
		echo "ragout was not installed. Please check file installation.log for details"
	fi


	echo "Installing mummer. Please wait...."
	./src/conda2/bin/conda install -c bioconda -y  mummer >> installation.log
	echo "Checking mummer installation...."
	./src/conda2/bin/conda list > condaList
	if grep -Fq mummer condaList; then
		echo "mummer was successfully installed"
	else
		echo "mummer was not installed. Please check file installation.log for details"
	fi

	echo "Installing bcftools. Please wait...."
	./src/conda2/bin/conda install -c bioconda -y  bcftools >> installation.log
	echo "Checking bcftools installation...."
	./src/conda2/bin/conda list > condaList
	if grep -Fq bcftools condaList; then
		echo "bcftools was successfully installed"
	else
		echo "bcftools was not installed. Please check file installation.log for details"
	fi

	echo "Installing lofreq. Please wait...."
	./src/conda2/bin/conda install -c bioconda -y  lofreq >> installation.log
	echo "Checking lofreq installation...."
	./src/conda2/bin/conda list > condaList
	if grep -Fq lofreq condaList; then
		echo "lofreq was successfully installed"
	else
		echo "lofreq was not installed. Please check file installation.log for details"
	fi

	else
	echo "miniconda2 and all its packages were not installed. Please check file installation.log for details"

	fi
fi

#Check if a c compiler is installed. If not:
# sudo apt-get install build-essential







