#--------- cd-hit:

git clone https://github.com/weizhongli/cdhit.git
cd cdhit
make

# ------- ANGSD: 

# install xz first from https://tukaani.org/xz/

cd
wget https://tukaani.org/xz/xz-5.2.4.tar.gz --no-check-certificate
tar vxf xz-5.2.4.tar.gz 
cd xz-5.2.4/
./configure --prefix=$HOME/xz-5.2.4/
make
make install

# edit .bashrc:
cd
nano .bashrc
   export LD_LIBRARY_PATH=$HOME/xz-5.2.4/lib:$LD_LIBRARY_PATH
   export LIBRARY_PATH=$HOME/xz-5.2.4/lib:$LIBRARY_PATH
   export C_INCLUDE_PATH=$HOME/xz-5.2.4/include:$C_INCLUDE_PATH
logout
# re-login

# now, install htslib:
cd
git clone https://github.com/samtools/htslib.git
cd htslib
make CFLAGS=" -g -Wall -O2 -D_GNU_SOURCE -I$HOME/xz-5.2.4/include"

# install ANGSD
cd
git clone https://github.com/ANGSD/angsd.git 
cd angsd
make HTSSRC=../htslib

# now adding ANGSD to $PATH
cd
nano .bashrc
# section 2:
   export PATH=$HOME/angsd:$PATH
   export PATH=$HOME/angsd/misc:$PATH
# save (Ctl-O, Ctl-X)

# -------  NGSadmix :

cd ~/bin/
wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
cd -

# ------- Moments: 
cd
git clone https://bitbucket.org/simongravel/moments.git 
cd moments
python setup.py build_ext --inplace

# add this to .bashrc, section 2:
  export PYTHONPATH=$PYTHONPATH:$HOME/moments
# re-login

# to see if it worked:
python
import moments
# if you get an error message something is wrong, if you just see >>> it is all fine
# quit()

# ------- PDGspider

cd ~/bin
wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.1.1.5.zip
unzip PGDSpider_2.1.1.5.zip

# ----- Bayescan :

cd ~/bin
wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
unzip BayeScan2.1.zip
cp BayeScan2.1/binaries/BayeScan2.1_linux64bits bayescan
chmod +x bayescan
rm -r BayeScan*

#-------  stairwayPlot 

# project page: https://sites.google.com/site/jpopgen/stairway-plot
cdw
# get version from June 2016 (v2beta2)
wget https://www.dropbox.com/s/toxnlvk8rhe1p5h/stairway_plot_v2beta2.zip
unzip stairway_plot_v2beta2.zip
mv stairway_plot_v2beta2 stairway_plot_v2beta

#------  EEMS (this one takes quite some time)

# Download and unpack boost from http://www.boost.org:
cd
wget https://dl.bintray.com/boostorg/release/1.70.0/source/boost_1_70_0.tar.gz
tar vxf boost_1_70_0.tar.gz

cd boost_1_70_0
./bootstrap.sh --prefix=$HOME/boost_1_70_0/libraries
./b2 install

# download and unpack Eigen from http://eigen.tuxfamily.org
cd
wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
tar vxf 3.3.7.tar.gz
mv eigen-eigen-323c052e1731 eigen

# git clone eems
cd
git clone https://github.com/dipetkov/eems.git

# compile eems
cd eems/runeems_snps/src
nano Makefile 
# edit these lines to match YOUR Eigen and boost locations!
EIGEN_INC=/home1/01211/cmonstr/eigen
BOOST_INC=/home1/01211/cmonstr/boost_1_70_0/libraries/include
BOOST_LIB=/home1/01211/cmonstr/boost_1_70_0/libraries/lib
# ctl-X 
make linux

cp runeems_snps ~/bin

# does it work?
cds
runeems_snps

#-------------- sanity checks, to see if everything works now
# re-login

cd
bowtie2
angsd
cutadapt

# if you get "command not found" something is wrong 
