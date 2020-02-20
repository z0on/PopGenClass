# INSTALLATIONS

------- ANGSD: 

# install xz first from https://tukaani.org/xz/
cd
wget https://tukaani.org/xz/xz-5.2.3.tar.gz --no-check-certificate
tar vxf xz-5.2.3.tar.gz 
cd xz-5.2.3/
./configure --prefix=$HOME/xz-5.2.3/
make
make install

# edit .bashrc:
nano .bashrc
   export LD_LIBRARY_PATH=$HOME/xz-5.2.3/lib:$LD_LIBRARY_PATH
   export LIBRARY_PATH=$HOME/xz-5.2.3/lib:$LIBRARY_PATH
   export C_INCLUDE_PATH=$HOME/xz-5.2.3/include:$C_INCLUDE_PATH
logout
# re-login

# now, install htslib:
cd
git clone https://github.com/samtools/htslib.git
cd htslib
make CFLAGS=" -g -Wall -O2 -D_GNU_SOURCE -I$HOME/xz-5.2.3/include"

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

-------  NGSadmix :
cd ~/bin/
wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
cd -

-------  ngsRelate :
cd 
git clone https://github.com/ANGSD/NgsRelate.git
cd NgsRelate
make HTSSRC=../htslib
cp ngs* ~/bin/
cd

------ ngsLD :

cd 
git clone https://github.com/fgvieira/ngsLD.git
cd ngsLD

nano Makefile
add -I${TACC_GSL_INC}  to CC and CXX macros (CFLAGS= ...);
and -L${TACC_GSL_LIB} to the 'LIB = ...' line.

module load gsl
export PKG_CONFIG_PATH=/opt/apps/intel18/gsl/2.2.1/lib/pkgconfig/
make
cp ngsLD ~/bin

-------  stairwayPlot :

# project page: https://sites.google.com/site/jpopgen/stairway-plot
cdw
# get version from June 2016 (v2beta2)
wget https://www.dropbox.com/s/toxnlvk8rhe1p5h/stairway_plot_v2beta2.zip
unzip stairway_plot_v2beta2.zip
mv stairway_plot_v2beta2 stairway_plot_v2beta


ls *.bam > bams


#===================== A  N  G  S  D =====================

# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).
# if your coverage is >10x, go to GATK section below

# install ANGSD (see InstallingSoftware.txt file... this can be a challenge, so let me know when/if you get stuck)

#----------- assessing base qualities and coverage depth

# entering interactive session, giving all node's memory to one process:
idev -tpn 1 -N 1

FILTERS="-uniqueOnly 1 -minMapQ 20 -maxDepth 10000"
# if only looking at high-confidence  SNPs
# FILTERS="-uniqueOnly 1 -minMapQ 20 -maxDepth 10000 -snp_pval 1e-5"

# T O   D O : 
 TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
# if only looking at high-confidence SNPs:
# TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2 -doMajorMinor 1 -doMaf 1"

# in the following line, -r argument is ~1 Mb (no need to do this for whole genome)
# (look up lengths of your contigs in the header of *.sam files if you need)
angsd -b bams -r chr10:1-1000000 -GL 1 $FILTERS $TODO -P 12 -out dd

# summarizing results (using modified script by Matteo Fumagalli)
Rscript ~/bin/plotQC.R dd
cat dd.info 
# scp dd.pdf to laptop to see distribution of base quality scores and fraction of sites in each sample depending on coverage threshold


#--------------- population structure (based on common polymorphisms, allele freq >0.05)

FILTERS="-uniqueOnly 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 54 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"

# Starting angsd with -P the number of parallel processesors. 
echo "angsd -b bams -GL 1 $FILTERS $TODO -P 12 -out OKall" >aa
ls5_launcher_creator.py -j aa -n aa -t 0:30:00 -a tagmap -e youremail@utexas.edu -w 1
sbatch aa.slurm 

# analyze OKall.ibsMat using OKall_ibs.R to identify clones to remove

# RERUNNING angsd after removal of clones:

# scp bams_noclones from laptop to here

FILTERS="-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 52 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2"
echo "angsd -b bams_noclones -GL 1 $FILTERS $TODO -P 12 -out OK" >ab
ls5_launcher_creator.py -j ab -n ab -t 0:30:00 -a tagmap -e youremail@utexas.edu -w 1
sbatch ab.slurm

#----------------------  ngs ADMIXTURE

# (must install NGSadmix)

for K in 2 3 4; 
do 
NGSadmix -likes OK.beagle.gz -K $K -P 10 -o ok_k${K};
done

# making a table of bams : population correspondence
cat bams_noclones | perl -pe 's/(.)(.+)/$1$2\t$1/' >inds2pops
 
 #==================================
# scp *Mat, *qopt, inds2pops files to laptop, use OK_ibs.R to plot PCA and admixturePlotting_v5.R to plot qopt
#==================================
