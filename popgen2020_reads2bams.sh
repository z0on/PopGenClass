Acropora millepora: Orpheus-Keppel comparison

READS (put them ito your $SCRATCH/RAD):
/corral-repl/utexas/tagmap/matz_shared/OK/OK_concatenated.tgz

GENOME (put it into your $WORK/db)::
/corral-repl/utexas/tagmap/matz_shared/OK/amilV2_chroms.fasta

BAMS (just in case, dont download them yet): 
/corral-repl/utexas/tagmap/matz_shared/OK/bams
SAFs, ibsMat, OK_dadi.data (just in case, dont download them yet):
/corral-repl/utexas/tagmap/matz_shared/OK/

BEFORE STARTING, replace, in this whole file:
	- matz@utexas.edu by your actual email;
	- yourusername with your TACC user name.

The idea is to copy the chunks separated by empty lines below and paste them into your cluster 
terminal window consecutively. 

The lines beginning with hash marks (#) are explanations and additional instructions - 
please make sure to read them before copy-pasting. 

==============================================

# login to TACC
ssh yourusername@ls5.tacc.utexas.edu


# ------------  INSTALLATIONS:

# downloading and installing  2bRAD scripts in $HOME/bin 
cd
mkdir bin 
cd ~/bin 
# cloning github repositories
git clone https://github.com/z0on/2bRAD_denovo.git
# move scripts to ~/bin from sub-directories
mv 2bRAD_denovo/* . 
# remove now-empty directories
rm -rf 2bRAD_denovo 

# designating all .pl, .R, and .py files (perl, R, and python scripts) as executable
chmod +x *.pl 
chmod +x *.py
chmod +x *.R
cd

# ----- configuring your environment

# adding ~/bin to your $PATH; loading gmodules
cd
nano .bashrc

# paste this where appropriate 
   export PATH=$HOME/bin:$PATH

# and this, into another section
	module load Rstats
	module load python
	module load samtools
	module load bedtools
	module load bowtie
	module load cutadapt
# press ctl-O, Enter, ctl-X


#==============================


# downloading the data (A.millepora from Orpheus and Keppels)
# switching to SCRATCH:
cds
# making a directory to work in, and going in there:
mkdir RAD
cd RAD

# downloading the gzipped data from ranch storage server
cp /corral-repl/utexas/tagmap/matz_shared/OK/OK_concatenated.tgz .
# press Ctl-Z (pause process), then say 
bg %
# to send the copying process to the background. Proceed with the walkthrough while the data are copied.

# ONLY IF if the previous chunk did not work ("permission denied" error):
scp cmonstr@ranch.tacc.utexas.edu:/samfs/fs1/01211/cmonstr/OK_concatenated.tgz .
# call Misha to punch in password and token
# press Ctl-Z (pause process), then say 
bg %
# to send the copying process to the background. Proceed with the walkthrough while the data are copied.

# downloading and installing 2bRAD scripts in $HOME/bin
cd
mkdir bin # make bin directory
cd bin # go there
# cloning github repository
git clone https://github.com/z0on/2bRAD_denovo.git
mv 2bRAD_denovo/* . # move it to here from sub-directory
rm -rf 2bRAD_denovo # remove now-empty directory

# designating all .pl and .py files (perl and python scripts) as executable
chmod +x *.pl 
chmod +x *.py

# adding ~/bin to your $PATH (plus invoking a bunch of pre-loaded modules)
cd
nano .bashrc
	#paste this in the appropriate section:
module load java
module load Rstats
module load cd-hit
module load python
module load bioperl
module load samtools
module load bowtie
module load cutadapt 

export PATH="~/bin:$PATH"

	# press ctl-O, Enter, ctl-X

# re-login

# does it work?
# try running:
2bRAD_trim_launch.pl
# if you get "command not found" something is wrong

# switch to where the data are:
cds
cd RAD
tar vxf OK_concatenated.tgz &

# how many  *.fq files we have? (should be 68)
ll *.fq | wc -l
# the command above is a typical linux one-liner script: is first lists all files with filenames ending with .fq, and then, instead of printing it all to screen as it normally would, it pipes it (|) into a command wc that, given option -l, counts lines => the output is the number of files ending with .fq

#=================== Pre-processing reads (removing adaptors and quality-trimming)

# how do the reads look before trimming (displaying top 25 DNA sequences in O9.fq):
head -100 K11.fq | grep -E "^[ATGCN]+$"
# in the line above, we meet grep, one of the most useful linux commands. It finds matches in text files and prints the matching lines. Patterns can be just some word, or a complicated "regular expression" like here: "^[ATGCN]+$", which matches lines composed of entirely A,T,G,C, or N symbols in any combination. Here is your cheat sheet for regular expressions: http://jkorpela.fi/perl/regexp.html 

# removing adaptor and trimming bad quality end-bases using cutadapt
# minimum length after trimming = 25, quality cutoff = PHRED Q15 on both ends.

# Here is a loop that writes a list of commands to run (same command for all files) to new file called filt:
# (this one is for reference-based analysis!)
>filt
for file in *.fq; do 
echo "cutadapt -a AGATCGGA --format fastq -q 15,15 -m 25 -o ${file/.fq/}.trim $file > ${file}_trimlog.txt" >> filt;
done

# creating job script based on commands in filt:
ls5_launcher_creator.py -j filt -n filt -t 0:15:00 -a mega2014 -e matz@utexas.edu -w 48 -N 3
# submitting job :
sbatch filt.slurm

# how is our job doing?
squeue -u yourusername

# Done! do we have the right number of output files?
ll *.fq | wc -l # this is how many fastq and fq files we started with
ll *.trim | wc -l  # this is how many trimmed files we got, should be the same number

# >>> DO IT YOURSELF: how do the reads look after trimming? (modify the 'head -100 ...' command above...)

# how many reads remained after trimming?
# let's pick one sample, say O9:
# how many reads were there originally?
grep @HWI O9.fq | wc -l
# >>> how many reads remained?

# now the read pre-processing is all done! time for genotyping.

# =========================== genome placement

# switching to $WORK, moving genome here

cds
mkdir db
cd db
#cdw db/
cp /corral-repl/utexas/tagmap/matz_shared/OK/amilV2_chroms.fasta .

# ONLY IF the above line did not work ("permission denied" error), run this:
# scp cmonstr@ls5.tacc.utexas.edu:/corral-repl/utexas/tagmap/matz_shared/OK/amilV2_chroms.fasta .
# call Misha to punch in password

# creating shortcut to the genome file
export GENOME_FASTA=$SCRATCH/db/amilV2_chroms.fasta

# indexing genome for bowtie2 mapper
echo "bowtie2-build $GENOME_FASTA $GENOME_FASTA" >btb
ls5_launcher_creator.py -j btb -n btb -l btbl -t 0:30:00 -a mega2014 -e matz@utexas.edu -w 1
sbatch btbl

# samtools index (this one is fast, can run on login node)
samtools faidx $GENOME_FASTA

#------------ Mapping and compressing into bam files
cds
cd RAD
export GENOME_FASTA=$SCRATCH/db/amilV2_chroms.fasta

# map with bowtie2 with soft-clipping (to avoid indel artifacts near read ends) and modified settings for better mapping of short 2bRAD tags
>maps2
for file in *.trim; do 
echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $GENOME_FASTA -U $file -S ${file/.trim/}.sam && \
samtools sort -O bam -o ${file/.trim/}.bam ${file/.trim/}.sam && samtools index ${file/.trim/}.bam " >> maps2;
done

ls5_launcher_creator.py -j maps2 -n maps2 -t 6:00:00 -w 24 -a mega2014 -e matz@utexas.edu -q normal
sbatch maps.slurm

# what are those sam files? 
less -S K10.sam
# (mapping quality - 5th column - uniqueness of the match in genome, the higher, the more unique)
# (mutations, if any, are listed in MD:Z: field  - scroll right)

cat maps.e*
# you will see a series of statements like this one:
5293541 reads; of these:
  5293541 (100.00%) were unpaired; of these:
    1056274 (19.95%) aligned 0 times
    1641138 (31.00%) aligned exactly 1 time
    2596129 (49.04%) aligned >1 times
80.05% overall alignment rate

# sanity check: number of bam files should be the same number as number of sam files 
ls *bam | wc -l  

# done! ready for fun stuff.





