# Apollo: A Sequencing-Technology-Independent, Scalable, and Accurate Assembly Polishing Algorithm

Apollo is an assembly polishing algorithm that attempts to correct the errors in an assembly. It can take multiple set of reads in a single run and polish the assemblies of genomes of any size.

## Installing Apollo

* Make sure you have a compiler that has support for C++14.
* Download the code from its GitHub repository.

```bash
git clone https://github.com/CMU-SAFARI/Apollo.git apollo
```
*  Change directory to `./apollo` and run the Makefile. If everything goes well, you will have a binary called `apollo` inside the `bin` folder.

```bash
cd ./apollo
make
cd ./bin
```
Now you can copy this binary wherever you want (preferably under a directory that is included in your `$PATH`). Assuming that you are in the directory that the binary is located, you may run the command below to display the help message.

```bash
./apollo -h
```

## Assembly polishing
Polishing using a single set of reads (i.e., non-hybrid):

Assume that you have 1) an assembly `assembly.fasta`, 2) a set of reads `reads.fasta`, 3) the alignment file `alignment.bam` that contains the alignment of the reads to the assembly, 4) and you would like to store polished assembly as `polished.fasta`. The command below uses `30` threads while polishing the assembly:

```bash
./apollo -a assembly.fasta -r reads.fasta -m alignment.bam -t 30 -o polished.fasta
```
Resulting fasta file `polished.fasta` will be the final output of Apollo.

Polishing using a hybrid set of reads:

Assume that you have 1) an assembly `assembly.fasta`, 2) a hybrid set of reads `reads1.fasta` and `reads2.fasta`, 3) the alignment of these reads to the assembly stored in `alignment1.bam` and `alignment2.bam`, respectively, 4) and you would like to store polished assembly as `polished.fasta`. The command below uses `30` threads while polishing the assembly:

```bash
./apollo -a assembly.fasta -r reads1.fasta -r reads2.fasta -m alignment1.bam -m alignment2.bam -t 30 -o polished.fasta
```
Resulting fasta file `polished.fasta` will be the final output of Apollo.

## Supported and Required Input Files
### Alignment File

* Apollo supports alignment files in BAM format. If you have a SAM file you can easily convert your `input.sam` to `input.bam` using the following command:

```bash
samtools view -hb input.sam > input.bam
```
* Apollo requires the input BAM file to be coordinate sorted. You can sort your `input.bam` file using the following command:

```bash
samtools view -h -F4 input.bam | samtools sort -m 16G -l0 > input_sorted.bam
```

* Apollo needs the BAM file to be indexed. You can index your `input.bam` file using the following command:

```bash
samtools index input.bam
```

### Set of Reads

* Apollo supports the reads set in FASTA format. For each read (i.e., sequence), the number of characters per line has to be the same, except for the last line. For example, a sequence of length 1000 can either be represented in a single line with 1000 characters or can be split into multiple lines where each line include the equal number of characters. Only exception here is the last line, which can have any number of characters but no more than the characters that the prior lines have. An illustration of a sequence with a length of 10 would be:

>\>read1  
>TAT  
>TAT  
>ATT  
>A

or in a single line:

>\>read1  
TATTATATTA

The restriction on the number of characters per line is required as Apollo constructs the index file (i.e., FAI file) for the input read set. Further information about indexing and the requirements can be found at: https://seqan.readthedocs.io/en/master/Tutorial/InputOutput/IndexedFastaIO.html

* If there are *too* long reads in a input read set, we recommend dividing these reads into smaller chunks to reduce the memory requirements. Apollo supports chunking during run time. One can simply use command to divide the reads into chunks of size 1000 (maximum) while polishing:

```bash
./apollo -a assembly.fasta -r reads1.fasta -r reads2.fasta -m alignment1.bam -m alignment2.bam -t 30 -o polished.fasta -c 1000
```

## Example run

You may use the following test run to check whether everything works as intended with Apollo. Note that you must have `curl` to download the required files and also `minimap2` to map the reads to the assembly.

```bash
#create a test folder
mkdir test; cd test
#download a read set that is publicly available by PacBio and only fetch small number of read set as this is a sanity check
curl -s http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/additional_data/2590969/0002/Analysis_Results/m140928_104939_ethan_c100699582550000001823139903261541_s1_p0.3.subreads.fasta | head -5000 > pacbio.fasta
#download the already constructed assembly
curl -L -o assembly.fasta http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/polished_assembly/polished_assembly.fasta
#generate read-to-assembly file
minimap2 -x map-pb -a assembly.fasta pacbio.fasta | samtools view -h -F4 | samtools sort -m 16G -l0 > alignment.bam
#indexing the alignment file
samtools index alignment.bam
#polishing. Here we assume that "apollo" is in your $PATH. If not you should specify the exact path to "apollo"
apollo -a assembly.fasta -r pacbio.fasta -m alignment.bam -o polished.fasta -c 1000
```

## Problems You May Encounter

### Input Format
* Apollo currently does not support reads in a compressed format such as `input.fasta.gz`. These FASTA files must be uncompressed. Future release may support compressed files as well.

* Apollo currently does not support paired-end reads. Those paired-end reads can be provided as multiple input read sets to the Apollo where they also should have multiple read-to-assembly alignment files. Another option is to merge the paired-end Illumina reads into one FASTA file but one should make sure that the read ids per sequence are unique.

* Apollo requires `samtools` to be preinstalled in order to generate the index files and sort the alignment file.

### Publication and citing Apollo

If you would like to cite Apollo, please cite the following publication:

> Can Firtina, Jeremie S. Kim, Mohammed Alser, Damla Senol Cali, A. Ercument Cicek, Can Alkan, and Onur Mutlu,
> "Apollo: A Sequencing-Technology-Independent, Scalable, and Accurate Assembly Polishing Algorithm"
> *Bioinformatics*, **btaa179**, 2020. [doi:10.1093/bioinformatics/btaa179][doi]
