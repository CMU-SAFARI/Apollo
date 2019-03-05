# Apollo: A Sequencing-Technology-Independent, Scalable, and Accurate Assembly Polishing Algorithm

## Installing Apollo

* Make sure you have a compiler that has support for C++14.
* Download the code from its GitHub repository.

```bash
git clone https://github.com/CMU-SAFARI/Apollo.git
```
*  Change directory to `Apollo/src/` and run the Makefile. If everything goes well, you will have a binary called `apollo` inside the `bin` folder.

```bash
cd Apollo/src/
make
cd ../bin/
```
Now you can copy this binary wherever you want (preferably under a directory that is included in your `$PATH`). Assuming that you are in the directory that the binary is located, you may run the command below to display the help message.

```bash
apollo -h
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
