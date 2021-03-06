NanoSim 1.0.0

------------------------------------------------------------------------------------
NanoSim is a fast and scalable read simulator that captures the technology-specific 
features of ONT data, and allows for adjustments upon improvement of nanopore
sequencing technology.
------------------------------------------------------------------------------------

------------------------------------------------------------------------------------
Dependencies:

LAST (Tested with version 581)
R (Tested with version 3.2.3)
Python (2.6 or above)
Numpy (Tested with version 1.10.1 or above)

------------------------------------------------------------------------------------
Usage

NanoSim is implemented using R for error model fitting and Python for read length 
analysis and simulation. The first step of NanoSim is read characterization, which 
provides a comprehensive alignment-based analysis, and generates a set of read
profiles serving as the input to the next step, the simulation stage. The simulation 
tool uses the model built in the previous step to produce in silico reads for a
given reference genome. It also outputs a list of introduced errors, consisting of 
the position on each read, error type and reference bases.

1. Characterization stage

Characterization stage takes a reference and a training read set in FASTA format as 
input. User can also provide their own alignment file in MAF format.

Usage:

./read_analysis.py <options>  
    [options]:  
    -h : print usage message  
    -i : training ONT real reads, must be fasta files  
    -r : reference genome of the training reads  
    -m : User can provide their own alignment file, with maf extension. Optional  
    -o : The prefix of output file, default = 'training'  

2. Simulation stage

Simulation stage takes reference genome and read profiles as input and outputs
simulated reads in FASTA fomat.

Usage:

./simulator.py [command] <options>  
   [command]:  
    circular | linear  
    # Do not choose 'circular' when there is more than one sequence in the reference  
    <options>:  
    -h : print usage message
    -r : reference genome in fasta file, specify path and file name. Required  
    -c : the prefix of training set profiles, same as the output prefix in 
         read_analysis.py, default = training  
    -o : The prefix of output file, default = 'simulated'  
    -n : Number of generated reads, default = 20,000 reads  
    --perfect: Output perfect reads, no mutations. Optional 
    --KmerBias: prohibits homopolymers with length >= 6 bases in output reads. Optional 

For example:
1 If you want to simulate E. coli genome, then circular command must be chosen 
  because it's a circular genome
./simulator.py circular -r Ecoli_ref.fasta -c ecoli

2 If you want to simulate only perfect reads, i.e. no snps, or indels, just simulate 
  the read length distribution
./simulator.py circular -r Ecoli_ref.fasta -c ecoli --perfect

3 If you want to simulate S. cerevisiae genome with kmer bias, then linear command 
  must be chosen because it's a linear genome
./simulator.py linear -r yeast_ref.fasta -c yeast --KmerBias

See more detailed example in example.sh

------------------------------------------------------------------------------------
Explaination of output files

1. Characterization stage
    training_aligned_length_ecdf: Length distribution of aligned regions on aligned 
                                  reads
    training_aligned_reads_ecdf: Length distribution of aligned reads
    training_align_ratio: Empirical distribution of align ratio of each read
    training_besthit.maf: The best alignment of each read based on length
    training_match.hist/training_mis.hist/training_del.hist/training_ins.hist: 
        Histogram of match, mismatch, and indels
    training_first_match.hist: Histogram of the first match length of each alignment
    training_error_markov_model: Markov model of error types
    training_ht_ratio: Empirical distribution of the head region vs total unaligned 
                       region
    training.maf: The output of LAST, alignment file in MAF format
    training_match_markov_model: Markov model of the length of matches (stretches 
                                 of correct base calls)
    training_model_profile: Fitted model for errors
    training_processed.maf: A re-formatted MAF file for user-provided alignment file
    training_unaligned_length_ecdf: Length distribution of unaligned reads

2. Simulation stage
    simulated.log: Log file for simulation process

    simulated_reads.fasta: FASTA file of simulated reads. Each reads has "unaligned",
        "aligned", or "perfect" in the header determining their error rate. 
        "unaligned" means that the reads have an error rate over 90% and cannot be 
        aligned. "aligned" reads have the same error rate as training reads. "perfect"
        reads have no errors.

    To explain the information in the header, we have two examples:
    a. >ref|NC-001137|-[chromosome=V]_468529_unaligned_0_F_0_3236_0
    All information before the first _ are chromosome information. 468529 is the start
    position and unaligned suggesting it should be unaligned to the reference. The 
    first 0 is the sequence index. F represents a forward strand. 0_3236_0 means that
    sequence length extracted from the reference is 3236 bases.
        
    b. >ref|NC-001143|-[chromosome=XI]_115406_aligned_16565_R_92_12710_2
    This is an aligned read coming from chromosome XI at position 115406. 16565 is 
    the sequence index. R represents a reverse complement strand. 92_12710_2 
    means that this read has 92-base head region (cannot be aligned), followed by 
    12710 bases of middle region, and then 2-base tail region.

    The information in the header can help users to locate the read easily.

    simulated_error_profile: Contains all the information of errors introduced into
    each reads, including error type, position, original bases and current bases.
