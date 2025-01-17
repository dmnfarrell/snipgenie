"""
    Simulate reads
    Created Sep 2022
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys, os, io, random
import numpy as np
import pandas as pd
import string
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO, Phylo
from . import tools

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__)) #path to module
config_path = os.path.join(home,'.config','snipgenie')
bin_path = os.path.join(config_path, 'binaries')

def run_phastsim(path, ref, newick):
    """Run phastsim """

    refseq = SeqIO.read(ref, 'fasta')
    scale = (10/len(refseq))
    outpath = os.path.join(path, 'phastsim_output')
    cmd = 'phastSim --outpath {o} --outputFile 1 --seed 1 --createFasta' \
             ' --createPhylip --treeFile {n}' \
             ' --scale {s} --invariable .1 --alpha 1.0 --omegaAlpha 1.0' \
             ' --reference {r}'.format(o=outpath,s=scale,r=ref,n=newick)
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    return

def create_meta_data(names):
    """Make fake meta data"""

    data = {}
    years = [2019,2020,2021]
    labels = ['A','B','C']
    for name in names:
        data[name] = {'year':random.choice(years),'label':random.choice(labels)}
    df = pd.DataFrame(data).T
    return df

def simulate_paired_end_reads(fasta_file, read_length, num_reads, output_path):
    """
    Simulates paired-end FASTQ files from a given FASTA file containing genomic sequences.

    Args:
        fasta_file (str): Path to the input FASTA file.
        read_length (int): Length of the simulated reads.
        num_reads (int): Number of reads to generate per sequence, evenly distributed across the genome.
        output_path (str): Directory where the output FASTQ files will be saved.

    Output:
        For each sequence in the FASTA file, generates two gzipped FASTQ files:
        - sample1_1.fastq.gz: Forward reads
        - sample1_2.fastq.gz: Reverse reads
    """

    import gzip
    def generate_quality_scores(length):
        # Generate random quality scores (Phred33 range: '!' to 'I')
        return ''.join(chr(random.randint(33, 73)) for _ in range(length))

    def reverse_complement(sequence):
        # Get the reverse complement of a DNA sequence
        complement = str.maketrans('ACGT', 'TGCA')
        return sequence.translate(complement)[::-1]

    # Ensure output directory exists
    os.makedirs(output_path, exist_ok=True)

    # Convert num_reads to an integer if it's not already
    num_reads = int(num_reads)

    for record in SeqIO.parse(fasta_file, "fasta"):
        sample_name = record.id
        print (sample_name)
        sequence = str(record.seq)
        seq_length = len(sequence)

        if seq_length < read_length:
            print(f"Warning: Sequence in {sample_name} is shorter than the read length ({read_length}). Skipping.")
            continue

        forward_reads = []
        reverse_reads = []

        # Calculate evenly distributed start positions
        step_size = max(1, (seq_length - read_length) // (num_reads - 1)) if num_reads > 1 else 0
        start_positions = [i for i in range(0, seq_length - read_length + 1, step_size)][:num_reads]

        # Generate reads
        for idx, start in enumerate(start_positions):
            forward_read = sequence[start:start+read_length]
            reverse_read = reverse_complement(sequence[start:start+read_length])

            # Generate random quality scores for each read
            forward_quality = generate_quality_scores(read_length)
            reverse_quality = generate_quality_scores(read_length)

            read_id = f"@{sample_name}_{idx+1}"

            # Append paired-end reads
            forward_reads.append(f"{read_id}/1\n{forward_read}\n+\n{forward_quality}\n")
            reverse_reads.append(f"{read_id}/2\n{reverse_read}\n+\n{reverse_quality}\n")

        # Write to gzipped paired-end FASTQ files
        with gzip.open(os.path.join(output_path, f"{sample_name}_1.fastq.gz"), "wt") as f1, \
             gzip.open(os.path.join(output_path, f"{sample_name}_2.fastq.gz"), "wt") as f2:
            f1.writelines(forward_reads)
            f2.writelines(reverse_reads)

def artificial_fastq_generator(ref, outfile, cmp=100):
    """Generate reads from reference"""

    jarfile = os.path.join(bin_path,'ArtificialFastqGenerator.jar')
    cmd = 'java -jar {j} -O {o} -R {r} -S ">temp" -RL 150 -CMP {cmp}'\
        ' -CSD 0.2 -SE true'.format(j=jarfile, r=ref, o=outfile,cmp=cmp)
         #-URQS true -F1 {f1} -F2 {f2}
    print (cmd)
    subprocess.check_output(cmd, shell=True)
    return

def generate_fastqs(infile, outpath, reads=1e5, overwrite=False):
    """Make multiple fastqs"""

    from joblib import Parallel, delayed
    import multiprocessing, time
    num_cores = 4

    simrecs = list(SeqIO.parse(infile,'fasta'))
    def my_func(rec):
        from tempfile import mkstemp
        x,tmp = mkstemp()
        SeqIO.write(SeqRecord(rec.seq,id='temp'), tmp, 'fasta')
        out = os.path.join(outpath,rec.id)
        if os.path.exists(out+'.1.fastq.gz') and overwrite == False:
            print ('found %s' %out)
            return
        cmp = reads*2*150/len(rec.seq)
        print (cmp)
        artificial_fastq_generator(tmp, out, cmp=cmp)
        cmd = 'pigz %s/*.fastq' %outpath
        subprocess.check_output(cmd, shell=True)

    st = time.time()
    #Parallel(n_jobs=num_cores)(delayed(my_func)(i) for i in simrecs)
    #print (time.time()-st)
    for i in simrecs:
        my_func(i)
