import yassi
import random
import time
from Bio.Alphabet import IUPAC
from Bio import Motif
from Bio.Seq import Seq

def random_site(n):
    return "".join([random.choice("acgt") for i in range(n)])

def reverse_complement(seq):
    base_complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T':'A',
                        'a': 't', 'c': 'g', 'g': 'c', 't':'a'}
    complement = ''.join(base_complements[l] for l in seq)
    return complement[::-1]

def get_genome():
    gbk_file = "/home/sefa/Dropbox/projects/comparative_genomics/data/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna"
    with open(gbk_file) as f:
        genome = "".join([line.strip() for line in f.readlines()[1:]])
    return "".join(g for g in genome if g in "ATGC")

def get_motif():
    motif_file = "/home/sefa/Dropbox/projects/comparative_genomics/data/LexA.sites"
    with open(motif_file) as f:
        lines = [line.strip() for line in f.readlines() if not line.startswith('>')]
    return lines

def biopython_search(motif, genome):
    """Biopython site search test"""
    # biopython doesn't like lowercase sequences, make them uppercase
    m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
    for s in motif:
        m.add_instance(Seq(s.upper(), m.alphabet))
    #print m.log_odds()
    start_time = time.time()
    x = sorted(list(m.search_pwm(Seq(genome.upper(), m.alphabet), both=False)), key=lambda x: x[1], reverse=True)
    for i in x[:5]:
        print i
    end_time = time.time()
    print "biopython search took %s secs." % (end_time - start_time)
    
def yassi_search(motif, genome):
    """Yassi site search test"""
    start_time = time.time()
    x = yassi.search(motif, genome)
    for i in x[:5]:
        print i
    end_time = time.time()
    print "yassi search took %s secs."  % (end_time - start_time)

if __name__ == '__main__':
    # example usage:
    motif = [random_site(8) for i in range(10)]
    genome = random_site(50000)
    #yassi.search(motif, genome)
    yassi.build_PSSM(motif)
    
    #motif = get_motif()
    #genome = get_genome()
    #genome = genome + reverse_complement(genome)
    #biopython_search(motif, genome)
    #yassi_search(motif, genome)
