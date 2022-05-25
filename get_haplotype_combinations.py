from Bio import SeqIO
import collections
import subprocess
import tempfile
import pprint
import argparse

parser = argparse.ArgumentParser(description='Take input arguments')
parser.add_argument('--dir', metavar = 'd', type=str, help='Working directory for analysis outputs', default='.')
parser.add_argument('--haplotypes', metavar='h', type=str, help='Fasta with haplotypes')
parser.add_argument('--marker', metavar='m', type=str, help='Marker name')
parser.add_argument('--indices', metavar='hid', type=str, help='Integers indicating haplotype ID comma-separates')


args = parser.parse_args()


def check_max_overlap(seq1, seq2):
    '''
    Check matches from left and right side of two strings
    '''
    Lc = 0
    Rc = 0
    for ind in range(len(seq1)):
        if seq1[ind] == seq2[ind]:
            Lc = Lc + 1
        else:
            break
    for ind in range(1, len(seq1)):
        if seq1[-ind] == seq2[-ind]:
            Rc = Rc + 1
        else:
            break
    return((Lc, Rc))

# Load args
hap_fasta = args.haplotypes
marker = args.marker
indicies = args.indices.split(',')
string_ints = [marker + "-" + str(id) for id in indicies]
tfile = args.dir + "tmp.txt"

# Create tmp txt file
with open(tfile, "w") as tempfile:
    for element in string_ints:
        tempfile.write(element + "\n")
     
# Run seqkit
log_seqkit = subprocess.call(["seqkit", "grep", "-f" , tfile, hap_fasta, "-o", args.dir + "temp.fasta"])

# Run snp-sites
log_snps = subprocess.call(["snp-sites", "-m" , "-o", args.dir + "temp.snps", args.dir + "temp.fasta"])

fastafile = args.dir + "temp.fasta"

input_file = open(fastafile)
my_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

# Check overlaps
myhashL = collections.defaultdict(dict)
myhashR = collections.defaultdict(dict)

for obj1 in list(my_dict.keys()):
    copy = list(my_dict.keys())
    try:        
        copy.remove(obj1)
    except ValueError:
        pass
    for obj2 in copy:
        seq1 = str(my_dict[obj1].seq)
        seq2 = str(my_dict[obj2].seq)
        left_overlap, right_overlap = check_max_overlap(seq1, seq2)
        myhashL[obj1][obj2] = left_overlap
        myhashL[obj2][obj1] = left_overlap

        myhashR[obj1][obj2] = right_overlap
        myhashR[obj2][obj1] = right_overlap

# Check pairs of combinations
combinations = {}

lenseq = len(list(my_dict.values())[0].seq)
for partition in range(lenseq):
    for obj1 in list(myhashL.keys()):
        potentialL = None
        potentialR = None
        potential=[]
        for key in list(myhashL[obj1].keys()):
            if myhashL[obj1][key] >= partition:
                potentialL = key
            if myhashR[obj1][key] >= lenseq - partition:
                potentialR = key
            if potentialL is not None and potentialR is not None:
                pot = potentialL + " + " + potentialR
                if pot not in potential:
                    potential.append(pot)
            if len(potential) != 0:
                combinations[obj1] = potential

pprint.pprint(combinations)
