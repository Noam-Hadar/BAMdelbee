from os import system
import sys

d = int(sys.argv[1])

seq = ''.join(open('fasta/IL10RA.fasta', 'r').read().split('\n')[1:])

def readFormat(coordinate, sequence):
    return '\n'.join(['@' + coordinate,sequence, '+', 'F' * len(sequence) + '\n'])

m = 0
for p in range(150 , len(seq)-150):
    mini_seq = seq[p:p+300]
    for n in range(150):
        m += 1
        del_seq = mini_seq[:n] + mini_seq[n + d:]
        del_seq = readFormat(str(m), del_seq[:150])
        open('del_' + str(d) + '.fastq', 'a').write(del_seq)
        ins_seq = mini_seq[:n] + 'A' * d + mini_seq[n:]
        ins_seq = readFormat(str(m), ins_seq[:150])
        open('ins_' + str(d) + '.fastq', 'a').write(ins_seq)

system('gzip del_' + str(d) + '.fastq')
system('gzip ins_' + str(d) + '.fastq')
    
