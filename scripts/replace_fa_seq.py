'''
Replace a sequence in a fasta file with new sequence
'''
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--fasta_fn',
        help='target fasta file'
    )
    parser.add_argument(
        '-s', '--seq_name',
        help='name of the sequence to be replaced'
    )
    parser.add_argument(
        '-n', '--newfasta_fn',
        help='new sequence'
    )
    args = parser.parse_args()
    return args

def replace_fa_seq(fasta_fn, seq_name, newfasta_fn):
    target_name = ''
    target_seq = ''
    is_target = False
    fasta_f = open(fasta_fn, 'r')
    for line in fasta_f:
        if line.startswith('>'+seq_name+' '):
            target_name = line.rstrip()
            is_target = True
        elif line.startswith('>'):
            is_target = False
            print (line.rstrip())
        elif is_target:
            target_seq += line.rstrip()
        else:
            print (line.rstrip())
    print (target_name + 'LO:1')

    newfasta_f = open(newfasta_fn, 'r')
    newseq = ''
    for line in newfasta_f:
        if line.startswith('>') == False:
            line = line.rstrip()
            print (line)
            newseq += line

    sys.stderr.write('len(newseq) = {0}, len(target_seq) = {1}\n'.format(len(newseq), len(target_seq)))
    assert len(newseq) == len(target_seq)
    sys.stderr.write('Successfully replace {0}, length={1}\n'.format(seq_name, len(newseq)))

if __name__ == '__main__':
    args = parse_args()
    fasta_fn = args.fasta_fn
    seq_name = args.seq_name
    newfasta_fn = args.newfasta_fn
    replace_fa_seq(fasta_fn, seq_name, newfasta_fn)
