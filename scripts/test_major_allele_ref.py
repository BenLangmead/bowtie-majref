'''
Test the accuracy of a major allele reference
'''
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m', '--ma_ref',
        help='the major allele ref file (fasta) to be tested'
    )
    parser.add_argument(
        '-r', '--ref',
        help='standard ref seq (fasta)'
    )
    parser.add_argument(
        '-c', '--check_points',
        default=None,
        help='sampled check points'
    )
    parser.add_argument(
        '-t', '--target',
        help='target sequence name, e.g. chr1, chr2, .... Leave blank to test all available seqs [None]'
    )
    args = parser.parse_args()
    return args


def read_genome(fn, target):
    '''
    Reads genome and supports more than one chromosomes

    Inputs:
        - fn: FASTA file name
    Output:
        genome in a list
    '''
    genome = {}
    f = open(fn, 'r')
    #: seq[0] is empty to fit vcf coordinate (1-based)
    seq = '^'
    name = None
    for line in f:
        if line.startswith('>'):
            if len(seq) > 1:
                if target == None or name == target:
                    genome[name] = seq
                seq = '^'
            #: update name
            name = line.split()[0][1:]
            if name.startswith('chr') == False:
                name = 'chr' + name
            continue
        seq += line[: line.find('\\')]
    if target == None or name == target:
        genome[name] = seq
    return genome


def print_diff(key, diff_dict, num_tp_dict, fn_dict):
    ''' Summarize differences and print. '''
    fp_dict = diff_dict[key]
    fn_dict = fn_dict[key]
    num_tp = num_tp_dict[key]

    if len(fp_dict) == 0 and len(fn_dict) == 0:
        if num_tp == 0:
            print (key, 'PASS (NO_DIFF)')
        else:
            print (key, 'PASS (NUM_EDITS=%d)' % num_tp)
        return 
    else:
        print (key, 'FAILED')

        print ('fp=', len(fp_dict))
        if len(fp_dict) > 0:
            print(fp_dict)
        print ('tp=', num_tp)
        print ('fn=', len(fn_dict))
        if len(fn_dict) > 0:
            print(fn_dict)


def test_major_allele_ref(fn_ma_ref, fn_ref, fn_check_points, target):
    if target != None and target.startswith('chr') == False:
        target = 'chr' + target

    genome_ref = read_genome(fn_ref, target)
    genome_ma_ref = read_genome(fn_ma_ref, target)

    #: checks keys
    if genome_ref.keys() != genome_ma_ref.keys():
        print ('ERROR: keys are not matched!')
        print (genome_ref.keys())
        print (genome_ma_ref.keys())
        exit()
    print ('SEQNAME TEST PASSED: major and ref have the same set of sequences')

    #: checks lengths of each chromosome
    length_check = True
    for k in genome_ref.keys():
        if len(genome_ref[k]) != len(genome_ma_ref[k]):
            length_check = False
            print ('ERROR: lengths are not matched: {0}'.format(k))
            exit()
    if length_check:
        print ('LENGTH TEST PASSED: lengths of all sequences are equal')
    
    #: stores diff between major and ref
    diff_dict = {}
    for k in genome_ref.keys():
        if (k != target) and (target != None):
            break
        diff_dict[k] = []
        ref_chrom_k = genome_ref[k]
        maj_chrom_k = genome_ma_ref[k]
        for j in range(len(ref_chrom_k)):
            if ref_chrom_k[j] != maj_chrom_k[j]:
                diff_dict[k].append([j, ref_chrom_k[j].upper(), maj_chrom_k[j].upper()])

    #: supports stdin for pipelined commands
    if fn_check_points == None:
        f_cp = sys.stdin
    else:
        f_cp = open(fn_check_points, 'r')

    num_tp_dict = {}
    fn_dict = {}
    for k in genome_ref.keys():
        if (k != target) and (target != None):
            break
        num_tp_dict[k] = 0
        fn_dict[k] = []

    for line in f_cp:
        line = line.split()
        [chrom, vid, pos, ref, alt, af] = line[:]
        pos = int(pos)
        #: case-insentivie comparison, all chars are casted to uppercase
        vcf_ref_allele = ref.upper()
        vcf_alt_allele = alt.upper()
        if chrom.startswith('chr') == False:
            chrom = 'chr' + chrom
        if target != None and chrom != target:
            continue
        diff = diff_dict[chrom]
        var_is_found = False
        for i, d in enumerate(diff):
            if d[0] > pos:
                break
            elif (d[0] == pos) and (d[1] == vcf_ref_allele) and (d[2] == vcf_alt_allele):
                var_is_found = True
                num_tp_dict[chrom] += 1
                diff_dict[chrom].pop(i)
                break
        if var_is_found == False:
            fn_dict[chrom].append([chrom, vid, pos, ref, alt, af])

    if target:
        print_diff(target, diff_dict, num_tp_dict, fn_dict)
    else:
        for k in genome_ref.keys():
            print_diff(k, diff_dict, num_tp_dict, fn_dict)


if __name__ == '__main__':
    args = parse_args()
    fn_ma_ref = args.ma_ref
    fn_ref = args.ref
    fn_check_points = args.check_points
    target = args.target
    test_major_allele_ref(fn_ma_ref, fn_ref, fn_check_points, target)
