'''
Fuctions to process vcf file
'''
import argparse
import sys
import random

class VCF:
    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
    v_chrom = ''
    v_pos = 0
    v_id = ''
    ref_allele = ''
    alt_allele = ''
    v_qual = 0.0
    v_filter = ''
    v_info = ''
    #v_format = ''
    v_type = ''
    v_af = None
    v_line = ''
    
    def __init__(self, row, alt_id, line):
        self.v_line = line
        self.v_chrom = row[0]
        self.v_pos = int(row[1])
        self.v_id = row[2]
        self.ref_allele = row[3]
        assert len(self.ref_allele.split(',')) == 1
        self.alt_allele = row[4].split(',')[alt_id]
        try:
            self.v_qual = float(row[5])
        except:
            self.v_qual = row[5]
        self.v_filter = row[6]
        self.v_info = row[7]
        #self.v_format = row[8]
        
        for info in self.v_info.split(';'):
            if info.startswith('AF='):
                af = info[info.find('=') + 1:]
                self.v_af = float(af.split(',')[alt_id])
                break

        # get var type
        if len(self.alt_allele) == len(self.ref_allele):
            self.v_type = 'SNP'
        # INS
        elif len(self.alt_allele) > len(self.ref_allele):
            self.v_type = 'INS'
        # DEL
        elif len(self.alt_allele) < len(self.ref_allele):
            self.v_type = 'DEL'
        else:
            print ('Error: unexpected variant type :', self.ref_allele, self.alt_allele)
            exit ()
    
    def print(self, show_chrm=True, show_id=True, show_pos=True, show_ref_allele=True, show_alt_allele=True, show_af=True):
        msg = ''
        if show_chrm:
            msg += self.v_chrom
            msg += '\t'
        if show_id:
            msg += self.v_id
            msg += '\t'
        if show_pos:
            msg += str(self.v_pos)
            msg += '\t'
        if show_ref_allele:
            msg += self.ref_allele
            msg += '\t'
        if show_alt_allele:
            msg += str(self.alt_allele)
            msg += '\t'
        if show_af:
            msg += str(self.v_af)
        print (msg)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        default=None,
        help='vcf file; leave blank for std.in [None]'
    )
    parser.add_argument(
        '--out_var_loc',
        default=None,
        help='output file storing filtered variant locations; no output if this argument is not specified [None]'
    )
    parser.add_argument(
        '--out_no_conflicting_vcf',
        default=None,
        help='output vcf file with conflicting variants removed; no output if this argument is not specified [None]'
    )
    parser.add_argument(
        '--min_af', type=float, default=0.0,
        help='[--out_var_loc and --out_no_conflicting_vcf] min allele frequency to be considered [0.0]'
    )
    parser.add_argument(
        '--indels', type=int, default=0,
        help='[--out_var_loc and --out_no_conflicting_vcf] set this to 1 to consider indels [0]'
    )
    parser.add_argument(
        '--mnps', type=int, default=0,
        help='[--out_var_loc and --out_no_conflicting_vcf] set this to 1 to consider mnps [0]'
    )
    parser.add_argument(
        '--rand_th', type=float, default=1.0,
        help='[--out_var_loc] fraction of variants should be kept, by default keeping all variants [1.0]'
    )
    parser.add_argument(
        '--max_allele_len', type=int, default=-1,
        help='[--out_var_loc and --out_no_conflicting_vcf] filter out variants with either ref or alt allele longer than this number; set a negative integer to allow all variants [-1]'
    )
    args = parser.parse_args()
    return args

def build_vcf(line, consider_indels, consider_mult, max_allele_len):
    if line.startswith('##'):
        return []
    #: header
    if line.startswith('#'):
        return []
    row = line.split()
    num_alt_alleles = len(row[4].split(','))
    #: skip if doesn't consider multi-allelic locus
    if consider_mult == 0 and num_alt_alleles > 1:
        return []
    list_vcf = []
    for i in range(num_alt_alleles):
        v = VCF(row, i, line)
        if consider_indels == 0 and v.v_type in ['INS', 'DEL']:
            continue
        #: skip max_allele_len check if it's a negative number
        if max_allele_len > 0:
            #: ignore when either ref_allele or alt_allele is longer than max_allele_len
            if len(v.ref_allele) > max_allele_len or len(v.alt_allele) > max_allele_len:
                continue
        list_vcf.append(v)
    return list_vcf

def comp_var_with_list(var, list_var):
    list_comp = [None] * len(list_var)
    for list_idx, v in enumerate(list_var):
        offset = var.v_pos - v.v_pos
        #: check if ref allele is the same
        ref_allele_check = True 
        for idx, ref_a in enumerate(var.ref_allele):
            try:
                if ref_a != v.ref_allele[offset + idx]:
                    ref_allele_check = False
            except:
                break
        assert ref_allele_check == True
        #: check if alt allele is the same
        alt_allele_check = True
        for idx, alt_a in enumerate(var.alt_allele):
            try:
                if alt_a != v.alt_allele[offset + idx]:
                    alt_allele_check = False
            except:
                break
        list_comp[list_idx] = alt_allele_check
    return list_comp

def specify_target_var(
        vcf_fn, 
        out_var_loc_fn, 
        min_af, 
        consider_indels, 
        consider_mnps, 
        rand_th,
        max_allele_len):
    list_pos = []
    list_ref_allele = []
    list_v = []
    if vcf_fn != None:
        vcf_f = open(vcf_fn, 'r')
    else:
        vcf_f = sys.stdin
    max_pos = 0
    import time
    for line in vcf_f:
        list_vcf = build_vcf(line, consider_indels, consider_mnps, max_allele_len)
        #: a header or an ignored variant
        if len(list_vcf) < 1:
            continue
        for vcf in list_vcf:
            ##
            ## if vcf.v_type != 'SNP':
            ##     print (line[:60])
            ##     time.sleep(0.5)
            ##
            if vcf.v_af <= min_af:
                continue
            if vcf.v_pos in list_pos:
                continue
            if random.random() <= rand_th:
                vcf.print()
            if vcf.v_pos > max_pos:
                max_pos = vcf.v_pos

def remove_conflicting_vars(vcf_fn, out_vcf_fn, min_af, consider_indels, consider_mnps, max_allele_len):
    vcf_f = open(vcf_fn, 'r')
    prev_range = range(0)
    prev_var = []
    for line in vcf_f:
        if line.startswith('#'):
            continue
        list_vcf = build_vcf(line, consider_indels, consider_mnps, max_allele_len)
        for vcf in list_vcf:
            #: ignore var with allele freq less than given threshold
            #: keep var if allele freq is not specified in vcf
            if vcf.v_af != None and vcf.v_af < min_af:
                continue
            if vcf.v_pos in prev_range:
                comp = comp_var_with_list(vcf, prev_var)
                for c in comp:
                    if c == False:
                        for pv in prev_var:
                            pv.print()
                        vcf.print()
                        print (comp)
                prev_var.append(vcf)
            if len(prev_range) == 0 or vcf.v_pos > max(prev_range):
                prev_range = range(vcf.v_pos, vcf.v_pos + len(vcf.ref_allele))
                prev_var = [vcf]


if __name__ == '__main__':
    args = parse_args()
    vcf_fn = args.vcf
    out_vcf_fn = args.out_no_conflicting_vcf
    out_var_loc_fn = args.out_var_loc
    min_af = args.min_af
    consider_indels = args.indels
    consider_mnps = args.mnps
    rand_th = args.rand_th
    max_allele_len = args.max_allele_len
    if out_vcf_fn != None:
        sys.stderr.write('Mode: remove conflicting variants\n')
        sys.stderr.write('Multi-allelic loci are ignored\n')
        #: current: report conflicts, no output yet
        remove_conflicting_vars(vcf_fn, out_vcf_fn, min_af, consider_indels, consider_mnps, max_allele_len)
    elif out_var_loc_fn != None:
        sys.stderr.write('Mode: specify target variant locations\n')
        sys.stderr.write('Multi-allelic loci are ignored\n')
        specify_target_var(vcf_fn, out_var_loc_fn, min_af, consider_indels, consider_mnps, rand_th, max_allele_len)
