#!/usr/bin/env python

'''Should take in a BAM+VCF from RNA-Seq along with a
VCF consisting of DNA mutations and compare them.
Returns the % of sites in the DNA VCF with sufficient
depth in the RNA data that have a matching genotype'''

import sys, signal, time
import itertools
import pysam
import vcf

RNA_DEPTH_THRESHOLD = 20
DNA_DEPTH_THRESHOLD = 10
NUM_ROWS_SUMMARY = 4096

class MyVCF(object):
    '''A wrapper around vcf that requires the annotations to be in sorted order'''
    def __init__(self, vcf_in):
        self.vcf_in = vcf_in
        self.cur = 1
        self.cur_chrom_num = None
        self.cur_pos = None
        self.next()

    def next(self):
        if self.cur is None:
            return None
        try:
            while True:
                cur = next(self.vcf_in)
                chrom = cur.CHROM
                if not chrom.startswith('chr') or not chrom[3:].isdigit():
                    continue
                self.cur = cur
                self.cur_chrom_num = int(chrom[3:])
                self.cur_pos = cur.POS
                return cur
        except StopIteration:
            self.cur = None
            return None

    def fetch(self, chrom_num, pos):
        if self.cur is None:
            return None
        while chrom_num > self.cur_chrom_num:
            self.next()
            if self.cur is None: return None
        if chrom_num < self.cur_chrom_num:
            # We already passed this, so it must not've existed
            return None
        # We should now be on the correct chromosome
        while chrom_num == self.cur_chrom_num and pos > self.cur_pos:
            self.next()
            if self.cur is None: return None
        if chrom_num != self.cur_chrom_num or pos < self.cur_pos:
            # We passed the requested site
            return None
        # We should be at the correct site
        if chrom_num != self.cur_chrom_num or pos != self.cur_pos:
            sys.stderr.write(
                "ERROR: found incorrect position for {}:{}".format(chrom_num, pos))
            return None
        return self.cur

def get_VCF_depth(sample):
    try:
        return sample['DP']
    except AttributeError:
        return None
    
def get_depth(bam, chrom, pos):
    depth = 0
    for read in bam.fetch(chrom, pos-1, pos):
        overlap = read.overlap(pos-1, pos)
        if overlap == 1:
            depth += 1
        elif overlap > 1:
            sys.stderr.write("ERROR: read.overlap > 1 = {}".format(overlap))
    return depth

def main(rna_bam, rna_vcf, dna_vcf):
    dna_vcf_in = vcf.Reader(open(dna_vcf, 'r'))
    rna_vcf_in = MyVCF(vcf.Reader(open(rna_vcf, 'r')))
    bam_in = pysam.Samfile(rna_bam, 'rb')

    # Initialize results
    sample_names = []
    sample_matches = None
    sample_overlap = []
    alt_matches = []
    alt_overlap = []
    num_overlap = 0
    total = 0
    alt_allele_matches = []
    alt_allele_overlap = []

    # Stats about types of match/mismatch: compares RNA_DNA
    not_autosomal = 0
    not_simple_SNP = 0
    no_record = 0
    ref_ref = 0
    alt_alt = 0
    ref_alt = 0
    alt_ref = 0

    def print_summary():
        print num_overlap, "/", total, time.strftime("%c")
        for i, name, count, overlap, alt_count, alt_over, alt_allele_match, alt_allele_over in zip(itertools.count(1), sample_names, sample_matches, sample_overlap, alt_matches, alt_overlap, alt_allele_matches, alt_allele_overlap):
            
            print i, name, count, overlap, alt_count, alt_over, alt_allele_match, alt_allele_over,

            if overlap == 0: overlap = 1   # Protect against divide by zero
            if alt_over == 0: alt_over = 1 
            if alt_allele_over == 0: alt_allele_over = 1 # Protect against divide by zero
            
            print "{:.2%}".format(count * 1.0 / overlap), "{:.2%}".format(alt_count * 1.0 / alt_over), "{:.2%}".format(alt_allele_match * 1.0 / alt_allele_over), i

    def print_stats():
        print "Stats: {} {} {} {} {} {} {}".format(not_autosomal, not_simple_SNP, no_record, ref_ref, alt_alt, ref_alt, alt_ref)
    
    # Set timer to shutdown
    def alarm_handler(signum, frame):
        print_summary()
        print_stats()
        print "Execution timed out at {}:{}".format(record.CHROM, record.POS)
        sys.exit(5)
    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(172 * 60 * 60)   # parameter in seconds
    
    for record in dna_vcf_in:
        chrom_prefix = True
        try:
            chrom_num = int(record.CHROM)
            chrom_prefix = False
        except ValueError:
            chrom_num = None
        if (not (chrom_num is not None and chrom_num >=1 and chrom_num <= 22)) and (not record.CHROM.startswith('chr') or not record.CHROM[3:].isdigit()):
            # Not a numbered chromosome -> IGNORE
            not_autosomal += 1
            continue
        # Skip non-SNPs and SNPs with more than 2 alt alleles
        if record.REF[1:] != '' or record.ALT[1:] != [] or record.ALT[0].sequence[1:] != '':
            not_simple_SNP += 1
            continue
        # First time initialization
        if sample_matches is None:
            sample_matches = []
            for sample in record.samples:
                sample_names.append(sample.sample)
                sample_matches.append(0)
                sample_overlap.append(0)
                alt_matches.append(0)
                alt_overlap.append(0)
                alt_allele_matches.append(0)
                alt_allele_overlap.append(0)
        total += 1
        # Check for it in VCF
        if chrom_num is None:
            chrom_num = int(record.CHROM[3:])
        pos = record.POS
        rna_record = rna_vcf_in.fetch(chrom_num, pos)
        if rna_record is None:
            no_record += 1
            # Doesn't appear in VCF, check the BAM for depth
            if chrom_prefix:
                depth = get_depth(bam_in, record.CHROM, pos)
            else:
                depth = get_depth(bam_in, 'chr'+record.CHROM, pos)
            # If we have a certain depth, then register this as a mismatch, otherwise ignore
            if depth > RNA_DEPTH_THRESHOLD:
                num_overlap += 1
                for i, sample in enumerate(record.samples):
                    sdp = get_VCF_depth(sample)
                    if sdp is not None and sdp <= DNA_DEPTH_THRESHOLD:
                        continue
                    sgt = sample['GT']
                    if sgt is None or '.' in sgt:
                        # Don't count overlap when no call or only partial call
                        continue
                    sample_overlap[i] += 1
                    if sgt == '0/0':
                        sample_matches[i] += 1
                        ref_ref += 1
                        # None of 4 compared alleles are alt -> no alt_allele_overlap
                    elif sgt == '0/1' or sgt == '1/0':
                        alt_overlap[i] += 1
                        ref_alt += 1
                        alt_allele_overlap[i] += 1
                        #if i == 2: sys.stderr.write('Mismatch at {}:{} RNA=0/0 DNA={}\n'.format(chrom_num, pos, sgt))
                    elif sgt == '1/1':
                        alt_overlap[i] += 1
                        ref_alt += 1
                        alt_allele_overlap[i] += 2
                        #if i == 2: sys.stderr.write('Mismatch at {}:{} RNA=0/0 DNA={}\n'.format(chrom_num, pos, sgt))
        else:
            # Position appears in VCF, see who the genotype matches
            if rna_record.INFO['DP'] <= RNA_DEPTH_THRESHOLD:
                # only count overlap when RNA sample has sufficient depth
                continue
            num_overlap += 1

            # make sure REF and ALT alleles match, before considering it a match
            if rna_record.REF != record.REF:
                sys.stderr.write("Warning: disagreement about REF allele at {}:{} RNA={} DNA={}\n".format(chrom_num, pos, rna_record.REF, record.REF))
                continue
            if rna_record.ALT[0].sequence != record.ALT[0].sequence:
                # Count as sample match for samples which called this site with depth
                for i, sample in enumerate(record.samples):
                    # FIXME: this could make the overlap slightly too large if multiple records at this pos
                    sgt = sample['GT']
                    sdp = get_VCF_depth(sample)
                    if sgt is not None and (sgt == '1/1' or sgt == '0/1' or sgt == '1/0') and (sdp is None or sdp > DNA_DEPTH_THRESHOLD):
                        sample_overlap[i] += 1
                        alt_overlap[i] += 1
                continue # Go to next record
            
            # Now we know the alt alleles match
            gt = rna_record.samples[0]['GT']
            # Skip if site is uncalled or half-called
            if gt is None or '.' in gt:
                continue
            rna_alt_allele_count = gt.count('1')
            for i, sample in enumerate(record.samples):
                sgt = sample['GT']
                sdp = get_VCF_depth(sample)
                if sgt is None or '.' in sgt or (sdp is not None and sdp <= DNA_DEPTH_THRESHOLD):
                    continue
                if sgt == '0/1' or sgt == '1/0':
                    # Alt allele for DNA
                    sample_overlap[i] += 1
                    alt_overlap[i] += 1
                    if sgt == gt:
                        sample_matches[i] += 1
                        alt_matches[i] += 1
                        alt_alt += 1
                    else:
                        ref_alt += 1
                    alt_allele_matches[i] += 1
                    alt_allele_overlap[i] += rna_alt_allele_count
                elif sgt == '1/1':
                    # Alt allele for DNA
                    sample_overlap[i] += 1
                    alt_overlap[i] += 1
                    if sgt == gt:
                        sample_matches[i] += 1
                        alt_matches[i] += 1
                        alt_alt += 1
                    else:
                        ref_alt += 1
                    alt_allele_matches[i] += rna_alt_allele_count
                    alt_allele_overlap[i] += 2
                elif sgt == '0/0':
                    # Ref allele for DNA
                    sample_overlap[i] += 1
                    if sgt == gt:
                        sample_matches[i] += 1
                        ref_ref += 1
                    else:
                        alt_ref += 1
                    alt_allele_overlap[i] += rna_alt_allele_count
                        
    # Final summary
    signal.alarm(0)  # Disable timeout
    if sample_matches is None:
        sample_matches = []
    print_summary()
    print_stats()
        
if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.stderr.write('Usage: correlate_SNPs.py rna_BAM rna_VCF dna_VCF\n')
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
    print "done"
    
