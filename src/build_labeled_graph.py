import os, sys
import re
import math
import numpy

from collections import defaultdict

import pysam

from grit.files.reads import RNAseqReads
import multiprocessing
import queue
from grit.lib.multiprocessing_utils import ProcessSafeOPStream, fork_and_wait

import gzip, io

import pickle

DATA_BASE_DIR = os.path.abspath(os.path.dirname(__file__) + "/../data/")

RNASEQ_SORT_INDICES = numpy.array((3,4,1,2,6,7)) 

class TFs(pysam.TabixFile):
    pass

def load_GENCODE_names(fname):
    gene_name_map = defaultdict(list)
    with io.TextIOWrapper(gzip.open(fname, 'rb')) as fp:
        for line in fp:
            if line.startswith("#"): continue
            data = line.split()
            if data[2] != 'gene': continue
            ensemble_id = re.findall("ID=(.*?)\.\d+;", line)[0]
            gene_name = re.findall("gene_name=(.*?);", line)[0]
            contig, start, stop = data[1], int(data[3]), int(data[4])
            gene_name_map[gene_name.upper()].append(ensemble_id)
    return gene_name_map

def load_GENCODE_genes(fname):
    genes = []
    with io.TextIOWrapper(gzip.open(fname, 'rb')) as fp:
        for line in fp:
            if line.startswith("#"): continue
            data = line.split()
            if data[2] != 'gene': continue
            ensemble_id = re.findall("ID=(.*?)\.\d+;", line)[0]
            gene_name = re.findall("gene_name=(.*?);", line)[0]
            contig, start, stop = data[0], int(data[3]), int(data[4])
            genes.append([contig, start, stop, gene_name, ensemble_id])
    
    return genes

def group_overlapping_intervals(intervals, max_size=10000):
    intervals = sorted(intervals)
    if len(intervals) == 0: return []
    
    curr_start, curr_stop = intervals[0][0], intervals[0][1]
    merged_intervals = [([curr_start, curr_stop], [intervals[0],]),]
    for interval in intervals[1:]:
        if curr_stop-curr_start > max_size or interval[0] > curr_stop:
            curr_start, curr_stop = interval[0], interval[1]
            merged_intervals.append(
                ([curr_start, curr_stop], [interval,]) )
        else:
            curr_stop = max(interval[1], curr_stop)
            merged_intervals[-1][0][1] = curr_stop
            merged_intervals[-1][1].append(interval)

    return merged_intervals

def load_enhancers(tfs):
    enhancers = {}
    for contig in tfs.contigs:
        intervals = []
        for contig, start, stop, tf in tfs.fetch(contig):
            intervals.append((int(start), int(stop)))
        merged_intevals = sorted( 
            x[0] for x in group_overlapping_intervals(intervals) )
        enhancers[contig] = merged_intevals
    
    return enhancers

def load_tf_gene_mapping(fname=os.path.join(
        DATA_BASE_DIR, "ENCODE_TFS.target.gene.map.txt")):
    tf_gene_map = defaultdict(list)
    with open(fname) as fp:
        for line in fp:
            data = line.split()
            if len(data) == 1:
                print( data )
                continue
            tf_gene_map[data[0]].extend(data[1].split(","))
    
    gene_tf_map = {}
    for tf, genes in tf_gene_map.items():
        for gene in genes:
            gene_tf_map[gene] = tf
    
    return tf_gene_map, gene_tf_map

def load_tf_sites(fname):
    def extract_pos_and_tfs(line):
        data = line.strip().split()
        chrm, start, stop = data[0], int(data[1]), int(data[2])
        return (chrm, start, stop), data[-1].split(",")
    
    tf_positions = defaultdict(list)
    fp = pysam.TabixFile(fname)
    for contig in fp.contigs:
        for i, record in enumerate(fp.fetch(contig)):
            region, tfs = extract_pos_and_tfs(record)
            for tf in tfs:
                tf_positions[tf].append(region)
            if i > 1000: break
    fp.close()
    return dict(tf_positions)


def load_expression(fname=os.path.join(
        DATA_BASE_DIR, "Het_Project.hg19_mm9_RSEM_gene_expression.txt")):
    header = None
    expression = {}
    with open(fname) as fp:
        for i, line in enumerate(fp):
            if i == 0: 
                header = line.split()[1:]
                continue
            data = line.split()
            gene = data[0].split('.')[0]
            expression[gene] = numpy.array(data[1:], dtype=float)[
                RNASEQ_SORT_INDICES]
    header = [header[i] for i in RNASEQ_SORT_INDICES]
    return header, expression

def load_tads(mouse_fname=os.path.join(
        DATA_BASE_DIR, "./called_TADS/MouseES.HIC.combined.domain.bed"),
              human_fname=os.path.join(
        DATA_BASE_DIR, "./called_TADS/IMR90.HIC.combined.domain.bed")):
              
    tads = defaultdict(set)
    with open(mouse_fname) as fp:
        for line in fp:
            contig, start, stop = line.split()
            tads['mm9_'+contig].add(int(start))
            tads['mm9_'+contig].add(int(stop))

    with open(human_fname) as fp:
        for line in fp:
            contig, start, stop = line.split()
            tads['hg19_'+contig].add(int(start))
            tads['hg19_'+contig].add(int(stop))

    for key, bndries in list(tads.items()):
        tads[key] = numpy.array(sorted(bndries))
        
    return dict(tads)

def load_tf_genes():
    try:
        with open('pickled_genes.obj', 'rb') as fp:
            return pickle.load(fp)
    except FileNotFoundError:
        pass
    
    m4_ann_fname = os.path.join(
        DATA_BASE_DIR, "gencode.vM4.annotation.gff3.gz")
    hg19_ann_fname = os.path.join(
        DATA_BASE_DIR, "gencode.v19.annotation.gff3.gz")
    tf_gene_map, gene_tf_map = load_tf_gene_mapping()

    all_genes = []
    for data in load_GENCODE_genes( hg19_ann_fname ):
        data[0] = 'hg19_' + data[0]
        try: 
            tf = gene_tf_map[data[-1]]
        except KeyError:
            tf = 'NONE'
            continue
        all_genes.append(data + [tf,])
    for data in load_GENCODE_genes( m4_ann_fname ):
        data[0] = 'mm9_' + data[0]
        try: 
            tf = gene_tf_map[data[-1]]
        except KeyError:
            tf = 'NONE'
            continue
        all_genes.append(data + [tf,])

    all_genes = sorted(all_genes)

    with open('pickled_genes.obj', 'wb') as ofp:
        pickle.dump(all_genes, ofp)

    return all_genes

class ATACSeq():
    def __init__(self):
        base = "/data/heterokaryon/ATAC-Seq/wigs/hg19_mm9/"
        all_samples = ["16hr_A", "16hr_rep2", "16hr_rep3", 
                       "3hr_A", "3hr_rep3", 
                       "48hr_A", "48hr_rep2",
                       "48hr_rep3", "CC_rep2 MRC5_rep2"]
        sample_prefixes = [
            '3hr_rep1', '3hr_rep3', 
            '16hr_rep2', '16hr_rep3', 
            '48hr_rep2', '48hr_rep3']

        self.all_signal_coverage = []
        for sample_prefix in sample_prefixes:
            fname = os.path.join(base, sample_prefix) + ".bedgraph.gz"
            self.all_signal_coverage.append(pysam.TabixFile(fname))

    def extract_signal_in_region(self, contig, start, stop):
        rv = []
        for signal_cov in self.all_signal_coverage:
            res = signal_cov.fetch(contig, start, stop)
            rv.append( sum( float(x.split()[-1]) for x in res ) )
        return numpy.array(rv)

    def build_signal_coverage_array(self, contig, start, stop):
        rv = numpy.zeros(
            (len(self.all_signal_coverage), stop-start), dtype=float)
        for i, signal_cov in enumerate(self.all_signal_coverage):
            res = signal_cov.fetch(contig, start, stop)
            for contig, r_start, r_stop, signal in (x.split() for x in res):
                rv[i, int(r_start)-start:int(r_stop)-start] = float(signal)
        return rv

def tf_bs_parser(line):
    print( line )
    data = line.split()
    return (data[0], int(data[1]), int(data[2]), data[3])

def cov_change(exp):
    def calc_zscore(x1s, y1s):
        x1_mu = sum(x1s)/len(x1s)
        x1_var = sum((x - x1_mu)**2 for x in x1s)/len(x1s)
        y1_mu = sum(y1s)/len(y1s)
        y1_var = sum((y - y1_mu)**2 for y in y1s)/len(y1s)
        return (y1_mu - x1_mu)/math.sqrt(
            x1_var/len(x1s) + y1_var/len(y1s) + 1)
        
    z1 = calc_zscore((exp[0], exp[1]), (exp[2], exp[3]))
    z2 = calc_zscore((exp[2], exp[3]), (exp[4], exp[5]))
    return z1, z2

def find_active_enhancers_in_tad(contig, tad_start, tad_stop, 
                                 tfs, tf_genes, all_atacseq,
                                 hg19_enhancers_ofp, mm9_enhancers_ofp):
    local_genes = [gene[4] for gene in tf_genes
                   if gene[0] == contig
                   and not( gene[2] < tad_start or gene[1] > tad_stop )]
    #if len(local_genes) == 0: return
    #print( contig, tad_start, tad_stop, file=sys.stderr )

    local_tfbs = [x.split()for x in tfs.fetch(
        contig, tad_start, tad_stop) ]
    local_tfbs = sorted((int(x[1]), int(x[2]), x[3]) 
                        for x in local_tfbs
                        if int(x[2]) - int(x[1]) < 1000)

    enhancers = group_overlapping_intervals(local_tfbs)

    filtered_enhancers = []
    for enhancer in enhancers:
        e_length = enhancer[0][1]-enhancer[0][0]+1
        cov = all_atacseq.build_signal_coverage_array(
            contig, enhancer[0][0], enhancer[0][1])
        if cov.sum(1).max()/e_length < 1e-3: continue
        #print(cov.sum(1)/e_length)
        #noisy_cov = numpy.random.random(6)/10 + cov.sum(1)/e_length
        z1, z2 = cov_change(cov.sum(1)) # /e_length
        score = max(abs(z1), abs(z2))
        #print( z1, z2, cov.sum(1) ) # /e_length
        #if score < 1:  continue
        #assert False
        filtered_enhancers.append((enhancer, cov, score))

    if len(filtered_enhancers) == 0: return

    if contig.startswith('hg19'):
        new_contig = contig[5:]
        ofp = hg19_enhancers_ofp
    else:
        new_contig = contig[4:]
        ofp = mm9_enhancers_ofp
    for enh, cov, score in filtered_enhancers:
        ofp.write("%s\t%i\t%i\t%s\t%i\t.\n" % (
            new_contig, enh[0][0], enh[0][1], 
            'enhancer', min(1000, int(score*50))))
        for tf_start, tf_stop, name in enh[1]:
            rel_start = tf_start - enh[0][0]
            rel_stop = tf_stop - enh[0][0]

            if cov[:,rel_start:rel_stop].sum() == 0: continue
            z1, z2 = cov_change(cov[:,rel_start:rel_stop].sum(1)) # /e_length
            #print( rel_start, rel_stop, cov.shape, cov[:,rel_start:rel_stop].sum(1) )
            score = max(abs(z1), abs(z2))
            if score < 1: continue
            ofp.write("%s\t%i\t%i\t%s\t%i\t.\n" % (
                new_contig, tf_start, tf_stop, 
                name, min(1000, int(score*50))))

    return

def worker( tads_queue, hg19_enhancers_ofp, mm9_enhancers_ofp):
    tfs = TFs(os.path.join(DATA_BASE_DIR, "ENCODE_TFS.bed.gz"))
    exp_header, expression = load_expression()
    tf_genes = load_tf_genes()
    all_atacseq = ATACSeq()
    initial_size = tads_queue.qsize()
    while tads_queue.qsize() > 0:
        try: 
            contig, tad_start, tad_stop = tads_queue.get(0.1)
        except queue.Empty: 
            break
        print( tads_queue.qsize(), initial_size, file=sys.stderr )
        find_active_enhancers_in_tad(
            contig, tad_start, tad_stop,
            tfs, tf_genes, all_atacseq,
            hg19_enhancers_ofp, mm9_enhancers_ofp)

    os._exit(0)

def main():    
    tads = load_tads()

    tfs = TFs(os.path.join(DATA_BASE_DIR, "ENCODE_TFS.bed.gz"))
    exp_header, expression = load_expression()
    tf_genes = load_tf_genes()
    all_atacseq = ATACSeq()
    
    hg19_enhancers_ofp = ProcessSafeOPStream(open("enhancers.hg19.bed", "w"))
    hg19_enhancers_ofp.write("track type=bed name=hg19_active_enhancers useScore=1\n")
    mm9_enhancers_ofp = ProcessSafeOPStream(open("enhancers.mm9.bed", "w"))
    mm9_enhancers_ofp.write("track type=bed name=mm9_active_enhancers useScore=1\n")

    tads_queue = multiprocessing.Queue()
    for contig, tad_bndrys in sorted(tads.items()):
        for tad_start, tad_stop in zip(tad_bndrys[:-1], tad_bndrys[1:]):
            tads_queue.put((contig, tad_start, tad_stop))
    args = [tads_queue, hg19_enhancers_ofp, mm9_enhancers_ofp]
            
    fork_and_wait(24, worker, args)
    hg19_enhancers_ofp.close()
    mm9_enhancers_ofp.close()

    #print( tads )
    assert False
    
    print("type\tname\tgene_ens_id\t%s" % "\t".join(exp_header))
    for tf_symbol, positions in list(tf_positions.items()):
        if tf_symbol not in gene_name_map:
            print(tf_symbol, file=sys.stderr) #, gene_name_map[tf]
        else:
            for gene_id in gene_name_map[tf_symbol]:
                if gene_id not in expression: 
                    print(("EXP", tf_symbol, gene_id, 
                                          gene_name_map[tf_symbol] ), file=sys.stderr)
                else: 
                    print("%s\t%s\t%s\t%s" % (
                        "GENE_EXP", gene_id, tf_symbol, expression[gene_id]))

if __name__ == '__main__':
    main()

#print bigbed_to_bed("chipseq_track_proximal.bb")
#CS = BigBedFile(open("chipseq_track_proximal.bb"))
#print dir(CS)
#print CS.get("chr1",0,1000000)
