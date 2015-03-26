import os, sys
import re

import numpy

from collections import defaultdict

import pysam

from idr.idr import load_bed, merge_peaks

import gzip, io


DATA_BASE_DIR = os.path.abspath(os.path.dirname(__file__) + "/../data/")

def load_GENCODE_names(fname):
    gene_name_map = defaultdict(list)
    with open(fname) as fp:
        for line in fp:
            if line.startswith("#"): continue
            data = line.split()
            if data[2] != 'gene': continue
            ensemble_id = re.findall("ID=(.*?)\.\d+;", line)[0]
            gene_name = re.findall("gene_name=(.*?);", line)[0]
            gene_name_map[gene_name.upper()].append(ensemble_id)
    return gene_name_map

def load_enhancers(tf_fnames):
    # we build enhancers from overlapping intervals of TF binding
    def extract_pos_and_tfs(line):
        data = line.strip().split()
        chrm, start, stop = data[0], int(data[1]), int(data[2])
        return (chrm, start, stop), set(x.split("(")[0] for x in data[-1].split(","))
    
    tf_positions = defaultdict(list)
    for fname in fnames:
        if fname.endswith("gz"): 
            fp = io.TextIOWrapper(gzip.open(fname, 'rb'))
        else: 
            fp = open(fname)

        for i, line in enumerate(fp):
            region, tfs = extract_pos_and_tfs(line)
            for tf in tfs:
                tf_positions[tf].append(region)
    
    return dict(tf_positions)

def load_tf_sites(fname):
    def extract_pos_and_tfs(line):
        data = line.strip().split()
        chrm, start, stop = data[0], int(data[1]), int(data[2])
        return ((chrm, start, stop), 
                set(x.split("(")[0] for x in data[-1].split(",")))
    
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


def load_expression(fname="all_quant.txt"):
    header = None
    expression = {}
    with open(fname) as fp:
        for i, line in enumerate(fp):
            if i == 0: 
                header = line.split()[1:]
                continue
            data = line.split()
            gene = data[0].split('.')[0]
            expression[gene] = "\t".join(data[1:])
    return header, expression

def load_gene_name_map():
    mm9_gene_name_map = load_GENCODE_names("gencode.vM4.annotation.gff3")
    print("Finished loading GENCODE mouse gene names", file=sys.stderr)
    
    hg19_gene_name_map = load_hugo_names()
    print("Finished loading hugo names", file=sys.stderr)
    #tf_positions = load_encode_merged_tf_sites()

    gene_name_map = defaultdict(list)
    for gene_name, genes in mm9_gene_name_map.items():
        gene_name_map[gene_name].extend(genes)
    for gene_name, genes in hg19_gene_name_map.items():
        gene_name_map[gene_name].extend(genes)
    return dict(gene_name_map)

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

def main():    
    all_peaks = []
    for file_index, fname in enumerate(sys.argv[1:]):
        #if file_index > 4: break
        print( file_index, fname )
        
        pk_fp = io.TextIOWrapper(gzip.open(fname, 'rb'))
        peaks = load_bed( pk_fp, signal_index=4 )
        all_peaks.append(peaks)
        pk_fp.close()

    merged_peaks = merge_peaks(all_peaks, max, use_nonoverlapping_peaks=True )
    print( merged_peaks[0] )
    return

    tads = load_tads()
    print(tads)
    
    exp_header, expression = load_expression()
    print("Finished loading expression", file=sys.stderr)

    tf_positions = load_tf_sites("merged_TF_peaks.hg19_mm9.bed.gz")
    print("Finished loading all TF peaks", file=sys.stderr)
    
    encode_merged_tf_sites = load_encode_merged_tf_sites()
    print(len(encode_merged_tf_sites))
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
