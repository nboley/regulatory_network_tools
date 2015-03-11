import os, sys
import re

from collections import defaultdict

import pysam

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

def load_hugo_names(fname="hgnc_complete_set.txt"):
    header = None
    gene_name_map = defaultdict(list)
    with open(fname) as fp:
        for i, line in enumerate(fp):
            if i == 0: 
                header = line.split("\t")
                continue
            data = line.split("\t")
            gene_symbols = [data[1],] + [x.strip() for x in data[6].split(",")]
            ensemble_id = data[36]
            if ensemble_id == '': continue
            assert ensemble_id.startswith("EN")
            for symbol in gene_symbols:
                #assert symbol not in gene_name_map
                gene_name_map[symbol].append( ensemble_id )

    return dict(gene_name_map)

def load_encode_merged_tf_sites(fnames = ["chipseq_track_distal.bed", 
                                          "chipseq_track_proximal.bed"]):
    def extract_pos_and_tfs(line):
        data = line.strip().split()
        chrm, start, stop = data[0], int(data[1]), int(data[2])
        return (chrm, start, stop), set(x.split("(")[0] for x in data[-1].split(","))
    
    tf_positions = defaultdict(list)
    for fname in fnames:
        with open(fname) as fp:
            for i, line in enumerate(fp):
                region, tfs = extract_pos_and_tfs(line)
                for tf in tfs:
                    tf_positions[tf].append(region)
    return dict(tf_positions)

def load_enhancers(tf_fnames):
    # we build enhancers from overlapping intervals of TF binding
    def extract_pos_and_tfs(line):
        data = line.strip().split()
        chrm, start, stop = data[0], int(data[1]), int(data[2])
        return (chrm, start, stop), set(x.split("(")[0] for x in data[-1].split(","))
    
    tf_positions = defaultdict(list)
    for fname in fnames:
        with open(fname) as fp:
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
    print >> sys.stderr, "Finished loading GENCODE mouse gene names"
    
    hg19_gene_name_map = load_hugo_names()
    print >> sys.stderr, "Finished loading hugo names"
    #tf_positions = load_encode_merged_tf_sites()

    gene_name_map = defaultdict(list)
    for gene_name, genes in mm9_gene_name_map.iteritems():
        gene_name_map[gene_name].extend(genes)
    for gene_name, genes in hg19_gene_name_map.iteritems():
        gene_name_map[gene_name].extend(genes)
    return dict(gene_name_map)

def main():
    gene_name_map = load_gene_name_map()
    print >> sys.stderr, "Finished loading gene name map"
    
    exp_header, expression = load_expression()
    print >> sys.stderr, "Finished loading expression"

    tf_positions = load_tf_sites("merged_TF_peaks.hg19_mm9.bed.gz")
    print >> sys.stderr, "Finished loading all TF peaks"
    print len(tf_positions)
    
    encode_merged_tf_sites = load_encode_merged_tf_sites()
    print len(encode_merged_tf_sites)
    assert False
    
    print "type\tname\tgene_ens_id\t%s" % "\t".join(exp_header)
    for tf_symbol, positions in tf_positions.items():
        if tf_symbol not in gene_name_map:
            print >> sys.stderr, tf_symbol #, gene_name_map[tf]
        else:
            for gene_id in gene_name_map[tf_symbol]:
                if gene_id not in expression: 
                    print >> sys.stderr, ("EXP", tf_symbol, gene_id, 
                                          gene_name_map[tf_symbol] )
                else: 
                    print "%s\t%s\t%s\t%s" % (
                        "GENE_EXP", gene_id, tf_symbol, expression[gene_id])

if __name__ == '__main__':
    main()

#print bigbed_to_bed("chipseq_track_proximal.bb")
#CS = BigBedFile(open("chipseq_track_proximal.bb"))
#print dir(CS)
#print CS.get("chr1",0,1000000)
