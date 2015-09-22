#!/usr/bin/env python3
import os, sys, time

import requests, json

import urllib3
import psycopg2
conn = psycopg2.connect("host=mitra dbname=cisbp")

from itertools import chain

from collections import defaultdict, namedtuple

import re

import gzip, io

from multiprocessing import Value

PeakFile = namedtuple('PeakFiles', [
    'exp_id', 'target_id', 
    'sample_type', 'rep_key', 'bsid', 
    'assembly',
    'file_format', 'file_format_type', 'output_type', 'file_loc'])

TargetInfo = namedtuple('TargetInfo', [
    'target_id', 
    'organism',
    'tf_name', 
    'uniprot_ids',
    'gene_name',
    'ensemble_ids',
    'cisbp_id'])

BASE_URL = "https://www.encodeproject.org/"

def group_overlapping_intervals(intervals):
    intervals = sorted(intervals)
    if len(intervals) == 0: return []
    
    curr_start, curr_stop = intervals[0][0], intervals[0][1]
    merged_intervals = [([curr_start, curr_stop], [intervals[0],]),]
    for interval in intervals[1:]:
        if interval[0] > curr_stop:
            curr_start, curr_stop = interval[0], interval[1]
            merged_intervals.append(
                ([curr_start, curr_stop], [interval,]) )
        else:
            curr_stop = max(curr_stop, interval[1])
            merged_intervals[-1][1].append(interval)

    return merged_intervals

def get_TF_name_and_label_from_fname(fname):
    TF_label = os.path.basename(fname).split("_")[0]
    TF_name = os.path.basename(fname).split("_")[1]
    return TF_name, TF_label

def load_peaks_from_fname(fname, peaks=None):
    if peaks == None: peaks = defaultdict(list)
    
    TF_name, TF_label = get_TF_name_and_label_from_fname(fname)
    
    if TF_name.endswith('human'):
        assembly_prefix = 'hg19_' 
    elif TF_name.endswith('mouse'):
        assembly_prefix = 'mm9_' 
    
    if fname.endswith("gz"): 
        fp = io.TextIOWrapper(gzip.open(fname, 'rb'))
    else: 
        fp = open(fname)

    for i, line in enumerate(fp):
        data = line.split()
        chrm = assembly_prefix + data[0]
        peaks[chrm].append((int(data[1]), int(data[2]), TF_label))

    return peaks

def load_peaks_merged_by_TF(fnames):
    tf_grpd_peaks = defaultdict(lambda: defaultdict(list))
    for i, fname in enumerate(fnames):
        print( "Loading %i/%i " % (i, len(fnames)) )
        TF_name, TF_label = get_TF_name_and_label_from_fname(fname)
        load_peaks_from_fname(fname, tf_grpd_peaks[TF_label])
    
    tf_merged_peaks = {}
    for i, (tf, tf_peaks) in enumerate(tf_grpd_peaks.items()):
        print( "Grouping %i/%i " % (i, len(tf_grpd_peaks)) )
        tf_merged_peaks[tf] = {}
        for contig, intervals in tf_peaks.items():
            intervals.sort()
            tf_merged_peaks[tf][contig] = [
                (x[0], tf) for x in group_overlapping_intervals(intervals) ]

    return tf_grpd_peaks

def flatten_peaks(peaks):
    merged_peaks = defaultdict(list)
    for contig, contig_peaks in peaks.items():
        merged_contig_peaks = group_overlapping_intervals(contig_peaks)
        for (start, stop), peaks in merged_contig_peaks:
            tfs = tuple(sorted(set(pk[2] for pk in peaks)))
            merged_peaks[contig].append( (start, stop, tfs) )
    
    return merged_peaks

def load_and_merge_peaks(fnames):
    TF_merged_peaks = load_peaks_merged_by_TF(fnames)

    all_peaks = defaultdict(list)

    all_peaks = []
    ofnames = []
    for tf, tf_peaks in TF_merged_peaks.items():
        print( "Merging ", tf )
        #for contig, contig_peaks in tf_peaks.items():
        #    all_peaks[contig].extend( contig_peaks )
        #continue
        ofname = "{}.mergedpeaks.bed".format(tf)
        with open(ofname, "w") as ofp:
            tf_peaks = flatten_peaks(tf_peaks)
            for contig, contig_peaks in sorted(tf_peaks.items()):
                for (start, stop, tfs) in contig_peaks:
                    assert len(tfs) == 1 and tfs[0] == tf
                    all_peaks.append((contig, start, stop, tf))
                    print("\t".join(str(x) for x in (
                        contig, start, stop, ",".join(tfs))),
                          file=ofp) 
        os.system("bgzip {fname}; tabix -p bed {fname}.gz".format(fname=ofname))
        ofnames.append(ofname)
    
    all_peaks.sort()
    ofname = "all_peaks.bed"
    with open(ofname, "w") as ofp:
        for peak in all_peaks:
            print( "\t".join(map(str, peak)), file=ofp )
    os.system("bgzip {fname}; tabix -p bed {fname}.gz".format(fname=ofname))
    
    ## merge all peaks
    #with open("mergedpeaks.bed".format(tf), "w") as ofp:
    #    all_peaks = flatten_peaks(all_peaks)
    #    for contig, contig_peaks in sorted(all_peaks.items()):
    #        for (start, stop, tfs) in contig_peaks:
    #            print("\t".join(str(x) for x in (
    #                contig, start, stop, ",".join(tfs))),
    #                  file=ofp) 

    return ofnames


def find_called_peaks(experiment_id, only_merged=True):
    URL = "https://www.encodeproject.org/experiments/{}/".format(experiment_id)
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()

    target_id = response_json_dict['target']['@id']
    target_label = response_json_dict['target']['label']

    # build the replicate mapping
    replicates = {}
    for rep_rec in response_json_dict['replicates']:
        # skip replicates without an associated library...
        if 'library' not in rep_rec:
            return
        # skip treated libraries
        if len(rep_rec['library']['treatments']) > 0:
            continue
        key = (rep_rec['biological_replicate_number'], 
               rep_rec['technical_replicate_number'])
        replicates[key] = rep_rec

    sample_type = response_json_dict['biosample_term_name']
    
    for file_rec in response_json_dict['files']:
        file_format = file_rec['file_format']
        if file_format not in ['bed', ]:
            continue

        file_format_type = file_rec['file_format_type']
        output_type = file_rec['output_type']
        
        if 'replicate' not in file_rec: 
            rep_key = 'merged'
            bsid = 'merged'
        else:
            if only_merged: continue
            rep_key = (file_rec['replicate']['biological_replicate_number'],
                       file_rec['replicate']['technical_replicate_number'] )
            bsid = replicates[rep_key]['library']['biosample']['accession']
        file_loc = "http://encodeproject.org" + file_rec['href']
        assembly = file_rec['assembly']
        yield PeakFile(experiment_id, target_id, 
                       sample_type, rep_key, bsid,
                       assembly,
                       file_format, file_format_type, output_type, file_loc)

    return 

def find_chipseq_experiments(assemblies=['mm9', 'hg19']): # 'hg19', 
    URL = "https://www.encodeproject.org/search/?type=experiment&assay_term_name=ChIP-seq&{}&target.investigated_as=transcription%20factor&limit=all&format=json".format(
        "&".join("assembly=%s"%x for x in assemblies) )
        
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()
    biosamples = set()
    for experiment in response_json_dict['@graph']:
        yield experiment['@id'].split("/")[-2]
    return 

def find_cisbp_tfids(species, tf_name, uniprot_ids, ensemble_ids):
    cur = conn.cursor()
    query = "SELECT tf_id FROM tfs WHERE dbid IN %s;"
    #print(cur.mogrify(
    #    query, [tuple(chain(uniprot_ids, chain(*ensemble_ids))),]))
    rv = []
    prot_and_gene_ids = tuple(chain(uniprot_ids, chain(*ensemble_ids)))
    if len(prot_and_gene_ids) > 0:
        cur.execute(query, [prot_and_gene_ids,]) 
        rv = [x[0] for x in cur.fetchall()]
    # if we can't find a reference from the uniprot or ensemble ids,
    # then try with the tf name
    if len(rv) != 1:
        query = "SELECT tf_id FROM tfs WHERE tf_species = %s and tf_name = %s;"
        res = cur.execute(
            query, (species, tf_name))
        rv = [x[0] for x in cur.fetchall()]
    
    # if we still an't find a match, try with the upper case
    if len(rv) != 1:
        query = "SELECT tf_id FROM tfs WHERE tf_species = %s and upper(tf_name) = upper(%s);"
        res = cur.execute(
            query, (species.replace(" ", "_"), tf_name))
        rv = [x[0] for x in cur.fetchall()]
        
    return rv

def find_target_info(target_id):
    URL = "https://www.encodeproject.org/{}?format=json".format(target_id)
    response = requests.get(URL, headers={'accept': 'application/json'})
    response_json_dict = response.json()
    organism = response_json_dict['organism']['scientific_name'].replace(" ", "_")
    tf_name = response_json_dict['label']
    uniprot_ids = [x[10:] for x in response_json_dict['dbxref']
                   if x.startswith("UniProtKB:")]
    gene_name = response_json_dict['gene_name']
    ensemble_ids = sorted(
        get_ensemble_genes_associated_with_uniprot_id(uniprot_id) 
        for uniprot_id in uniprot_ids)
    cisbp_ids = find_cisbp_tfids(organism, tf_name, uniprot_ids, ensemble_ids)
    if len(cisbp_ids) == 0:
        cisbp_id = None
    else:
        assert len(cisbp_ids) == 1
        cisbp_id = cisbp_ids[0]
    rv = TargetInfo(target_id, organism, 
                    tf_name, uniprot_ids, 
                    gene_name, ensemble_ids, 
                    cisbp_id)
    return rv

http = urllib3.PoolManager()
def get_ensemble_genes_associated_with_uniprot_id(uniprot_id):
    ens_id_pat = '<property type="gene ID" value="(ENS.*?)"/>'
    res = http.request(
        "GET", "http://www.uniprot.org/uniprot/%s.xml" % uniprot_id)
    #print( res.read().decode('utf-8') )
    gene_ids = set(re.findall(ens_id_pat, res.read().decode('utf-8')))
    return sorted(gene_ids)

def find_peaks_and_group_by_target(
        only_merged=True, prefer_uniformly_processed=True):
    # find all chipseq experiments and group them by target
    targets = defaultdict(list)
    chipseq_exps = list(find_chipseq_experiments())
    for i, exp_id in enumerate(chipseq_exps):
        if i > 1: break
        print( i, len(chipseq_exps), exp_id )
        for rep_i, res in enumerate(find_called_peaks(exp_id, only_merged)):
            targets[res.target_id].append(res)
            print( res.file_format, res.file_format_type, "-", res.output_type )
        print()
    
    # remove redudant experiments
    for i, (target, file_data) in enumerate(targets.items()):
        any_uniformly_processed = any(
            f.output_type == 'UniformlyProcessedPeakCalls'
            for f in file_data )

        new_file_data = [
            f for f in file_data 
            if (not prefer_uniformly_processed 
                or not any_uniformly_processed 
                or f.output_type=='UniformlyProcessedPeakCalls')
            and (
                    not only_merged
                    or f.bsid == 'merged'
            )
        ]
        print( i, len(targets), target )
        yield target, new_file_data
    return

def download_sort_and_index_tfs():
    res = []
    for target, files in find_peaks_and_group_by_target():
        tf_label, tf_name, uniprot_ids = find_target_info(target)
        all_genes = [ get_ensemble_genes_associated_with_uniprot_id(uid)
                      for uid in uniprot_ids ]
        for file_data in files:
            ofname = file_data.file_loc.split("/")[-1]

            human_readable_ofname = (
                "_".join((tf_label, tf_name, 
                          file_data.output_type, 
                          file_data.sample_type.replace(" ", "-"))) \
                + ".EXP-%s.%s.%s.%s.%s.bgz" % (
                    file_data.exp_id, 
                    file_data.rep_key, 
                    file_data.output_type,
                    file_data.file_format,
                    file_data.file_type))
            human_readable_ofname = human_readable_ofname.replace(
                "/", "#FWDSLASH#")

            #print "Downloading %i/%i" % (i+1, len(chipseq_exps))
            download_cmd = "wget --quiet {URL}"
            sort_and_compress_cmd = \
               "zcat {FNAME} | sort -k1,1 -k2n -k3n | bgzip -c > {HR_FNAME}"
            mv_cmd = "rm {FNAME}"
            index_cmd = "tabix -p bed {HR_FNAME}"
            #print cmd
            if "/" in human_readable_ofname:
                cmd = " ; ".join(
                    (download_cmd, sort_and_compress_cmd, mv_cmd, index_cmd)
                ).format(URL=BASE_URL+file_data.file_loc, 
                         FNAME=ofname, 
                         HR_FNAME=human_readable_ofname )
                print( cmd )
                os.system(cmd)
                res.append(
                    (file_data.exp_id, tf_name, tf_label, 
                     ",".join(uniprot_ids), 
                     ",".join(",".join(genes) for genes in all_genes), 
                     file_data.sample_type, 
                     file_data.output_type,
                     BASE_URL[:-1]+file_data.file_loc, 
                     human_readable_ofname))

def find_or_insert_experiment_from_called_peaks(called_peaks):
    cur = conn.cursor()

    # first check to see if the experiment is already in the DB. If
    # so, we are done
    encode_exp_ids = set(
        called_peak.exp_id for called_peak in called_peaks)
    assert len(encode_exp_ids) == 1, str(encode_exp_ids)
    encode_exp_id = encode_exp_ids.pop()
    query = "SELECT encode_experiment_id FROM encode_chipseq_experiments WHERE encode_experiment_id = %s;"
    cur.execute(query, [encode_exp_id,])
    # if the experiment is already there, we are done
    if len(cur.fetchall()) == 1:
        return
    # otherwise, insert everything
    encode_target_ids = set(
        called_peak.target_id for called_peak in called_peaks)
    assert len(encode_target_ids) == 1
    encode_target_id = encode_target_ids.pop()
    # find our associated target id 
    query = "SELECT chipseq_target_id FROM chipseq_targets WHERE encode_target_id = %s"
    cur.execute(query, [encode_target_id,])
    res = cur.fetchall()
    # if we can't find a matching tf id, insert it
    if len(res) == 0:
        target_info = find_target_info(encode_target_id)
        query = "INSERT INTO chipseq_targets (encode_target_id, tf_id, organism, tf_name, uniprot_ids, ensemble_gene_ids) VALUES (%s, %s, %s, %s, %s, %s) RETURNING chipseq_target_id"
        cur.execute(query, [encode_target_id, 
                            target_info.cisbp_id, 
                            target_info.organism,
                            target_info.tf_name,
                            target_info.uniprot_ids, 
                            target_info.ensemble_ids])
        res = cur.fetchall()
    assert len(res) == 1
    target_id = res[0][0]

    sample_types = set(
        called_peak.sample_type for called_peak in called_peaks)
    assert len(sample_types) == 1
    sample_type = sample_types.pop()
    # add the experiment data
    query = "INSERT INTO encode_chipseq_experiments " \
          + "(encode_experiment_id, target, sample_type) " \
          + "VALUES (%s, %s, %s)"
    cur.execute(query, [
        encode_exp_id, target_id, sample_type])
    return

def encode_exp_is_in_db(exp_id):
    cur = conn.cursor()
    query = "SELECT encode_experiment_id FROM encode_chipseq_experiments WHERE encode_experiment_id = %s;"
    cur.execute(query, [exp_id,])
    # if the experiment is already there, we are done
    if len(cur.fetchall()) == 1:
        return True
    return False

def insert_chipseq_experiment_into_db(exp_id):
    """
    CREATE TABLE chipseq_targets (
        id SERIAL PRIMARY KEY,
        encode_id text UNIQUE,
        tf_id text NOT NULL,
        organism text NOT NULL,
        tf_name text NOT NULL,
        uniprot_ids text[] NOT NULL,
        ensemble_gene_ids text[][] NOT NULL
    );

    CREATE TABLE encode_chipseq_experiments (
        id text PRIMARY KEY,
        target int NOT NULL REFERENCES chipseq_targets(id),
        sample_type text NOT NULL
    );

    CREATE TABLE encode_chipseq_peak_files (
        experiment_id text NOT NULL REFERENCES encode_chipseq_experiments(id),
        bsid text NOT NULL,
        rep_key text NOT NULL,
        file_format text NOT NULL,
        file_format_type text NOT NULL,
        file_output_type text NOT NULL,
        remote_filename text NOT NULL,
        local_filename text NOT NULL
    );
    """
    if encode_exp_is_in_db(exp_id):
        return
    
    called_peaks = list(find_called_peaks(exp_id, only_merged=False))
    if len(called_peaks) == 0: return
    
    # insert the experiment and target into the DB if necessary
    num_inserted = find_or_insert_experiment_from_called_peaks(called_peaks)
    cur = conn.cursor()
    for called_peak in called_peaks:
        # add the peak data
        query = "INSERT INTO encode_chipseq_peak_files " \
              + "(encode_experiment_id, bsid, rep_key, file_format, file_format_type, file_output_type, remote_filename) " \
              + "VALUES (%s, %s, %s, %s, %s, %s, %s)"
        try: 
            cur.execute(query, [
                called_peak.exp_id, called_peak.bsid, called_peak.rep_key,
                called_peak.file_format, called_peak.file_format_type, called_peak.output_type,
                called_peak.file_loc])
        except psycopg2.IntegrityError:
            print( "ERROR" )
            raise
            pass
    conn.commit()
    return

if __name__ == '__main__':
    all_chipseq_exps = list(find_chipseq_experiments())
    for i, exp_id in enumerate(all_chipseq_exps):
        print( i, len(all_chipseq_exps), exp_id )
        try:
            insert_chipseq_experiment_into_db( exp_id )
        except:
            print( "ERROR with %s" % exp_id )
            raise
    #for target, files in find_peaks_and_group_by_target(
    #        only_merged=False, prefer_uniformly_processed=False):
    #    insert_chipseq_data_into_db(target, files)
    #res = download_sort_and_index_tfs()
    #with open("tfs.txt", "w") as ofp:
    #    for entry in res:
    #        print( "\t".join(entry), file=ofp )
    #call_peaks_for_experiment( sys.argv[1] )
    #load_and_merge_peaks( sys.argv[1:] )
