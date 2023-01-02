'''Utilities.'''


# --- Workflow setup --- #


import glob
import gzip
import os
from os import makedirs, symlink
from os.path import abspath, basename, exists, join
import pandas as pd
import shutil


def ingest_samples(samples, tmp):
    df = pd.read_csv(samples, header = 0, index_col = 0) # name, mag_dir, fwd, rev
    s = list(df.index)
    lst = df.values.tolist()
    for i,l in enumerate(lst):
        if not exists(join(tmp, s[i])): # Make a temporary directory for all of the MAGs in the sample
            makedirs(join(tmp, s[i]))
            with open(join(tmp, s[i] + '.out'), 'w') as f_out: # Enables the CheckM rule to run
                for m in glob.glob(l[0] + '/*.fa*'):
                    if 'unbinned' not in m:
                        prefix = basename(m).split('.')[1]
                        symlink(abspath(m), join(tmp, s[i], prefix + '.fa'))
                        f_out.write(str(prefix) + '\n')
            symlink(abspath(l[1]), join(tmp, s[i] + '_1.fastq.gz'))
            symlink(abspath(l[2]), join(tmp, s[i] + '_2.fastq.gz'))
    return s


def check_make(d):
    if not exists(d):
        makedirs(d)


class Workflow_Dirs:
    '''Management of the working directory tree.'''
    OUT = ''
    TMP = ''
    LOG = ''

    def __init__(self, work_dir, module):
        self.OUT = join(work_dir, module)
        self.TMP = join(work_dir, 'tmp') 
        self.LOG = join(work_dir, 'logs') 
        check_make(self.OUT)
        out_dirs = ['0_checkm', '1_cmseq', '2_gtdbtk', '3_dnadiff', '4_quast', 'final_reports']
        for d in out_dirs: 
            check_make(join(self.OUT, d))
        # Add a subdirectory for symlinked-in input files
        check_make(self.TMP)
        # Add custom subdirectories to organize rule logs
        check_make(self.LOG)
        log_dirs = ['checkm', 'gtdbtk', 'cmseq', 'dnadiff', 'quast']
        for d in log_dirs: 
            check_make(join(self.LOG, d))


def cleanup_files(work_dir, df):
    smps = list(df.index)
    for s in smps: 
        os.system('rm  ' + join(dirs.OUT, '1_cmseq', s, '*bam*'))


def print_cmds(log):
    fo = basename(log).split('.')[0] + '.cmds'
    lines = open(log, 'r').read().split('\n')
    fi = [l for l in lines if l != '']
    write = False
    with open(fo, 'w') as f_out:
        for l in fi:
            if 'rule' in l:
                f_out.write('# ' + l.strip().replace('rule ', '').replace(':', '') + '\n')
            if 'wildcards' in l: 
                f_out.write('# ' + l.strip().replace('wildcards: ', '') + '\n')
            if 'resources' in l:
                write = True 
                l = ''
            if '[' in l: 
                write = False 
            if write:
                f_out.write(l.strip() + '\n')
            if 'rule make_config' in l:
                break


# --- Workflow functions --- #


from cmseq import CMSEQ_DEFAULTS, BamFile
import numpy as np
from os import getenv
from os.path import getsize


def add_bin_num(fi, bin_num, fo):
    ctg_num = 0
    with open(fi,'r') as f_in, open(fo, 'w') as f_out:
        for line in f_in:
            if line[0] == '>':
                ctg_name = line.strip('\n').replace('>','')
                new_ctg_name = ">%s_%i\t%s" % (bin_num, ctg_num, ctg_name) 
                f_out.write(new_ctg_name + '\n')
                # print(new_ctg_name)
                ctg_num += 1
            else:
                f_out.write(line)


def get_bin_nums(s, d):
    tmp = open(join(d, str(s) + '.out'), 'r').readlines()
    # print(*[i.strip() for i in tmp])
    return [i.strip() for i in tmp]


def pair_mag_refs(row, out_dir):
    r = row['fastani_reference']
    r_path = 'None'
    if str(r) != 'nan':
        GTDB = getenv('GTDBTK_DATA_PATH')
        parts = r.split('_')
        r_path = join(GTDB, 'fastani/database', parts[0], parts[1][0:3], parts[1][3:6], parts[1][6:9], r + '_genomic.fna.gz')
    with open(join(out_dir, str(row['user_genome']) + '.ref'), 'w') as f_out:
        f_out.write(r_path + '\n')


def parse_dnadiff(fi, fo):
    if getsize(fi) != 0: # If the MAG was classified as a species
        first_line = open(fi, 'r').readlines()[0].split()
        ref = first_line[0]
        quer = first_line[1]
        with open(fi, 'r') as f_in:
            for line in f_in:
                if "TotalBases" in line:
                    cols = line.strip().split()
                    lenref = int(cols[1])
                    lenquer = int(cols[2])
                if "AlignedBases" in line:
                    cols = line.strip().split()
                    aliref = cols[1].split("(")[-1].split("%")[0]
                    alique = cols[2].split("(")[-1].split("%")[0]
                if "AvgIdentity" in line:
                    cols = line.strip().split()
                    ident = float(cols[1])
            output = "%s\t%s\t%i\t%.2f\t%i\t%.2f\t%.2f" % (quer, ref, lenref, float(aliref), lenquer, float(alique), float(ident))
    else:
        quer = fi.split('/')[-1].split('.')[0]
        output = quer + '\tNone\t0\t0.00\t0\t0.00\t0.00'
    with open(fo, 'w') as f_out:
        f_out.write(output + '\n')


def polymut_from_cmseq(bam, mag, gff, fo, minqual, mincov, dominant_frq_thrsh):
    bf = BamFile(bam, filterInputList = mag, sort = False, index = False, stepper = 'nofilter', \
        minlen = CMSEQ_DEFAULTS.minlen)
    bf.parse_gff(gff)
    # Calculate the strain heterogeneity for each contig in the MAG
    outputDicts = []
    for i in bf.get_contigs_obj():
        dominanceArray, mutationStats = i.easy_polymorphism_rate(minqual = minqual, mincov = mincov,dominant_frq_thrsh = dominant_frq_thrsh)
        outputDicts.append({'Ref':i.name, 'DN':mutationStats['DN'],'DS':mutationStats['DS'],'D?':mutationStats['D?'], "consid_pos":len([x for x in dominanceArray if not np.isnan(x)])})
    out_df = pd.DataFrame.from_dict(outputDicts).set_index('Ref')
    # Print the average strain heterogeneity across the entire MAG
    line = "%.4f\t%.4f\t%.4f\t" % (float(np.sum(out_df["DN"])), float(np.sum(out_df["DS"])), float(sum(out_df["consid_pos"])))
    prop_DN = float(np.sum(out_df["DN"]))*100/float(sum(out_df["consid_pos"])) if float(sum(out_df["consid_pos"])) != 0 else 'NA'
    line += str(prop_DN)
    with open(fo, 'w') as f_out:
        f_out.write(line + '\n')


def aggregate_quast(fi_lst, fo):
    df_lst = []
    unc_mags = []
    for fi in fi_lst:
        if getsize(fi) != 0:
            df_lst.append(pd.read_csv(fi, header = 0, sep = '\t').transpose())
        else: # If QUAST report is empty, then no classification
            unc_mags.append(fi.split('/')[-2])
    if len(df_lst):
        raw_df = pd.concat(df_lst)
        raw_df.reset_index(inplace = True)
        df = raw_df[['index', '# contigs', 'Total length', 'Genome fraction (%)', 'NG50', 'NA50', '# misassemblies', '# misassembled contigs', 'Misassembled contigs length', '# unaligned contigs', 'Unaligned length']]
        col_names = ['mag', 'num_ctgs', 'size', 'genome_fraction', 'NG50', 'NA50', 'num_misassemb', 'num_misassemb_ctgs', 'misassemb_ctg_len', 'num_unaln_ctgs', 'unaln_len']
        df.columns = col_names
        df['prop_misassemb_ctgs'] = df.apply(lambda row : float(row[7]/row[1]), axis = 1) 
        df['prop_misassemb_len'] = df.apply(lambda row : float(row[8]/row[2]), axis = 1) 
        df['prop_unaln_ctgs'] = df.apply(lambda row : float(row[9]/row[1]), axis = 1) 
        df['prop_misassemb_ctgs'] = df.apply(lambda row : float(row[10]/row[2]), axis = 1) 
        for m in unc_mags: 
            unc_mag_row  = {} # Create empty rows for unclassified MAGs
            for c in col_names + ['prop_misassemb_ctgs', 'prop_misassemb_len', 'prop_unaln_ctgs', 'prop_unaln_len']:
                unc_mag_row[c] = 0
            unc_mag_row['mag'] = m
            df = df.append(unc_mag_row, ignore_index = True)
        fin_df = df[['mag', 'genome_fraction', 'NG50', 'NA50', 'num_misassemb', 'prop_misassemb_ctgs', 'prop_misassemb_len', 'prop_unaln_ctgs', 'prop_unaln_len']]
        df.to_csv(fo, header = True, index = False)
    else:
        open(str(fo), 'w').close()


def summarize_reports(checkm, n50_sz, cmseq, gtdb, diff, quast, output):
    # Load completeness and contamination (mag,completeness,contamination)
    summ_df = pd.read_csv(checkm, header = 0)
    summ_df['mag'] = summ_df.apply(lambda row : str(row[0]).split('.')[0], axis = 1) # Reshape X.fa into X
    # Load strain heterogeneity
    cmsq_df = pd.read_csv(cmseq, sep = '\t', header = None)
    cmsq_df[0] = cmsq_df.apply(lambda row : str(row[0]).split('/')[-1].split('.')[0], axis = 1) # Reshape /path/to/X.* into X
    cmsq_df.columns = ['mag', 'strain_heterogeneity']
    # Load taxonomic classification (all levels)
    raw_df = pd.read_csv(gtdb, sep = '\t')
    gtdb_df = raw_df[['user_genome', 'classification']]
    gtdb_df.columns = ['mag', 'classification']
    gtdb_df['mag'] = gtdb_df['mag'].astype(str) # Otherwise, interpreted as int
    # Load MAG-reference aligned length, reference genome coverage, ANI
    raw_df = pd.read_csv(diff, sep = '\t', header = None)
    raw_df[0] = raw_df.apply(lambda row : str(row[0]).split('/')[-1].split('.')[0], axis = 1) # Reshape /path/to/X.* into X
    raw_df.columns = ['mag', 'ref', 'ref_len', 'ref_cov', 'bin_size', 'bin_cov', 'ANI']
    diff_df = raw_df[['mag', 'bin_cov', 'ref_cov', 'ANI']]
    # Load N50, MAG size (in terms of bp), GC
    size_dct = {}
    for line in open(n50_sz):
        info = line.strip().split('\t')
        size_dct[info[0]] = eval(info[1])
    raw_df = pd.DataFrame(data = size_dct).transpose()
    size_df = raw_df[['GC', 'N50 (contigs)', '# contigs', 'Genome size']]
    size_df = size_df.reset_index()
    size_df.columns = ['mag', 'GC', 'N50', 'num_ctgs', 'size']
    # Load reference genome-based completion, misassembly, and unaligned statistics from QUAST report
    quas_df = pd.read_csv(quast, header = 0)
    # Put all of the dataframes together and output
    for df in [cmsq_df, gtdb_df, diff_df, size_df, quas_df]:
        summ_df = pd.merge(summ_df, df, on = 'mag', how = 'outer')
    # Add in standard quality thresholds
    summ_df['Quality'] = np.where((summ_df['completeness'] >= 90) & (summ_df['contamination'] <= 5), 'High quality', 'Medium quality') # Quality labels for DNAdiff figure
    summ_df.to_csv(output, header = True, index = False)

