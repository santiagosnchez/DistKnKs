#!/usr/bin/python env

import re
import urllib2
import sys
import math
import argparse

def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="DistKnKs.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""
    Quickly calculates Kn/Ks distances for a pair of aligned sequences based on
    Nei and Gojobori's 1986 algorithm.\n""",
    epilog="""
    Examples:
    python DistKnKs.py -f myAlignment.fas -c 1

    It assumes that the codons in the alignment are aligned. \"-c\" is the number code of codon table.
    See: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1

    The output is tabulated and printed to the screen.\n""")
    parser.add_argument(
    '--file', '-f', nargs="*", type=str,
    help='one or more alignments in FASTA format.')
    parser.add_argument(
    '--codon_table', '-c', type=str, default="1",
    help='the codon table number (default is 1, which stands for the standard genetic code.')
    args = parser.parse_args()
    gc = get_genetic_codes(args.codon_table)
    cc = codon_block_counts(gc)
    print "file,Kn_K80,Ks_K80,Kn_JC,Ks_JC"
    for f in args.file:
        file = f.split("/")[-1]
        file = re.sub("\.[fF][aA][sS]{0,1}[tT]{0,1}[aA]{0,1}$","",file)
        d = readfasta(f)
        if len(d.keys()) != 2:
            parser.error("Only two sequences per file are allowed.")
        var, total_syn, total = do_counts(d, cc)
        kn_k2p,ks_k2p,kn_jc,ks_jc = get_knks(var, total_syn, total)
        print "{0},{1:5f},{2:5f},{3:5f},{4:5f}".format(file,kn_k2p,ks_k2p,kn_jc,ks_jc)

# functions

def readfasta(file):
    data = {}
    with open(file, 'r') as f:
        lines = f.readlines()
        lines = map(lambda x: x.rstrip(), lines)
        ihead = map(lambda i: lines.index(i), filter(lambda k: ">" in k, lines))
        for i in range(len(ihead)):
            if ihead[i] != ihead[-1]:
                data[lines[ihead[i]][1:]] = ''.join(lines[ihead[i]+1:ihead[i+1]]).upper()
            else:
                data[lines[ihead[i]][1:]] = ''.join(lines[ihead[i]+1:]).upper()
    return data

def get_genetic_codes(numeric_code):
    # import urllib2
    # print "# Fetching NCBI for genetic code table: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1"
    link = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1'
    req = urllib2.Request(link)
    site = urllib2.urlopen(req).read().splitlines()
    code_numbers = [ x[8:x.index('.')] for x in filter(lambda x: '<h2>' in x, site) ]
    AAs = [ x[11:] for x in filter(lambda x: 'AAs' in x, site) ]
    Starts = [ x[11:] for x in filter(lambda x: 'Starts' in x, site) ]
    Base1 = [ x[11:] for x in filter(lambda x: 'Base1' in x, site) ]
    Base2 = [ x[11:] for x in filter(lambda x: 'Base2' in x, site) ]
    Base3 = [ x[11:] for x in filter(lambda x: 'Base3' in x, site) ]
    gen_codes = {}
    for i in range(len(code_numbers)):
        gen_codes[code_numbers[i]] = {}
        for j in range(len(AAs[i])):
            gen_codes[code_numbers[i]][Base1[i][j]+Base2[i][j]+Base3[i][j]] = AAs[i][j]
    return gen_codes[numeric_code]

def codon_block_counts(gc):
    cs = sorted(gc.keys())
    codon_blocks = [ cs[i:i+4] for i in range(0,len(cs),4) ]
    aa_by_block = []
    for b in codon_blocks:
        aa = [ gc[i] for i in b ]
        aa_by_block += [ i+str(aa.count(i)) for i in aa ]
    aa_by_block_dict = {}
    for i in range(len(cs)):
        aa_by_block_dict[cs[i]] = aa_by_block[i]
    return aa_by_block_dict

def get_codon_diff(cod):
    res = []
    n = 0
    for i in range(3):
        if cod[0][i] != cod[1][i]:
            n += 1
            if ''.join(sorted([cod[0][i],cod[1][i]])) == 'AG' or ''.join(sorted([cod[0][i],cod[1][i]])) == 'CT':
                res.append((i,"ts"))
            else:
                res.append((i,"tv"))
    res.append(n)
    return res

def do_counts(d, cc):
    seqlen = len(d[d.keys()[0]])
    var = {'S':{'ts':0,'tv':0},'N':{'ts':0,'tv':0}}
    total_syn = 0
    total = 0
    tmp = {'N':0,'S':0}
    for i in range(0,seqlen,3):
        cod = []
        for k in d.keys():
            cod.append(d[k][i:i+3])
        cod = list(set(cod))
        if len(cod) == 1:
            try:
                cc[cod[0]]
            except KeyError:
                pass
            else:
                if cc[cod[0]][0] == "*":
                    pass
                else:
                    total += 3
                    if cc[cod[0]][1] == '1':
                        pass
                    elif cc[cod[0]][1] == '2' or cc[cod[0]][1] == '3' or cc[cod[0]][1] == '4':
                        total_syn += 1
        else:
            try:
                cc[cod[0]]
                cc[cod[1]]
            except KeyError:
                pass
            else:
                if cc[cod[0]][0] == "*" or cc[cod[1]][0] == "*":
                    pass
                else:
                    total += 3
                    diff = get_codon_diff(cod)
                    if cc[cod[0]][0] == cc[cod[1]][0]: # check synonymous
                        if diff[-1] == 1:
                            var['S'][diff[0][1]] += 1
                            if (cc[cod[0]][1] == '4' and cc[cod[1]][1] == '4') or (cc[cod[0]][1] == '3' and cc[cod[1]][1] == '3'):
                                total_syn += 1
                            else:
                                total_syn += 2
                        else:
                            if cc[cod[0]][0] == 'L' or cc[cod[0]][0] == 'R':
                                var['S'][diff[0][1]] += 1
                                var['S'][diff[1][1]] += 1
                                total_syn += 2
                            elif gc[cod[0]] == 'S':
                                if diff[-1] == 2:
                                    var['N'][diff[0][1]] += 1
                                    var['N'][diff[1][1]] += 1
                                elif diff[-1] == 3:
                                    var['N'][diff[0][1]] += 1
                                    var['N'][diff[1][1]] += 1
                                    var['S'][diff[2][1]] += 1
                                    total_syn += 1
                    else: # nonsynonymous
                        if diff[-1] == 1:
                            var['N'][diff[0][1]] += 1
                        else:
                            if diff[-1] == 2:
                                if diff[1][0] == 2 and \
                                (cc[cod[0]][1] == '4' or cc[cod[1]][1] == '4'):
                                    var['S'][diff[1][1]] += 1
                                    var['N'][diff[0][1]] += 1
                                elif diff[1][0] == 2 and \
                                (cc[cod[0]][1] == '2' or cc[cod[1]][1] == '2') and \
                                diff[1][1] == 'ts':
                                    var['S'][diff[1][1]] += 1
                                    var['N'][diff[0][1]] += 1
                                else:
                                    var['N'][diff[0][1]] += 1
                                    var['N'][diff[1][1]] += 1
                            elif diff[-1] == 3:
                                if diff[2][0] == 2 and \
                                (cc[cod[0]][1] == '4' or cc[cod[1]][1] == '4'):
                                    var['S'][diff[2][1]] += 1
                                    var['N'][diff[0][1]] += 1
                                    var['N'][diff[1][1]] += 1
                                elif diff[1][0] == 2 and \
                                (cc[cod[0]][1] == '2' or cc[cod[1]][1] == '2') and \
                                diff[2][1] == 'ts':
                                    var['S'][diff[2][1]] += 1
                                    var['N'][diff[0][1]] += 1
                                    var['N'][diff[1][1]] += 1
                                else:
                                    print cod[0]
                                    var['N'][diff[0][1]] += 1
                                    var['N'][diff[1][1]] += 1
                                    var['N'][diff[2][1]] += 1
    return var, total_syn, total

def get_knks(var, total_syn, total):
    total_nonsyn = total - total_syn
    Pn = float(var['N']['ts'])/total_nonsyn
    Qn = float(var['N']['tv'])/total_nonsyn
    Ps = float(var['S']['ts'])/total_syn
    Qs = float(var['S']['tv'])/total_syn
    pn = float(sum(var['N'].values()))/total_nonsyn 
    ps = float(sum(var['S'].values()))/total_syn 
    kn_k2p = -0.5 * math.log((1-(2*Pn)-Qn)*math.sqrt(1-(2*Qn)))
    ks_k2p = -0.5 * math.log((1-(2*Ps)-Qs)*math.sqrt(1-(2*Qs)))
    kn_jc = -(3./4.) * math.log(1-((4./3.)*pn))
    ks_jc = -(3./4.) * math.log(1-((4./3.)*ps))
    return kn_k2p,ks_k2p,kn_jc,ks_jc

if __name__ == "__main__":
    main()













