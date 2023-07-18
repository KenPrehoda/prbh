#!/usr/bin/env python
import gzip
import argparse
from dataclasses import dataclass

HYD_RES = {'A':0.05,'F':0.2,'I':0.2,'L':0.2,'M':0.25,'V':0.1,'W':0.1,'Y':0.1}
BASIC_RES = {'R':0.1,'K':0.07,'H':0.04,'D':-0.06,'E':-0.06,'Y':0.0}
M3_RES = {'F':0.25,'L':0.25,'M':0.25,'I':0.15,'V':0.1,'W':0.1,
        'Y':0.1,'A':0.05,'R':0.1,'K':0.1,'H':0.05}
M2_RES = {'R':0.35,'K':0.15,'H':0.05}
BH_DICT = {'O':-.17, 'U':0,'X':0,'A': -.17, 'C':0.24, 'D':-1.23, 'E': -2.02, 'F':1.13,
            'G':-.01, 'H':-.17, 'I':.31, 'K':2,'L':.56,'M':.23,'N':-.42,'P':- .45,
            'Q':-.58,'R':2, 'S':-.13,'T':-.14,'V':-.07,'W':1.85,'Y':.94}

@dataclass
class fasta:
    desc:str
    seq:str

def fasta_file(it): 
    """generator that returns an iterator over the records in a FASTA format file """
    d = fasta('','')
    for i in it:
        if i[0] == '>':
            if len(d.seq) > 0:
                yield d
                d = fasta('','')
            d.desc = i[1:].strip()
        else:
            d.seq += i.strip()
    if len(d.seq) > 0: 
        yield d

class prbh_motif:
    def __init__(self,bh,apkcs,score):
        self.bh, self.apkcs, self.score = (bh,apkcs,score)
    def __eq__(self,other):
        return self.bh == other.bh and (self.score -
                other.score) < 0.00001 and self.apkcs[0] == other.apkcs[0]
class bh_motif:
    def __init__(self,start,end,height,area):
        self.start, self.end, self.height, self.area = (start,end,height,area)
    def __eq__(self,other):
        return self.start == other.start and (self.height - 
            other.height) < 0.00001 and (self.area - other.area) < 0.00001
class aPKC_site:
    def __init__(self,pos,score,seq):
        self.pos, self.score, self.seq = (pos,score,seq)
    def __eq__(self,other):
        return self.pos == other.pos and (self.score - other.score) < 0.00001 and self.seq == other.seq
def max_seq(seq,score_dict): 
    hyd = 0.0
    for i in seq:
        hyd = max(hyd,score_dict.get(i,0.0))
    return hyd
def calc_aPKC_sites(seq, threshold=0.6): 
    sites = []
    for pos,char in enumerate(str(seq).upper()): 
        if char == 'S' or char == 'T':
            score = 0.0
            # pos + 1 hydrophobic if pos+1 < len(seq):
            if pos+1 < len(seq):
                score += HYD_RES.get(seq[pos+1],-0.4)
            basicTemp = 0.0
            for i in range(2,7):
                if pos+i+1 < len(seq):
                    basicTemp += BASIC_RES.get(seq[pos+i],0.0)
            score += basicTemp * 0.8
            if pos-2 > 0:
                score += M2_RES.get(seq[pos-2],0.)*1.4
            if pos-3 > 0:
                score += M3_RES.get(seq[pos-3],0.)
            hyd = max_seq(seq[max(0,pos-9):max(0,pos-4)],HYD_RES)
            bas = max_seq(seq[max(0,pos-9):max(0,pos-4)],BASIC_RES) 
            score += (hyd + bas) * 0.5
            if score > threshold:
                sites.append(aPKC_site(pos,score,seq[max(0,pos-
                                        7):min(pos+17,len(seq))]))
    return sites
def running_bh_sum(seq, window_size): 
    pos = 0
    if len(seq) < window_size: 
        return
    running_sum = sum([BH_DICT.get(x,0) for x in seq[:window_size]]) 
    yield running_sum
    for c in seq[window_size:]:
        running_sum += BH_DICT.get(c,0) - BH_DICT.get(seq[pos],0) 
        pos += 1
        yield running_sum

def calc_bh(seq, window_size=19, threshold=0.6): 
    motifs = []
    height = 0.0
    area = 0.0
    startRes = -1
    for pos,bsum in enumerate(running_bh_sum(seq,window_size)): 
        bh_score = bsum/window_size
        if bh_score > threshold:
            if startRes == -1: #new BH 
                startRes = pos
            height = max(bh_score,height)
            area += bsum
        elif startRes != -1: #close out BH
            motifs.append(bh_motif(startRes,pos+window_size,height,area)) 
            startRes = -1
            height = 0.0
            area = 0.0
    return motifs
def calc_bh_score(seq,window_size=19): 
    scores = []
    for pos,bsum in enumerate(running_bh_sum(seq,window_size)): 
        bh_score = bsum/window_size
        scores.append(bh_score)
    return scores
def calc_prbh(seq,threshold=1.6): 
    prbhs = []
    bhs = []
    sites = calc_aPKC_sites(seq) 
    if sites:
        bhs = calc_bh(seq) 
    for bh in bhs:
        phos = []
        for apkc in sites:
            if apkc.pos > bh.start and apkc.pos < bh.end: 
                phos.append(apkc)
        if len(phos):
            score = bh.height + max([x.score for x in phos]) + len(phos)*0.2 
            if bh.area < len(phos) * 400.:
                if bh.area/90. > 1: 
                    if score > threshold:
                        prbhs.append(prbh_motif(bh,phos,score))
    return prbhs

class tagged_prbh(object):
    def __init__(self,prbh,fasta):
        self.prbh, self.fasta = (prbh,fasta) 
    def __eq__(self,other):
        return self.prbh.score-other.prbh.score < 0.00001 and self.prbh.apkcs[0].seq == other.prbh.apkcs[0].seq
    def __hash__(self):
        return hash((self.prbh.score, self.prbh.apkcs[0].seq))
def fasta_prbh(fasta): 
    if fasta[-2:] == 'gz':
        fi = gzip.GzipFile(fasta,'r')
    else:
        fi = open(fasta,'r')
    hits = [] 
    with fi as f:
        for rec in fasta_file(f): 
            prbhs = calc_prbh(rec.seq) 
            for i in prbhs:
                hits.append(tagged_prbh(i,rec)) 
    return hits

def format_prbh(prbh):
    return f'{prbh.fasta.desc[:30]}\t{prbh.prbh.bh.start}\t{prbh.prbh.bh.end}\t{prbh.prbh.score:.3f}'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='prbh',
                    description='implements the PhosphoRegulated Basic Hydrophobic (PRBH) motif finder as described in Bailey and Prehoda, Dev Cell (2018)')
    parser.add_argument('fasta_file',help="File in FASTA format (can be gzipped) that contains protein sequences")
    parser.add_argument('-o', '--output', help="Optional output file (otherwise output it written to stdout)")
    args = parser.parse_args()
    prbhs = sorted(fasta_prbh(args.fasta_file),
                   key=lambda rec: rec.prbh.score, reverse=True)
    output = open(args.output,"w+") if args.output else None
    print("FASTA Description\tBH start\tBH end\tPRBH score",file=output)
    for prbh in prbhs:
        print(format_prbh(prbh),file=output)