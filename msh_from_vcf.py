import sys
import numpy as np
import argparse


def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen",dest="genname")
    parser.add_argument("--sub",dest="subname")
    parser.add_argument("--gen-idx",dest="genidx",type=int,default=0)
    parser.add_argument("--squish",dest="squish",action="store_true")
    return parser

def splitAllelesAll(la):
    alleles = []
    ac = 0
    an = 0
    for i in range(9,len(la)):
        try:
            g1 = int(la[i][0])
            g2 = int(la[i][2])
        except ValueError:
            return []
        alleles.append(g1)
        alleles.append(g2)
        ac += g1
        ac += g2
        an += 2
    if ac == 1 or ac == an - 1:
        return []
    return alleles

def splitAllelesSub(la,idx_list):
    alleles = []
    ac = 0
    an = 0
    for i in range(len(idx_list)):
        reg = (idx_list[i]%2)*2
        f_idx = (idx_list[i]//2)+9
        try:
            geno = int(la[f_idx][reg])
        except ValueError:
            return []
        alleles.append(geno)
        ac += geno
        an += 1
    if ac <= 1 or ac >= an - 1:
        return []
    return alleles

def splitAlleles(la,idx_list=None):
    if idx_list is None:
        return splitAllelesAll(la)
    else:
        return splitAllelesSub(la,idx_list)

def getVectors(a_prev,d_prev,cur_row,k):
    p = k+1
    q = k+1
    a = []
    d = []
    b = []
    e = []
    for i in range(len(a_prev)):
        p = max(d_prev[i],p)
        q = max(d_prev[i],q)
        if cur_row[a_prev[i]] == 0:
            a.append(a_prev[i])
            d.append(p)
            p = 0
        else:
            b.append(a_prev[i])
            e.append(q)
            q = 0
    a.extend(b)
    d.extend(e)
    return a,d

def msh(a,d,k,pos_list):
    l = len(a)
    y_msh = [0 for i in range(l)]
    site_msh = [0 for i in range(l)]
    for i in range(l):
        if i == l-1:
            c_idx = d[l-1]
            #y_msh[i] = abs(pos_list[-1] - pos_list[d[l-1]])
        elif i == 0:
            c_idx = d[1]
            #y_msh[i] = abs(pos_list[-1] - pos_list[d[1]])
        else:
            c_idx = min(d[i],d[i+1])
            #y_msh[i] = abs(pos_list[-1] - pos_list[min(d[i],d[i+1])])
        if c_idx == 0:
            y_msh[i] = -2
        elif c_idx != len(pos_list) - 1:
            y_msh[i] = abs(pos_list[-1] - pos_list[c_idx])
        else:
            y_msh[i] = abs(pos_list[-1] - pos_list[-2])
    for a_i,a_v in enumerate(a):
        site_msh[a_v] = y_msh[a_i]
    return site_msh

def msh_gen(a,d,k,gen_list):
    l = len(a)
    g_msh = [0.0 for i in range(l)]
    site_msh = [0.0 for i in range(l)]
    for i in range(l):
        if i == l-1:
            c_idx = d[l-1]
            #g_msh[i] = abs(gen_list[-1] - gen_list[d[l-1]])
        elif i == 0:
            c_idx = d[1]
            #g_msh[i] = abs(gen_list[-1] - gen_list[d[1]])
        else:
            c_idx = min(d[i],d[i+1])
            #g_msh[i] = abs(gen_list[-1] - gen_list[min(d[i],d[i+1])])
        if c_idx == 0:
            g_msh[i] = -2
        elif c_idx != len(pos_list) - 1:
            g_msh[i] = abs(gen_list[-1] - gen_list[c_idx])
        else:
            g_msh[i] = abs(gen_list[-1] - gen_list[-2])
    for a_i,a_v in enumerate(a):
        site_msh[a_v] = g_msh[a_i]
    return site_msh

def parseGenLine(l,offset):
    la = l.strip().split()
    a = float(la[0])
    b = float(la[1+offset])
    return a,b

def getGenMap(f,idx=0,squish=False):
    l1 = []
    l2 = []
    for line in f:
        if line[0] not in [0,1,2,3,4,5,6,7,8,9]:
            continue
        a,b = parseGenLine(line,idx)
        #if b != 0:
        if not squish or (len(l2) > 0 and l2[-1] != b):
            l1.append(a)
            l2.append(b)
    return l1,l2

def getGenPos(pos,l1,l2,cm=True):
    rate = 1.0
    if cm:
        rate = .01
    if float(pos) < l1[0]:
        gr = (l2[1]-l2[0])*rate/(l1[1]-l1[0])
        return gr*(pos-l1[0])
    i = 1
    while i < len(l1)-1 and float(pos) > float(l1[i]):
        i += 1
    gr = (l2[i]-l2[i-1])*rate/float(l1[i]-l1[i-1])
    return rate*l2[i-1]+gr*float(pos-l1[i-1])

parser = createParser()
args = parser.parse_args()

gen_flag = False
sub_flag = False
if args.genname is not None:
    gf = open(args.genname,'r')
    l1,l2 = getGenMap(gf,idx=args.genidx,squish=args.squish)
    gen_flag = True
    sys.stderr.write(str(len(l1))+'\n')

idx_list = None
if args.subname is not None:
    sf = open(args.subname,'r')
    idx_list = [int(i.strip()) for i in sf]
    sub_flag = True

a_mat = []
d_mat = []
m_mat = []

k = 0
pos_list = []
gen_list = None
if gen_flag:
    gen_list = []
for line in sys.stdin:
    if line[0] == '#':
        continue
    la = line.strip().split()
    if 'CPG_TAG' in la[7]:
        continue
    if len(la[3]) != 1 or len(la[4]) != 1:
        continue
    alleles = splitAlleles(la,idx_list)
    if len(alleles) == 0:
        continue
    if k == 0:
        sample_count = len(alleles)
        sys.stderr.write(str(sample_count)+'\n')
        a_prev = [i for i in range(sample_count)]
        d_prev = [0 for i in range(sample_count)]
    #sys.stderr.write(str(k)+'\n')
    pos_list.append(int(la[1]))
    if gen_flag:
        gen_list.append(float(getGenPos(int(la[1]),l1,l2)))
    a,d = getVectors(a_prev,d_prev,alleles,k)
    msh_vec = msh(a,d,k+1,pos_list)
    if gen_flag:
        g_vec = msh(a,d,k+1,gen_list)
    #tv = [str(ll) for ll in sorted(msh_vec)]
    sys.stdout.write(str(pos_list[-1]))
    if gen_flag:
        sys.stdout.write('\t'+str(gen_list[-1]))
    for i in range(len(msh_vec)):
        sys.stdout.write('\t'+str(msh_vec[i]))
        if gen_flag:
            sys.stdout.write(':'+str(g_vec[i]))
    sys.stdout.write('\n')

    a_prev = a
    d_prev = d
    k += 1
