import sys
from AlleleAgeEstimator import data,datalist,fitmodel
from aae_parse import p, intw, floatw
import numpy as np
import gzip
import argparse
import math

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("left_file")
    parser.add_argument("right_file")
    parser.add_argument("--start",dest="start",type=int,default=-1)
    parser.add_argument("--end",dest="end",type=int,default=-1)
    parser.add_argument("--n",dest="n_model",type=int,default=100)
    parser.add_argument("--n0",dest="n0_model",type=int,default=1000)
    parser.add_argument("--mut",dest="mut_rate",type=float,default=1e-8)
    parser.add_argument("--rec",dest="rec_rate",type=float,default=1e-8)
    parser.add_argument("--no-cache",dest="cache_est",action="store_false")
    parser.add_argument("--mod-gen",dest="mod_gen",action="store_true")
    parser.add_argument("--full-output",dest="full_out",action="store_true")
    return parser

def genFunction(val):
    return float(1-math.exp(-2*val))/2

def makeData(pos,l1,l2,gen_flag,rec,mu,mod,mod_gen=False):
    d = None
    if gen_flag:
        d1 = intw(l1.split(':')[0])
        d2 = intw(l2.split(':')[0])
        m1 = floatw(l1.split(':')[1])
        m2 = floatw(l2.split(':')[1])
    else:
        d1 = intw(l1)
        d2 = intw(l2)
        m1 = (d1*rec if d1 is not None else None)
        m2 = (d2*rec if d2 is not None else None)
    if mod_gen and m1 is not None:
        m1 = genFunction(m1)
    if mod_gen and m2 is not None:
        m2 = genFunction(m2)
    if d1 is None and d2 is None:
        return d
    d = data(dis1 = d1,dis2 = d2, morgans1 = m1, morgans2 = m2, rho1 = rec, rho2 = rec, mu=mu, model=mod)
    #if d1 != None and d2 != None:
    #    d = data(dis1 = d1, dis2 = d2, morgans1 = m1, morgans2 = m2, rho1 = rec, rho2 = rec, mu=mu, model=mod)
    #elif d1 != None:
    #    d = data(dis1 = d1, morgans1 = m1, rho1 = rec, mu=mu, model=mod)
    #elif d2 != None:
    #    d = data(dis1 = d2, morgans1 = m2, rho1 = rec, mu=mu,model=mod)
    return d

def getDataLists(pos,l1,l2,gen_flag,rec,mu,mod):
    dl1 = datalist()
    #pos = int(la[1])
    ann = ''
    #d1 = intw(la[2])
    #d2 = intw(la[3])
    #m1 = floatw(la[4])
    #m2 = floatw(la[5])
    #if len(la) >= 7:
    #	ann = la[6]
    if gen_flag:
        d1 = intw(l1.split(':')[0])
        d2 = intw(l2.split(':')[0])
        m1 = floatw(l1.split(':')[1])
        m2 = floatw(l2.split(':')[1])
    else:
        d1 = intw(l1)
        d2 = intw(l2)
        m1 = floatw(d1)*rec
        m2 = floatw(d2)*rec

    if d1 != None and d2 != None:
        #m1 = float(d1)*rec
        #m2 = float(d2)*rec
        d = data(dis1 = d1, dis2 = d2, morgans1 = m1, morgans2 = m2, rho1 = rec, rho2 = rec, mu=mu, model=mod)
        dl1.append(d)
    elif d1 != None:
        #m1 = float(d1)*rec
        d = data(dis1 = d1, morgans1 = m1, rho1 = rec, mu=mu, model=mod)
        dl1.append(d)
    elif d2 != None:
        #m2 = float(d2)*rec
        d = data(dis1 = d2, morgans1 = m2, rho1 = rec, mu=mu,model=mod)
        dl1.append(d)

    return dl1


def parseFilename(fn):
    fwd = fn.split('/')[-1]
    nec = 0
    n = 0
    ac = int(fwd.split('.')[-1])
    return ac

def changeModel(dl,m):
    for d in dl:
        d.model = m

def parseGenLine(l,offset):
    la = line.strip().split()
    a = float(la[0])
    b = float(la[1+offset])
    return a,b

def getGenMap(f):
    l1 = []
    l2 = []
    for line in f:
        a,b = parseGenLine(line,0)
        if b != 0 and (len(l1) < 2 or b != l1[-1]):
            l1.append(a)
            l2.append(b)
    return l1,l2

def getGenRate(pos,l1,l2):
    if pos < l1[0]:
        return (l2[1]-l2[0])*.01/(l1[1]-l1[0])
    i = 1
    while i < len(l1)-1 and pos < l1[i]:
        i += 1
    return (l2[i]-l2[i-1])*.01/(l1[i]-l1[i-1])

def geomean(l):
    ll = [float(1 if iiii == 0 else iiii) for iiii in l]
    b = np.array(ll).astype(np.float)
    a = np.log(b)
    return np.exp(a.sum()/len(a))

def fullStr(dl):
    d = dl[0]
    left_side = ''
    if d.side != 2:
        left_side = str(d.dis1)+','+str(d.morgans1)
    else:
        left_side = '-1,-1'
    if d.side != 1:
        right_side = str(d.dis2)+','+str(d.morgans2)
    else:
        right_side = '-1,-1'
    return left_side+','+right_side+','+str(dl.tcest)

parser = createParser()
args = parser.parse_args()

#fnl = str(sys.argv[1])
if args.left_file[-3:] == '.gz':
    fl = gzip.open(args.left_file,'r')
else:
    fl = open(args.left_file,'r')



#fnr = str(sys.argv[2])
if args.right_file[-3:] == '.gz':
    fr = gzip.open(args.right_file,'r')
else:
    fr = open(args.right_file,'r')


start_idx = args.start
end_idx = args.end



mu = args.mut_rate
rec = args.rec_rate

n_model = args.n_model
n0_model = args.n0_model

#fl = f.readlines()
#tot = len(fl)

header="pos\testimate-ml\testimate-bayes\n"

#outf.write(header)
ii = 0
has_genetic_positions = False
dl1 = datalist()
for l1 in fl:

    l2 = fr.readline()
    if start_idx != -1 and ii < start_idx:
        ii += 1
        continue
    elif start_idx != -1 and ii >= end_idx:
        break
    la1 = l1.strip().split()
    la2 = l2.strip().split()
    est_list_ml = []
    est_list_bayes = []
    if la1[0] != la2[0]:
        raise Exception('positions dont match')
    pos = int(la1[0])
    if not has_genetic_positions and ':' in la1[2]:
        has_genetic_positions = True
    start_inds = 1
    if has_genetic_positions:
        gen = float(la1[1])
        start_inds = 2
    sys.stdout.write(str(pos))
    for i in range(start_inds,len(la1)):
        mc = fitmodel(n=n_model,N0=n0_model,popmodel="c")
        #dl1 = getDataLists(pos,la1[i],la2[i],has_genetic_positions,rec,mu,mc)
        d = makeData(pos,la1[i],la2[i],has_genetic_positions,rec,mu,mc,mod_gen=args.mod_gen)
        if d is None:
            if args.full_out:
                sys.stdout.write('\t-1,-1,-1,-1,-1')
            continue
        if len(dl1) == 0:
            dl1.append(d)
        else:
            dl1[0] = d
        if len(dl1) != 0:
            dl1.estimate_tc_cache(cache_est=args.cache_est)
            est_list_ml.append(dl1.tcest)
            if args.full_out:
                sys.stdout.write('\t'+fullStr(dl1))
    est_all_ml = geomean(est_list_ml)
    #if not args.full_out:
    sys.stdout.write('\t'+str(est_all_ml))
    sys.stdout.write('\n')
    ii += 1
