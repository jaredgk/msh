from AlleleAgeEstimator import datalist,data


def intw(v):
    if '*' in v:
        return None
    if v == '-2':
        return None
    return int(v)

def floatw(v):
    if v is None or '*' in v or '-' in v:
        return None
    if v == '-2':
        return None
    return float(v)

def getDataLists(tl,ac,rec,mu,mod):
    dl1 = datalist()
    dl2 = datalist()
    al = datalist()
    fl = []
    for i in xrange(ac):
        la = tl[i].strip().split()
        fl.append(la)
    pos = int(fl[0][1])
    for i in xrange(ac):
        la = fl[i]
        d1 = intw(la[2])
        d2 = intw(la[3])
        if d1 != None and d2 != None:
            d = data(dis1 = d1, dis2 = d2, rho1 = rec, rho2 = rec, mu=mu, model=mod)
            dl1.append(d)
        elif d1 != None:
            d = data(dis1 = d1, rho1 = rec, mu=mu, model=mod)
            dl1.append(d)
        elif d2 != None:
            d = data(dis1 = d2, rho1 = rec, mu=mu,model=mod)
            dl1.append(d)
        d1 = intw(la[4])
        d2 = intw(la[5])
        if d1 != 0 or d2 != 0:
            if d1 is not None and d2 is not None:
                d = data(dis1 = d1, dis2 = d2, rho1 = rec, rho2 = rec, mu=mu,model=mod)
                dl2.append(d)
            elif d1 is not None:
                d = data(dis1 = d1, rho1 = rec, mu=mu,model=mod)
                dl2.append(d)
            elif d2 != None:
                d = data(dis1 = d2, rho1 = rec, mu=mu,model=mod)
                dl2.append(d)


        if ac >= 2:
            for j in xrange(i):
                d1 = intw(la[6+j])
                d2 = intw(la[6+ac+j])
                if d1 != None and d2 != None:
                    a = data(dis1 = d1, dis2 = d2, rho1 = rec, rho2 = rec, mu=mu, model=mod)
                    al.append(a)
                elif d1 != None:
                    a = data(dis1 = d1, rho1 = rec, mu=mu,model=mod)
                    al.append(a)
                elif d2 != None:
                    a = data(dis1 = d2, rho1 = rec, mu=mu,model=mod)
                    al.append(a)
    return dl1,dl2,al,pos

def p(val,f):
    f.write(str(val)+'\t')	
