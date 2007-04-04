#!/usr/bin/python

def loadup(fname) :
    f = file(fname)
    r= " "
    found = False
    while(found==False and r!="") :
        r = f.readline()
       
        try :
            if r.index("TIME")==0 :
                found = True
        except ValueError :
            continue
       
    headings = r.split()
    if headings[1]=="[YR]" :
        del headings[1]
    tindex = 0
    vindex = headings.index("M_V")
    vrindex = headings.index("(V-R)")
    viindex = headings.index("(V-I)")
    bvindex = headings.index("(B-V)")
    vkindex = headings.index("(V-K)")
    
    time = []
    b_m = []
    i_m = []
    r_m = []
    v_m = []
    k_m = []
    
    while(r!="") :
        r = f.readline()
        vs = r.split()
        try:
            time.append(float(vs[tindex]))
            b_m.append(float(vs[vindex])+float(vs[bvindex]))
            r_m.append(float(vs[vindex])-float(vs[vrindex]))
            i_m.append(float(vs[vindex])-float(vs[viindex]))
            v_m.append(float(vs[vindex]))
            k_m.append(float(vs[vindex])-float(vs[vkindex]))
        except IndexError :
            continue
        except ValueError :
            continue
        
    return time,b_m,i_m,r_m,v_m,k_m

def writefile(fname,time,metals,vals) :
    f = file(fname,'w')
    f.write("#METAL ")
    for i in range(0,len(metals)) :
        f.write(str(metals[i][0])+" ")
    f.write("\n")

    for t in xrange(0,len(time)) :
        f.write(str(time[t])+" ")
        for i in range(0,len(metals)) :
            f.write(str(vals[i][t])+" ")
        f.write("\n")
    
import glob, operator;

files = glob.glob("sb99_z*")
metals = []

for f in files :
    g = f.split("_z")
    try:
        metals.append([float(g[1]),f])
    except ValueError:
        continue
    
b = []
i = []
r = []
v = []
k = []
time_loc = []

metals.sort(key=operator.itemgetter(1))

for n in range(0,len(metals)) :
    time,b_loc,i_loc,r_loc,v_loc,k_loc = loadup(metals[n][1])
    if time!=time_loc and time_loc!=[] :
        raise "TIME MISMATCH",metals[n][1]
        
    time_loc=time
    b.append(b_loc)
    i.append(i_loc)
    r.append(r_loc)
    v.append(v_loc)
    k.append(k_loc)

writefile("LUM_B",time_loc,metals,b)
writefile("LUM_I",time_loc,metals,i)
writefile("LUM_R",time_loc,metals,r)
writefile("LUM_V",time_loc,metals,v)
writefile("LUM_K",time_loc,metals,k)
