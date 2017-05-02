import sys
import os
import argparse
import random

progs = {
    "ssw":"bin/ssw_test",
    "ssearch36":"bin/ssearch36",
    "swipe":"bin/swipe",
    "swalign":"bin/swalign.x",
    "seqan":"bin/align.x"
    }

scoremat = "BLOSUM62.txt"

def parse_opts():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input fasta file (query and target are generated from it)", type=str, required=True)
    #parser.add_argument("-t", help="input target file (fasta)", type=str)
    #parser.add_argument("-q", help="input query file (fasta)", type=str)
    parser.add_argument("-r", help="number of repeats to perform", type=int, default=1)
    parser.add_argument("-x", help="comparison, formatted as [query]x[target]", type=str, required=True)
    parser.add_argument("-fast", help="skip non-vectorized benchmarks", action="store_true")
    args = parser.parse_args()
    return args

def parse_comparison(comp):
    comp = comp.split("x")
    return (int(comp[0]), int(comp[1]))

def parse_fasta(infile):
    seqs = {}
    header = -1
    with open(infile, "r") as f:
        for i, n in enumerate(f.readlines()):
            if n.startswith(">"):
                header += 1
                seqs.setdefault(header, "")
                continue
            elif n == "\n":
                continue
            else:
                seqs[header] += n.strip("\n")
    return seqs

def clean_fasta(fasta):
    gone = []
    for i in fasta:
        if len(fasta[i]) < 30:
            gone.append(i)
    for x in gone:
        del fasta[x]
    return fasta

def split_fasta(fasta, num):
    seqs = {}
    keys = list(fasta)
    for i in range(num):
        roll = random.randint(0, len(keys) - 1)
        seqs[keys[roll]] = fasta[keys[roll]]
    return seqs

def write_fasta(fasta, fname):
    with open(fname, "w") as f:
        for i in fasta.keys():
            f.write(">{key}\n{val}\n".format(key=i, val=fasta[i]))
    return fname

def get_call(prog, target, query, out):
    call = ""
    if prog == "ssw":
        call = "time {loc} -p -a {mat} {t} {q} > {o}".format(loc=progs[prog], mat=scoremat, 
                                                             t=target, q=query, o=out)
    elif prog == "ssearch36":
        call = "time {loc} -s {mat} {q} {t} > {o}".format(loc=progs[prog], mat=scoremat, 
                                                             t=target, q=query, o=out)
    elif prog == "swipe":
        call = "time {loc} -M {mat} -G 11 -E 1 -i {q} -d {t} > {o}".format(loc=progs[prog], 
                                                                           mat=scoremat, 
                                                                           t=target, 
                                                                           q=query, o=out)
    elif prog == "swalign":
        call = "time {loc} -s {mat} -q {q} -d {t} > {o}".format(loc=progs[prog], mat=scoremat, 
                                                             t=target, q=query, o=out)
    elif prog == "seqan":
        call = "time {loc} {t} {q} > {o}".format(loc=progs[prog], t=target, q=query, o=out)
    else:
        assert(False)
    return call

def main():
    opts = parse_opts()

    comp = parse_comparison(opts.x)

    if opts.i:
        assert(os.path.exists(opts.i))
        infile = clean_fasta(parse_fasta(opts.i))
    #else:
        #if not opts.t and not opts.q:
            #print("if not using -i, -t and -q are required")
            #sys.exit(1)
        #assert(os.path.exists(opts.t))
        #assert(os.path.exists(opts.q))
        #targ = opts.t
        #query = opts.q
        
    calls = {}
    for i in range(opts.r):
        query = write_fasta(split_fasta(infile, comp[0]), "sup/query_{d}_{suff}".format(d=i, suff=opts.i))
        targ = write_fasta(split_fasta(infile, comp[1]), "sup/target_{d}_{suff}".format(d=i, suff=opts.i))
        os.system("makeblastdb -in {t} -dbtype prot > /dev/null".format(t=os.path.abspath(targ)))
        for program in progs:
            if (program == "seqan" and opts.fast) or (program == "swalign" and opts.fast):
                continue
            calls.setdefault(program, [])
            calls[program].append(get_call(program, targ, query, "out/{p}_{d}".format(p=program, d=i)))
    
    with open("calls", "w") as f:
        for i in calls:
            f.write("{{ {exp}; }} 2> {p} \n\n".format(exp=" && ".join(calls[i]), p=i))
        f.write("python graphs.py {i} {op}".format(i=" ".join(list(calls)), op=opts.x))
    

if __name__ == '__main__':
    main()




























