import sys
import os
import argparse
import random

progs = {
    "ssw":"bin/ssw_test",
    "ssearch36":"bin/ssearch36",
    "swipe":"bin/swipe",
    "swalign":"bin/swalign.x",
    "seqan":"bin/align.x",
    "diagonalsw":"bin/diagonalsw",
    "parasail":"bin/parasail_aligner"
    }

scoremat = "BLOSUM62.txt" # NOTE: parasail matrix option must be specified

def parse_opts():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input fasta file (query and target are generated from it)", type=str, required=True)
    #parser.add_argument("-t", help="input target file (fasta)", type=str)
    #parser.add_argument("-q", help="input query file (fasta)", type=str)
    parser.add_argument("-r", help="number of repeats to perform", type=int, default=1)
    parser.add_argument("-x", help="comparison, formatted as [query]x[target]", type=str, required=True)
    parser.add_argument("-fast", help="skip non-vectorized benchmarks", action="store_true")
    parser.add_argument("-resample", help="skip non-vectorized benchmarks", action="store_true")
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
    used = []
    for i in range(num):
        roll = random.randint(0, len(keys) - 1)
        if roll in used:
            while True:
                roll = random.randint(0, len(keys) - 1)
                if roll not in used:
                    break
        seqs[keys[roll]] = fasta[keys[roll]]
        used.append(roll)
    return seqs

def write_fasta(fasta, fname):
    with open(fname, "w") as f:
        for i in fasta.keys():
            f.write(">{key}\n{val}\n".format(key=i, val=fasta[i]))
    return fname

def get_call(prog, target, query, out):
    call = ""
    if prog == "ssw":
        call = "time {loc} -c -p -a {mat} {t} {q} > {o}".format(loc=progs[prog], mat=scoremat,
                                                             t=target, q=query, o=out)
    elif prog == "ssearch36":
        call = "time {loc} -d 99999 -T 1 -f '-11' -g '-1' -s {mat} {q} {t} > {o}".format(loc=progs[prog],
                                                                                        mat=scoremat,
                                                                                t=target, q=query, o=out)
        #consider -E 9999
    elif prog == "swipe":
        call = "time {loc} -a 1 -b 99999 -M {mat} -G 11 -E 1 -i {q} -d {t} > {o}".format(loc=progs[prog],
                                                                           mat=scoremat,
                                                                           t=target,
                                                                           q=query, o=out)
    elif prog == "swalign":
        call = "time {loc} -h -t 1 -s {mat} -q {q} -d {t} > {o}".format(loc=progs[prog], mat=scoremat,
                                                             t=target, q=query, o=out)
    elif prog == "seqan":
        call = "time {loc} {t} {q} > {o}".format(loc=progs[prog], t=target, q=query, o=out)
    elif prog == "diagonalsw":
        call = "time {loc} -s {mat} -d {t} -q {q} -i '-11' -e '-1' -u -t 1 > {o}".format(loc=progs[prog],
                                                                           mat=scoremat,
                                                                           t=target,
                                                                           q=query, o=out)
    elif prog == "parasail":
        call = "time {loc} -d -c 1 -m {mat} -f {t} -q {q} -o 11 -e 1 -t 1 -g /dev/stdout > {o}".format(loc=progs[prog],
                                                                           mat=scoremat,
                                                                           t=target,
                                                                           q=query, o=out)
    else:
        assert(False)
    return call

def main():
    #clean old files
    filelist_sup = [ f for f in os.listdir(os.path.abspath("sup/")) if f != "placeholder" ]
    filelist_out = [ f for f in os.listdir(os.path.abspath("out/")) if f != "placeholder" ]
    for f in filelist_sup:
        os.remove(os.path.abspath("sup/{ff}".format(ff=f)))
    for f in filelist_out:
        os.remove(os.path.abspath("out/{ff}".format(ff=f)))

    opts = parse_opts()

    comp = parse_comparison(opts.x)

    if opts.i:
        assert(os.path.exists(opts.i))
        parsed = parse_fasta(opts.i)
        infile = clean_fasta(parsed)

    if not opts.resample:
        query = write_fasta(split_fasta(infile, comp[0]), "sup/query_{suff}".format(suff=opts.i))
        targ = write_fasta(split_fasta(infile, comp[1]), "sup/target_{suff}".format(suff=opts.i))
        os.system("makeblastdb -in {t} -dbtype prot > /dev/null".format(t=os.path.abspath(targ)))

    calls = {}
    for i in range(opts.r):
        if opts.resample:
            query = write_fasta(split_fasta(infile, comp[0]), "sup/query_{d}_{suff}".format(d=i, suff=opts.i))
            targ = write_fasta(split_fasta(infile, comp[1]), "sup/target_{d}_{suff}".format(d=i, suff=opts.i))
            os.system("makeblastdb -in {t} -dbtype prot > /dev/null".format(t=os.path.abspath(targ)))

        for program in progs:
            if (program == "seqan" and opts.fast) or (program == "swalign" and opts.fast):
                continue
            if (program == "diagonalsw" and comp[0] != 1):
                continue
            calls.setdefault(program, [])
            # calls[program].append(get_call(program, targ, query, "/dev/null"))
            calls[program].append(get_call(program, targ, query, "out/{p}_{d}".format(p=program, d=i)))
    with open("calls", "w") as f:
        for i in calls:
            f.write("{{ {exp}; }} 2> {p} \n\n".format(exp=" && ".join(calls[i]), p=i))
        f.write("python graphs.py {i} {op}".format(i=" ".join(list(calls)), op=opts.x))


if __name__ == '__main__':
    main()
