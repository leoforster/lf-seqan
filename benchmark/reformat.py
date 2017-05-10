import numpy as np
import matplotlib.pyplot as plt
import sys
import os


def reset(d):
    for i in d:
        d[i] = None

def reformat_ssw(path):
    with open(floc + path, "r") as f:
        data = f.read()#.replace('\n', '')
        
    data = data.split("target_name")
    
    if outfmt == "seed":
        with open(floc + path + "_formatted", "w") as w:
            w.write("# Fields: s.len, s.seqnum, s.start, strand, "\
                    "q.len, q.seqnum, q.start, score\n")
            for i in data:
                reset(fields)
                if i == "":
                    continue
                i = i.split("\n")
                tq = i[0:2]
                rest = i[2].split("\t")
                                
                fields["s.seqnum"] = int(tq[0][2:])
                fields["q.seqnum"] = int(tq[1].split(":")[1])
                
                fields["score"] = int(rest[0].split(":")[1])
                if rest[1].startswith("suboptimal"):
                    fields["strand"] = 'P' if rest[2].split(":")[1] == "-" else 'F'
                    fields["s.start"] = int(rest[3].split(":")[1])
                    fields["s.len"] = int(rest[4].split(":")[1]) - fields["s.start"] + 1
                    fields["q.start"] = int(rest[5].split(":")[1])
                    fields["q.len"] = int(rest[6].split(":")[1]) - fields["q.start"] + 1
                else:
                    fields["strand"] = 'P' if rest[1].split(":")[1] == "-" else 'F'
                    fields["s.start"] = int(rest[2].split(":")[1])
                    fields["s.len"] = int(rest[3].split(":")[1]) - fields["s.start"] + 1
                    fields["q.start"] = int(rest[4].split(":")[1])
                    fields["q.len"] = int(rest[5].split(":")[1]) - fields["q.start"] + 1
                
                w.write("{sl}\t{ss}\t{st}\t{s}\t{ql}\t{qs}\t{qt}\t{c}\n".format(
                    sl = fields["s.len"],
                    ss = fields["s.seqnum"],
                    st = fields["s.start"],
                    s = fields["strand"],
                    ql = fields["q.len"],
                    qs = fields["q.seqnum"],
                    qt = fields["q.start"],
                    c = fields["score"]))
                
    if outfmt == "align":
        print("TODO: implement align output")
        
def reformat_ssearch(path):
    with open(floc + path, "r") as f:
        data = f.read()#.replace('\n', '')
        
    data = data.split(">>")
    out = [0]
    
    for n, i in enumerate(data):
        if i.startswith(">"):
            out.append(n)
    for i in out[::-1]:
        del data[i]
    
    if outfmt == "seed":
        with open(floc + path + "_formatted", "w") as w:
            w.write("# Fields: s.len, s.seqnum, s.start, strand, "\
                    "q.len, q.seqnum, q.start, score, identity\n")
            for i in data:
                reset(fields)
                i = i.split("\n")
                
                rest = i[2].split("overlap (")[1][:-1].split("-")
                
                fields["q.seqnum"] = int(i[5].split(" ")[0])
                fields["s.seqnum"] = int(i[7].split(" ")[0])
                fields["score"] = int(i[2].split(";")[0].split(":")[1])
                fields["strand"] = 'F'
                fields["q.start"] = int(rest[0])
                fields["q.len"] = int(rest[1].split(":")[0]) - fields["q.start"] + 1
                fields["s.start"] = int(rest[1].split(":")[1])
                fields["s.len"] = int(rest[2]) - fields["s.start"] + 1
                fields["identity"] = float(i[2].split(";")[1].split("%")[0])
                
                w.write("{sl}\t{ss}\t{st}\t{s}\t{ql}\t{qs}\t{qt}\t{c}\t{i}\n".format(
                        sl = fields["s.len"],
                        ss = fields["s.seqnum"],
                        st = fields["s.start"],
                        s = fields["strand"],
                        ql = fields["q.len"],
                        qs = fields["q.seqnum"],
                        qt = fields["q.start"],
                        c = fields["score"],
                        i = fields["identity"]))

    if outfmt == "align":
        print("TODO: implement align output")
        
def reformat_seqan(path):
    with open(floc + path, "r") as f:
        data = f.read()
        
    data = data.split("Score = ")
    del data[0]
    
    if outfmt == "seed":
        with open(floc + path + "_formatted", "w") as w:
            w.write("# Fields: s.len, s.seqnum, s.start, strand, "\
                    "q.len, q.seqnum, q.start, score\n")
            for i in data:
                reset(fields)
                i = i.split("\n")
                rest = i[-3].split("to")
                
                fields["s.seqnum"] = int(rest[0].split("[")[0])
                fields["q.seqnum"] = int(rest[1].split("[")[0])
                fields["score"] = int(i[0])
                fields["strand"] = 'F'
                fields["s.start"] = int(rest[0].split("[")[1].split(":")[0])
                fields["s.len"] = int(rest[0].split("[")[1].split(":")[1][:-2]) - fields["s.start"] + 1
                fields["q.start"] = int(rest[1].split("[")[1].split(":")[0])
                fields["q.len"] = int(rest[1].split("[")[1].split(":")[1][:-1]) - fields["q.start"] + 1
                
                w.write("{sl}\t{ss}\t{st}\t{s}\t{ql}\t{qs}\t{qt}\t{c}\n".format(
                        sl = fields["s.len"],
                        ss = fields["s.seqnum"],
                        st = fields["s.start"],
                        s = fields["strand"],
                        ql = fields["q.len"],
                        qs = fields["q.seqnum"],
                        qt = fields["q.start"],
                        c = fields["score"]))

    if outfmt == "align":
        print("TODO: implement align output")

def reformat_swipe(path):
    with open(floc + path, "r") as f:
        data = f.read()
    
    if outfmt == "seed":
        with open(floc + path + "_formatted", "w") as w:
            w.write("# Fields: s.len, s.seqnum, s.start, strand, "\
                    "q.len, q.seqnum, q.start, score, identity\n")
            q = data.split("Database file")
            del q[0]
            for subq in q:
                query = int(subq.split("Query description")[1][:100].split("\n")[0][1:])

                run = subq.split(">")
                del run[0]
                for i in run:
                    reset(fields)
                    i = i.split("\n")
                                    
                    fields["s.seqnum"] = int(i[0][i[0].index(" "):])
                    fields["q.seqnum"] = query
                    fields["score"] = int(i[3].split("=")[1])
                    fields["strand"] = 'F'
                    fields["s.start"] = int("".join([x for x in i[8].split(":")[1][:10] if x.isdigit()]))
                    fields["q.start"] = int("".join([x for x in i[6].split(":")[1][:10] if x.isdigit()]))
                    try:
                        fields["q.len"] = int(i[-6].split(" ")[-1]) - fields["q.start"] + 1
                        fields["s.len"] = int(i[-4].split(" ")[-1]) - fields["s.start"] + 1
                    except:
                        fields["q.len"] = int(i[-5].split(" ")[-1]) - fields["q.start"] + 1
                        fields["s.len"] = int(i[-3].split(" ")[-1]) - fields["s.start"] + 1
                    fields["identity"] = int(i[4].split("%")[0].split("(")[1])
                    
                    w.write("{sl}\t{ss}\t{st}\t{s}\t{ql}\t{qs}\t{qt}\t{c}\t{i}\n".format(
                            sl = fields["s.len"],
                            ss = fields["s.seqnum"],
                            st = fields["s.start"],
                            s = fields["strand"],
                            ql = fields["q.len"],
                            qs = fields["q.seqnum"],
                            qt = fields["q.start"],
                            c = fields["score"],
                            i = fields["identity"]))
                
    if outfmt == "align":
        print("TODO: implement align output")
        
def reformat_swalign(path):
    with open(floc + path, "r") as f:
        data = f.read()
    
    data = data.split("\n")
    
    with open(floc + path + "_formatted", "w") as w:
        w.write("# Fields: s.len, s.seqnum, s.start, strand, "\
                "q.len, q.seqnum, q.start, score\n")
        for i in data:
            reset(fields)
            if i.startswith("#") or i == "":
                continue
            i = i.split("\t")
        
            fields["s.seqnum"] = int(i[0])
            fields["q.seqnum"] = int(i[1])
            fields["score"] = int(i[6].strip("\n"))
            fields["strand"] = 'F'
            fields["s.start"] = int(i[2])
            fields["s.len"] = int(i[3])
            fields["q.start"] = int(i[4])
            fields["q.len"] = int(i[5])
            
            w.write("{sl}\t{ss}\t{st}\t{s}\t{ql}\t{qs}\t{qt}\t{c}\n".format(
                    sl = fields["s.len"],
                    ss = fields["s.seqnum"],
                    st = fields["s.start"],
                    s = fields["strand"],
                    ql = fields["q.len"],
                    qs = fields["q.seqnum"],
                    qt = fields["q.start"],
                    c = fields["score"]))
            
    if outfmt == "align":
        pass

def main():
    files = [ f for f in os.listdir(os.path.abspath(floc)) if f != "placeholder" and
                                                            not f.startswith(".") and 
                                                            not f.endswith("_formatted")]
    print("attempting to reformat: ")
    print(files)
    
    for i in files:
        try:
            if "ssw" in i:
                reformat_ssw(i)
                pass
            elif "ssearch36" in i:
                reformat_ssearch(i)
                pass
            elif "seqan" in i:
                reformat_seqan(i)
                pass
            elif "swipe" in i:
                reformat_swipe(i)
                pass
            elif "swalign" in i:
                reformat_swalign(i)
                pass
            else:
                print("skipping", i)
        except:
            print("error in ", i)
        
fields = {
    "s.len":None,
    "s.seqnum":None,
    "s.start":None,
    "strand":None,
    "q.len":None,
    "q.start":None,
    "q.seqnum":None,
    "score":None,
    "edist":None,
    "identity":None
    }

if sys.argv[1] == "-seed" or sys.argv[1] == "-s":
    outfmt = "seed"
elif sys.argv[1] == "-swalign" or sys.argv[1] == "-a":
    outfmt = "align"
else:
    print("first argument should be -seed/-s or -swalign/-a")
    sys.exit(1)

try:
    floc = sys.argv[2]
except:
    print("second argument is location of output files")
    sys.exit(1)
        
if __name__ == '__main__':
    main()