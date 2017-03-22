import sys
import os
import subprocess
import argparse

#todo: mincov in seqan
#todo: -maxmat and -l -> maximal matches of minimum length l
#todo: chaining
#todo: input an encoded file

#todo: ouput option
#todo: seedfile option
#todo: option to compare successful seeds, failed seeds, or both
def parse_opts():
  parser = argparse.ArgumentParser()
  parser.add_argument("-s", help="compiled seqan script location", type=str, required=True)
  parser.add_argument("-i", help="input fasta file", type=str, required=True)
  parser.add_argument("-q", help="input query fasta file for multiple comparison", type=str)
  parser.add_argument("-seedlen", help="seedlength", type=int)
  parser.add_argument("-mincov", help="minimum coverage for seed_extend", type=int)
  parser.add_argument("-x", help="use xdrop algorithm for seed extension", action="store_true")
  parser.add_argument("-g", help="use greedy algorithm for seed extension", action="store_true")
  parser.add_argument("-xcut", help="xdrop cutoff score", type=int)
  parser.add_argument("-minid", help="minimum identity of matches", type=int)
  args = parser.parse_args()
  return args
  
def check_opts(args):
  params = {
    "seqan":None,
    "infile":None,
    "qfile":None,
    "seedlen":None,
    "mincov":None,
    "minid":None,
    "do_xdrop":False,
    "do_greedy":False,
    "xcutoff":None 
    }
  
  #seqan
  if args.s[0] == "/":
    params["seqan"] = args.s
  else:
    params["seqan"] = os.getcwd() + "/" + args.s
  #check if executable?
  if not os.path.isfile(args.s):
    print("failed to find %s - does it exist?" %args.s)
    sys.exit()
  
  #infile
  if args.i[0] == "/":
    params["infile"] = args.i
  else:
    params["infile"] = os.getcwd() + "/" + args.i
  if not os.path.isfile(args.i):
    print("failed to find %s - does it exist?" %args.i)
    sys.exit()
      
  #query file
  if args.q:
    if args.q[0] == "/":
      params["qfile"] = args.q
    else:
      qfile = os.getcwd() + "/" + args.q
    if not os.path.isfile(args.q):
      print("failed to find %s - does it exist?" %args.q)
      sys.exit()

  params["seedlen"] = args.seedlen if args.seedlen else 8
  params["minid"] = args.minid if args.minid else 80
  params["do_xdrop"] = True if args.x else False
  params["do_greedy"] = True if args.g else False
  if args.x:
      params["xcutoff"] = args.xcut if args.xcut else 6

  #xdrop, greedy exclusion
  if args.x and args.g:
    print("cannot have both -x and -g")
    sys.exit()
    
  if args.mincov:
    params["mincov"] = args.mincov
  else:
    params["mincov"] = params["seedlen"] * 2.5 #float or int?

  return params
  
def sequences_from_fasta(infile):
  seqs = []
  idx = -1
  with open(infile, "r") as f:
    for i, n in enumerate(f.readlines()):
      if n.startswith(">"):
        idx += 1
        seqs.append("")
        continue
      elif n == "\n":
        continue
      else:
        seqs[idx] += n.strip("\n")
  return seqs

def do_gt_extend(params):
  pstr = ""
  
  if params["do_xdrop"]: pstr += " -extendxdrop 100" 
  if params["do_greedy"]: pstr += " -extendgreedy 100" 
  if params["qfile"] != None: pstr += " -qii %s" %(params["qfile"]) 
  if params["seedlen"]: pstr += " -seedlength %s" %(params["seedlen"])
  #if params["xcutoff"]: pstr += " -xdropbelow %s" %(params["xcutoff"])
  #if params["mincov"]: pstr += " -mincoverage %s" %(params["mincov"])
  #if params["minid"]: pstr += " -minidentity %s" %(params["minid"])
  #temporary:
  pstr += " -no-reverse"
  
  call = "gt seed_extend -ii %s %s -outfmt seed failed_seed" %(params["infile"], pstr)
  #print(call)
  try:
    seeds = subprocess.check_output([call], shell=True, stderr=subprocess.STDOUT)
  except:
    print("failed to get seeds -- call was:\n ", call)
    sys.exit()
  seeds = seeds.decode("utf-8").split("\n")
  return seeds
  
def run_with_seqan(seqan, infile, seedfile):
  call = "%s %s %s" %(seqan, infile, seedfile) #todo: how to handle other options in seqan-file?
  #print(call)
  try:
    extend = subprocess.check_output([call], shell=True, stderr=subprocess.STDOUT)
  except:
    print("failed to run seqan script -- call was:\n " %call)
    sys.exit()
  extend = extend.decode("utf-8").split("\n")[:-2] #slice due to irregularities
  return extend
  
def encode(infile):
  #encseq encode places encoded files in script calling directory -- change it?
  try:
    subprocess.call(["gt encseq encode %s" %infile], shell=True, stderr=subprocess.STDOUT)
  except:
    print("failed to encode %s" %infile)
    sys.exit()
    
def decode(infile):
  #todo
  pass
  
def parse_seeds(seeds):
  realseeds = []
  failseeds = []
  gt_extend = []
  for s in seeds:
    if s != "":
      if s.startswith("#"):
        if s.startswith("# failed_seed"):
          s = s.split("\t")
          failseeds.append(Seed(True, s[4], s[2], s[3], s[5], s[6], s[1]))
      else:
        gt_extend.append(ExtendedSeed.from_string(s))
        s = s.split(" ")
        realseeds.append(Seed(False, s[3], s[1], s[11], s[5], s[12], s[10]))
  return (realseeds, failseeds, gt_extend)
  
def seeds_to_file(seedfile, seeds):
  with open(seedfile, "w") as f:
    for s in seeds:
      #assert(type(s) == Seed)
      f.write(s.to_line())

class ExtendedSeed:
  def __init__(self, sl, s, spos, strand, ql, q, qpos, score, edist, ident):
    #check spos <= slen && qpos <= qlen ?
    self.slen = sl
    self.s = s
    self.spos = spos
    self.strand = strand
    self.qlen = ql
    self.q = q
    self.qpos = qpos
    self.score = score
    self.edist = edist
    self.ident = ident
    
  @classmethod
  def from_string(cls, line):
    sl, s, spos, strand, ql, q, qpos, score, edist, ident = line.split(" ")[:10]
    return cls(sl, s, spos, strand, ql, q, qpos, score, edist, ident)
    
  def compare_seeds(self, other, skip=[]):
    if (self.slen == other.slen and 
        self.s == other.s and
        self.spos == other.spos and
        self.strand == other.strand and
        self.qlen == other.qlen and
        self.q == other.q and
        self.qpos == other.qpos): #and
       #self.score == other.score and
       #self.edist == other.edist and
       #self.ident == other.ident:
         return True
    return False
    
  def to_line(self):
    return "%s %s %s %s %s %s %s %s %s %s" %(self.slen, self.s, self.spos, 
                                             self.strand, self.qlen, self.q,
                                             self.qpos, self.score, self.edist,
                                             self.ident)

class Seed:
  def __init__(self, fail, strand, s, spos, q, qpos, l):
    self.fail = fail
    self.strand = strand
    self.s = s
    self.spos = spos
    self.q = q
    self.qpos = qpos
    self.l = l

  def to_line(self):
    return "%s,%s,%s,%s,%s,%s,%s\n" %("-" if self.fail else "+", 
                                      self.strand, self.s, self.spos, 
                                      self.q, self.qpos, self.l)


def main():
  opts = parse_opts()
  params = check_opts(opts)
  
  encode(params["infile"])
  seeds = do_gt_extend(params)
  succ, fail, gt_extend = parse_seeds(seeds)
  
  seeds_to_file("seeds.txt", succ)
  sq_out = run_with_seqan(params["seqan"], params["infile"], "seeds.txt")
  
  sq_extend = []
  for s in sq_out:
    sq_extend.append(ExtendedSeed.from_string(s))
  
  print("got %s seeds from gt (%s success, %s fail)" %(len(succ + fail), len(succ), len(fail)))
    
  matchcount = 0
  print("matching extensions:")
  for gt in gt_extend:
    for sq in sq_extend:
      if gt.compare_seeds(sq, ["score", "edist", "ident"]):
        matchcount += 1
        print("\tgt: %s\n\tsq: %s\n" %(gt.to_line(), sq.to_line()))
  if matchcount == 0:
    print("none")
  else:
    print("%s seeds match (%.2f percent of total)" %(matchcount, 
                                        (float(matchcount)/len(succ))*100))
  
  return 0
  
  

if __name__ == '__main__':
    main()














