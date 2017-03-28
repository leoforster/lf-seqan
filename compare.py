#!/usr/bin/env python

import sys
import os
import subprocess
import argparse

verbose = False

#todo: mincov in seqan
#todo: -maxmat and -l -> maximal matches of minimum length l
#todo: chaining
#todo: input an encoded file

#todo: option to compare successful seeds, failed seeds, or both
def parse_opts():
  parser = argparse.ArgumentParser(
    description="compare gt seed_extend seeds with those from the seqan library.")
  parser.add_argument(
    "-s", "--script", help="compiled seqan script location (required)",
    type=str, required=True)
  parser.add_argument(
    "-i", "--input", help="input fasta file (required)", type=str, required=True)
  parser.add_argument(
    "-q", "--query", help="input query fasta file for multiple comparison", type=str)
  parser.add_argument(
    "-o", "--output", help="base name for saving output of gt and seqan seed extensions", type=str)
  parser.add_argument(
    "-l", "--seedlen", help="seedlength", type=int)
  parser.add_argument(
    "-p", "--seedfile", help="path of file to write seed information", type=str)
  parser.add_argument(
    "-n", "--mincov", help="minimum coverage for seed_extend", type=int)
  parser.add_argument(
    "-y", "--noxdrop", help="dont use xdrop algorithm for seed extension", action="store_true")
  parser.add_argument(
    "-g", "--greedy", help="use greedy algorithm for seed extension", action="store_true")
  parser.add_argument(
    "-c", "--xcut", help="xdrop cutoff score", type=int)
  parser.add_argument(
    "-m", "--minid", help="minimum identity of matches for seed_extend", type=int)
  parser.add_argument(
    "--compare", help="use internal match comparison", action="store_true")
  parser.add_argument(
    "-v", "--verbose", help="generate additional output (verbose mode)", action="store_true")
  args = parser.parse_args()
  return args
  
def check_opts(args):
  params = {
    "seqan":None,
    "infile":None,
    "qfile":None,
    "outfile":None,
    "seedlen":None,
    "mincov":None,
    "minid":None,
    "do_xdrop":True, #default!
    "do_greedy":False,
    "xcutoff":None, 
    "compare":False,
    "seedfile":None
    }

  #seqan
  if args.script[0] == "/":
    params["seqan"] = args.script
  else:
    params["seqan"] = os.getcwd() + "/" + args.script
  #check if executable?
  if not os.path.isfile(args.script):
    print("failed to find %s - does it exist?" %args.script)
    sys.exit()
  
  #infile
  if args.input[0] == "/":
    params["infile"] = args.input
  else:
    params["infile"] = os.getcwd() + "/" + args.input
  if not os.path.isfile(args.input):
    print("failed to find %s - does it exist?" %args.input)
    sys.exit()
      
  #outfile
  if args.output:
    outpath = os.path.split(args.output)[0]
    if os.path.exists(outpath) or outpath == "":
      params["outfile"] = args.output
  #else:
    #params["outfile"] = "output_"
      
  #query file
  if args.query:
    if args.query[0] == "/":
      params["qfile"] = args.query
    else:
      params["qfile"] = os.getcwd() + "/" + args.query
    if not os.path.isfile(args.query):
      print("failed to find %s - does it exist?" %args.query)
      sys.exit()
      
  #seedfile
  if args.seedfile:
    seedpath = os.path.split(args.seedfile)[0]
    if os.path.exists(seedpath) or seedpath == "":
      params["seedfile"] = args.seedfile
  else:
    params["seedfile"] = "seeds.txt"

  params["seedlen"] = args.seedlen if args.seedlen else 8
  params["do_xdrop"] = False if args.noxdrop else True #default: do xdrop !
  params["do_greedy"] = True if args.greedy else False
  params["compare"] = True if args.compare else False
  if args.xcut:
    params["xcutoff"] = args.xcut
  if args.minid:
    params["minid"] = args.minid
  if args.mincov:
    params["mincov"] = args.mincov

  #xdrop, greedy exclusion
  #if args.xdrop and args.greedy:
    #print("cannot have both -x and -g")
    #sys.exit()

  #verbose
  if args.verbose:
    global verbose 
    verbose = True

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
  if params["do_xdrop"]: pstr += " -extendxdrop" 
  if params["do_greedy"]: pstr += " -extendgreedy" 
  if params["qfile"] != None: pstr += " -qii %s" %(params["qfile"]) 
  if params["seedlen"]: pstr += " -seedlength %s" %(params["seedlen"])
  if params["xcutoff"]: pstr += " -xdropbelow %s" %(params["xcutoff"])
  if params["mincov"]: pstr += " -mincoverage %s" %(params["mincov"])
  if params["minid"]: pstr += " -minidentity %s" %(params["minid"])
  
  #with -outfmt: failed_seed alignment
  call = "gt seed_extend -ii %s %s -outfmt seed" %(params["infile"], pstr)
  
  try:
    seeds = subprocess.check_output([call], shell=True, stderr=subprocess.STDOUT)
  except:
    print("failed to get seeds -- call was:\n " + call)
    sys.exit()
  if verbose: print(call)
  
  seeds = seeds.decode("utf-8").split("\n")
  return seeds
  
def run_with_seqan(seqan, seedfile, infile, qfile = None):
  if qfile:
    call = "%s 2 %s %s %s" %(seqan, seedfile, infile, qfile)
  else:
    call = "%s 1 %s %s" %(seqan, seedfile, infile)
  
  try:
    extend = subprocess.check_output([call], shell=True, stderr=subprocess.STDOUT)
  except:
    print("failed to run seqan script -- call was:\n " + call)
    sys.exit()
  if verbose: print(call)
  
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

def compare_matches(gt_matches, sq_matches, total, thresh=0.8):
  assert(total > 0)
  matchcount = 0
  
  print("matching extensions:")
  for gt in gt_matches:
    for sq in sq_matches:
      if gt.compatible(sq):
        if gt.equivalent(sq):
          matchcount += 1
          if verbose: print("\tgt: %s\n\tsq: %s\n" %(gt.to_line(), sq.to_line()))
        else:
          olap = ExtendedSeed.overlap_score(gt, sq)
          if olap > thresh:
            matchcount += 1
            if verbose: print("\tgt: %s\n\tsq: %s\n\t%.2f%% overlapping\n" 
                              %(gt.to_line(), sq.to_line(), olap * 100))
  if matchcount == 0:
    print("\tnone")
  else:
    print("\t%s seeds match (%.2f percent of total)" %(matchcount, 
                                        (float(matchcount)/total)*100))

class ExtendedSeed:
  forward_seeds = []
  reverse_seeds = []
  
  def __init__(self, sl, s, spos, strand, ql, q, qpos, score, edist, ident):
    #check spos <= slen && qpos <= qlen ?
    self.slen = int(sl)
    self.s = int(s)
    self.spos = int(spos)
    self.send = self.spos + self.slen
    self.qlen = int(ql)
    self.q = int(q)
    self.qpos = int(qpos)
    self.qend = self.qpos + self.qlen
    self.score = int(score)
    self.edist = int(edist)
    self.ident = float(ident)
    if strand == "P":
      self.rev = True
      ExtendedSeed.reverse_seeds.append(self)
    else:
      self.rev = False
      ExtendedSeed.forward_seeds.append(self)
    
  @classmethod
  def from_string(cls, line):
    sl, s, spos, strand, ql, q, qpos, score, edist, ident = line.split(" ")[:10]
    return cls(sl, s, spos, strand, ql, q, qpos, score, edist, ident)
    
  def get_forward_seeds():
    return forward_seeds
    
  def get_reverse_seeds():
    return reverse_seeds
    
  def equivalent(self, other, skip=[]):
    if (self.slen == other.slen and 
        self.s == other.s and
        self.spos == other.spos and
        self.rev == other.rev and
        self.qlen == other.qlen and
        self.q == other.q and
        self.qpos == other.qpos): #and
       #self.score == other.score and
       #self.edist == other.edist and
       #self.ident == other.ident:
         return True
    return False
  
  def compatible(self, other, lenthresh=10, posthresh=25):
    return ((self.s == other.s and self.q == other.q) and
            (abs(self.slen - other.slen) <= lenthresh and 
             abs(self.qlen - other.qlen) <= lenthresh) and
            (abs(self.spos - other.spos) <= posthresh and
             abs(self.qpos - other.qpos) <= posthresh))

  @staticmethod
  def get_overlap_length(spos, send, opos, oend):
    olap = 0
    if spos > opos:
      if send > oend:
        olap = oend - spos
      else:
        olap = send - spos
    else:
      if oend > send:
        olap = send - opos
      else:
        olap = oend - opos
    return abs(olap)
  
  def overlap_score(self, other):
    solap = ExtendedSeed.get_overlap_length(self.spos, self.send, other.spos, other.send)
    qolap = ExtendedSeed.get_overlap_length(self.qpos, self.qend, other.qpos, other.qend)
    #biased due to seedlength being counted (it always overlaps)
    return ((float(solap)/self.slen + 
             float(solap)/other.slen + 
             float(qolap)/self.qlen + 
             float(qolap)/other.qlen) / 4)

  def to_line(self):
    return "%s %s %s %s %s %s %s %s %s %s" %(self.slen, self.s, self.spos,
                                             "P" if self.rev else "F", 
                                             self.qlen, self.q, self.qpos, 
                                             self.score, self.edist,
                                             self.ident)

class Seed:
  def __init__(self, fail, strand, s, spos, q, qpos, l):
    self.fail = fail
    self.s = s
    self.spos = spos
    self.q = q
    self.qpos = qpos
    self.l = l
    self.rev = (strand == "P")
    
  def is_failed(self):
    return self.fail

  def to_line(self):
    return "%s,%s,%s,%s,%s,%s,%s\n" %("-" if self.fail else "+",
                                      "P" if self.rev else "F",
                                      self.s, self.spos, self.q, 
                                      self.qpos, self.l)


def main():  
  opts = parse_opts()
  params = check_opts(opts)
  
  encode(params["infile"])
  seeds = do_gt_extend(params)
  succ, fail, gt_extend = parse_seeds(seeds)
  
  if len(succ) == 0:
    print("no seeds were successfully extended")
    return 1
  
  seeds_to_file(params["seedfile"], succ)
  if params["qfile"]:
    sq_out = run_with_seqan(params["seqan"], params["seedfile"], 
                            params["infile"], params["qfile"])
  else:
    sq_out = run_with_seqan(params["seqan"], params["seedfile"], params["infile"])

  sq_extend = []
  for s in sq_out:
    sq_extend.append(ExtendedSeed.from_string(s))
  
  if params["compare"]:
    print("got %s seeds from gt (%s success, %s fail)" %(len(succ + fail), len(succ), len(fail)))
    compare_matches(gt_extend, sq_extend, len(succ))

  fields = "# Fields: s.len, s.seqnum, s.start, strand, q.len, q.seqnum, "\
             "q.start, score, editdist, identity" #, seedlen, s.seedstart, q.seedstart\n"
  
  if params["outfile"]:
    for suff in ["_gt", "_sq"]:
      with open(params["outfile"] + suff, "w") as f:
        f.write(fields + "\n")
        extend = sq_extend if suff == "_sq" else gt_extend
        for s in extend:
          f.write(s.to_line() + "\n")
  elif not params["compare"]:
    for i in range(2):
      print(fields)
      extend = sq_extend if i == 1 else gt_extend
      for s in extend:
        print(s.to_line())
      if i == 0:
        print("\n\n")
  #2nd fields entry before reverse seeds isnt printed here

  return 0
  
  

if __name__ == '__main__':
    main()














