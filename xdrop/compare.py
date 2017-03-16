import sys
import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input fasta file", type=str)
args = parser.parse_args()

if args.i != None:
  if args.i[0] == "/":
    infile = args.i
  else:
    infile = os.getcwd() + "/" + args.i
    
  if not os.path.isfile(infile): #and not os.path.isfile(infile + ".esq"):
    print "fail %s" %infile
    sys.exit()

try:
  subprocess.call(["gt encseq encode %s" %infile], shell=True, stderr=subprocess.STDOUT)
except:
  print "failed to encode"
  sys.exit()

#gt seed_extend -ii longseeded1 -outfmt seed failed_seed -maxmat 2 -l 100 -extendxdrop -seedlength 20
try:
  out = subprocess.check_output(["gt seed_extend -ii %s -seedlength 8 -outfmt seed failed_seed" %infile], shell=True, stderr=subprocess.STDOUT)
except:
  print "failed to get seeds"
  sys.exit()

def makeseed(sid, qid, seed):
  #s_seed = (int(sid), int(seed[1]))
  #q_seed = (int(qid), int(seed[2]))
  #seedlen = int(seed[0])
  s_seed = (sid, seed[1])
  q_seed = (qid, seed[2])
  seedlen = seed[0]
  return (s_seed, q_seed, seedlen)

#s.len, s.seqnum, s.start, strand, 
#q.len, q.seqnum, q.start, score, editdist, identity, 
#seedlen, s.seedstart, q.seedstart
seedlist = []
for seed in out.split("\n"):
  if seed != "":
    if seed.startswith("#"):
      if seed.startswith("# failed_seed"):
        seed = seed.split("\t")
        seedlist.append(makeseed(seed[2], seed[5], [seed[1], seed[3], seed[6]]))
      else:
        pass
    else:
      #todo: parse rest for match comparison
      #todo: forward/reverse strand
      seed = seed.split(" ")
      seedlist.append(makeseed(seed[1], seed[5], seed[10:]))

#print len(seedlist)
outfile = "seeds.txt"
with open(outfile, "w") as f:
  for seed in seedlist:
    f.write(seed[0][0] + ", " + seed[0][1] + ", " + 
            seed[1][0] + ", " + seed[1][1] + ", " + 
            seed[2] + "\n")

try:
  out = subprocess.check_output(["~/Desktop/projects/hiwi/lf-seqan/xdrop/seqan/test.x"], shell=True, stderr=subprocess.STDOUT)
except:
  print "failed to extend seeds"
  sys.exit()
