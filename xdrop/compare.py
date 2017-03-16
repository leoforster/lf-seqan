import sys
import os
import subprocess
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input fasta file", type=str)
parser.add_argument("-s", help="seedlength", type=int)
args = parser.parse_args()

if args.s:
  seedlen = args.s
else:
  seedlen = 8

if args.i != None:
  if args.i[0] == "/":
    infile = args.i
  else:
    infile = os.getcwd() + "/" + args.i
    
  if not os.path.isfile(infile): #and not os.path.isfile(infile + ".esq"):
    print "fail %s" %infile
    sys.exit()
else:
  print "fail %s" %infile
  sys.exit()

#encseq encode places encoded files in script calling directory -- change it?
try:
  subprocess.call(["gt encseq encode %s" %infile], shell=True, stderr=subprocess.STDOUT)
except:
  print "failed to encode"
  sys.exit()

#gt seed_extend -ii longseeded1 -outfmt seed failed_seed -maxmat 2 -l 100 -extendxdrop -seedlength 20
try:
  seed_out = subprocess.check_output(["gt seed_extend -ii %s -seedlength %i -outfmt seed failed_seed" %(infile, seedlen)], shell=True, stderr=subprocess.STDOUT)
  seed_out = seed_out.split("\n")
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
for seed in seed_out:
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
#outfile = "seeds.txt"
#with open(outfile, "w") as f:
  #for seed in seedlist:
    #f.write(seed[0][0] + ", " + seed[0][1] + ", " + 
            #seed[1][0] + ", " + seed[1][1] + ", " + 
            #seed[2] + "\n")
  
fasta_sequences = SeqIO.parse(open(infile),'fasta')
  
seqan_path = "~/Desktop/projects/hiwi/lf-seqan/xdrop/seqan/extend.x"
seqlist = []; count = 0
for fasta in fasta_sequences:
  count += 1
  seqlist.append(str(fasta.seq))

seqan_path += " %i" %count
seqan_path += " %i" %len(seedlist)

for seq in seqlist:
  seqan_path += " " + seq
  
for seed in seedlist:
  seqan_path += " " + seed[0][0] + " " + seed[1][0] + " " + seed[0][1] + " " + seed[1][1] + " " + seed[2]

try:
  extend_out = subprocess.check_output([seqan_path], shell=True, stderr=subprocess.STDOUT)
  extend_out = extend_out.split("\n")
except:
  print "failed to extend seeds"
  sys.exit()

def parts(match):
  try:
    [slen, sid, spos, strand, qlen, qid, qpos, score, edist, ident] = match.split(" ")[:10]
    return [slen, sid, spos, strand, qlen, qid, qpos, score, edist, ident]
  except:
    return None
  
def compare_match(seqan, gt):
  p1 = parts(seqan)
  p2 = parts(gt)
  
  if p1 == None or p2 == None:
    return False
  
  for i in range(len(p1)):
    if i in [3, 7, 8, 9]: #skip unsure lines
      continue
    if p1[i] != p2[i]:
      return False
  return True

#integrate above?
cand = []
for seed in seed_out:
  if seed != "":
    if not seed.startswith("#"):
      cand.append(seed)

matches = []
for s in extend_out:
  for c in cand:
    if compare_match(s, c):
      matches.append(s)
      
print matches














