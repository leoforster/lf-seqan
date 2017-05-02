import numpy as np
import matplotlib.pyplot as plt
import sys
import os

files = sys.argv[1:-1]
times = {}
for i in files:
    times.setdefault(i, [])

avg = []
std = []
normed = {}
for p in files:
    with open(p, "r") as f:
        for i in f.readlines():
            if i.startswith("real"):
                times[p].append(float(i.split("m")[1].strip("s\n")))
                
    avg.append(np.average(times[p]))
    std.append(np.std(times[p]))
    
    normed[p] = [(x - min(times[p])) / (max(times[p]) - min(times[p])) for x in times[p]]

fig, (ax1, ax2) = plt.subplots(2)
axis = np.arange(len(times[list(times)[0]]))

for p in normed:
    ax1.plot(axis, normed[p])
ax1.set_title("normalized times for {x} measurements".format(x=len(times[list(times)[0]])))

ticks = np.arange(len(list(times)))
ax2.bar(ticks, avg, yerr=std, alpha=0.5, capsize=3, align="center")
ax2.set_xticks(ticks)
ax2.set_xticklabels(files)
ax2.set_title("aligning " + sys.argv[-1])

plt.tight_layout()
plt.show()