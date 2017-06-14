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
            if i.startswith("user"):
                t = 60 * float(i.split("m")[0].split("\t")[1]) + float(i.split("m")[1].strip("s\n"))
                times[p].append(t)

    avg.append(np.average(times[p]))
    std.append(np.std(times[p]))

    if len(times[p]) == 1 or min(times[p]) == max(times[p]):
        normed[p] = times[p]
    else:
        normed[p] = [(x - min(times[p])) / (max(times[p]) - min(times[p])) for x in times[p]]

fig, (ax1, ax2) = plt.subplots(2)
axis = np.arange(len(times[list(times)[0]]))

for p in normed:
    ax1.plot(axis, normed[p], label=p)
ax1.set_title("normalized times for {x} measurements".format(x=len(times[list(times)[0]])))
ax1.legend(loc=1)

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax2.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%.3f' % float(height),
                ha='center', va='bottom', alpha=0.7)

ticks = np.arange(len(list(times)))
rects = ax2.bar(ticks, avg, yerr=std, alpha=0.5, capsize=3, align="center")
ax2.set_xticks(ticks)
ax2.set_xticklabels(files)
ax2.set_title("aligning " + sys.argv[-1])
autolabel(rects)

plt.tight_layout()
plt.show()
