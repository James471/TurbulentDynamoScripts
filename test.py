import matplotlib.pyplot as pl
import designParams

fig1 = pl.figure()
fig2 = pl.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)

ax1.plot([1, 2, 3], [1, 2, 3])
ax2.plot([1, 2, 3], [3, 2, 1])

ax1.set_xlabel("Test")
ax2.set_ylabel("Test")

fig1.savefig("fig1.png")
fig2.savefig("fig2.png")
