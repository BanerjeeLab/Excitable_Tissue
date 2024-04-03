import numpy as np
import matplotlib.pyplot as plt
import sys

def save(name):
    plt.savefig(f"{name}.pdf", bbox_inches="tight")

def getJunctionData(v1, v2, nSteps):
    length = []
    length0 = []
    gammaA = []
    timeActive = []

    for i in range(nSteps):
        data = np.loadtxt(f"./lengths_tensions/{i}.txt")
        v1s = data[:,0].astype(int)
        v2s = data[:,1].astype(int)
        j = np.where((v1s == v1) & (v2s == v2))[0][0]

        length.append(data[:,2][j])
        length0.append(data[:,3][j])
        gammaA.append(data[:,4][j])
        timeActive.append(data[:,5][j])

    length = np.array(length)
    length0 = np.array(length0)
    strain = length / length0 - 1
    gammaA = np.array(gammaA)
    timeActive = np.array(timeActive)
    return length, length0, strain, gammaA, timeActive

def getActivationTimes(nSaves, timeActive, tau, time):
    activationTimes = []
    for i in range(nSaves):
        tAct = timeActive[i]
        t = time[i]
        if tAct != -1 and (len(activationTimes) == 0 or abs(t - tAct - activationTimes[-1]) > tau / 100):
            activationTimes.append(t - tAct)
    return activationTimes

def plotJunctionData(junctions, totalTime, nSaves, tau, name="junction_data"):
    time = np.linspace(0, totalTime, nSaves, endpoint=False)
    fig, (ax1, ax0, ax2, ax3) = plt.subplots(4, 1, figsize=(5,5.5), sharex=True)

    # Plot timeseries data.
    for j in junctions:
        v1, v2 = sorted(j)
        length, length0, strain, gammaA, timeActive = getJunctionData(v1, v2, nSaves)

        ax1.scatter(time, length, marker=".",s=12, color="k",zorder=10)
        ax1.plot(time, length, linewidth=1,color="k",zorder=10)
        ax2.scatter(time, strain, marker=".",s=12,color="k",zorder=10)
        ax2.plot(time, strain, linewidth=1,color="k",zorder=10)
        ax3.scatter(time, gammaA, marker=".",s=12,color="k",zorder=10)

        ax0.scatter(time, length0, marker=".",s=12, color="k",zorder=10)
        ax0.plot(time, length0, linewidth=1,color="k",zorder=10)

        ax2.axhline(y=0.1, linestyle="--", color="red",label="$\\varepsilon_{\mathrm{on}}$")
        ax2.legend()

        # Get activation times.
        activationTimes = getActivationTimes(nSaves, timeActive, tau, time)

        # Plot junction state.
        for t in activationTimes:
            ax3.fill_between([0,t], [0.6,0.6],[-0.5,-0.5],zorder=1,color="#C1CDCD",alpha=0.3)
            ax3.fill_between([t+2*tau,8], [0.6,0.6],[-0.5,-0.5],zorder=1,color="#C1CDCD",alpha=0.3)
            if t + tau <= totalTime:
                #ax3.axvspan(t, t+tau, alpha=0.3, color='red')
                ax3.fill_between([t, t+tau], [0.6,0.6],[-0.5,-0.5],zorder=1,color="maroon",alpha=0.3)
            if t + 2 * tau <= totalTime:
                #ax3.axvspan(t+tau, t+2*tau, alpha=0.3, color='blue')
                ax3.fill_between([t+tau, t+2*tau], [0.6,0.6],[-0.5,-0.5],zorder=1,color="dodgerblue",alpha=0.3)

    ax2.set_ylim((-0.18, 0.18))
    ax3.set_ylim((-0.1, 0.55))
    ax0.set_ylim((0.5,0.7))
    ax1.set_ylim((0.5,0.7))

    ax1.set_xlim((0, 8))
    ax2.set_xlim((0, 8))
    ax3.set_xlim((0, 8))

    ax1.set_ylabel("Junction \n Length ($\sqrt{A_0}$)")
    ax0.set_ylabel("Junction Rest \n Length ($\sqrt{A_0}$)")
    ax2.set_ylabel("Junction \n Strain")
    ax3.set_ylabel("Active Contractility \n ($KA_0$)")

    ax3.set_xlabel("Time (min)")

    save(name)

# Setup

plt.rcParams['text.usetex'] = True
dt = 0.01

stepsPerSave = 25 // 5
nSaves = 300
totalTime = dt * stepsPerSave * nSaves
tau = 2
junction = [(79, 229)]
plt.rcParams.update({'font.size': 12})
plotJunctionData(junction, totalTime, nSaves, tau,name="single_junction_data")
