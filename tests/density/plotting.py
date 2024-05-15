import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

inputs = [];
for i in range(12):
    temp = [];
    inputs.append(temp);

fp = open("out1.txt");
for i in range(12):
    line = fp.readline();
    inputs[i] = line.rstrip().split();

data = [];
for i in range(12):
    temp = [];
    data.append(temp);
for i in range(12):
    for j in range(100):
        data[i].append(float(inputs[i][j]));

"""
for i in range(12):
    for j in range(100):
        
        print(data[i][j], end = " ")
    print("");
"""
expected = [0.009174, 0.007812, 0.0009911, 0.0009766, 0.25, 0.2222, 0.1111, 0.1176, 0.2857, 0.25, 0.1176, 0.125]
titles = ["ModMin: Mod is 109", "ModMin: Mod is 128", "ModMin: Mod is 1009", "ModMin: Mod is 1024", "WindowMin: Large Window is 7", "WindowMin: Large Window is 8", "WindowMin: Large Window is 17", "WindowMin: Large Window is 16", "Syncmer: Large Window is 7", "Syncmer: Large Window is 8", "Syncmer: Large Window is 17", "Syncmer: Large Window is 16"]
for i in range(12):
    fig, axs = plt.subplots(1, 1, figsize =(10, 7), tight_layout = True)
    
    axs.xaxis.set_ticks_position('none')
    axs.yaxis.set_ticks_position('none')
    
    axs.xaxis.set_tick_params(pad = 5)
    axs.yaxis.set_tick_params(pad = 10)
    
    N, bins, patches = axs.hist(data[i], bins = 10)
    plt.axvline(expected[i], color='k', linestyle='dashed', linewidth=1)
    fracs = ((N**(1 / 5)) / N.max())
    norm = colors.Normalize(fracs.min(), fracs.max())
     
    for thisfrac, thispatch in zip(fracs, patches):
        color = plt.cm.viridis(norm(thisfrac))
        thispatch.set_facecolor(color)
    plt.xlabel("X-bar")
    plt.ylabel("Count")
    legend_str = " ,E[X] = " + str(expected[i])
    plt.title(titles[i] + legend_str)
    plt.show()

