#!/usr/bin/env python
import matplotlib.pyplot as plt

basedir="caso9/"
#filenames = [ "Caso10_L2_0.2.dat", "Caso10_L2_0.1.dat","Caso10_L2_0.05.dat","Caso10_L2_0.025.dat","Caso10_L2_0.0125.dat","Caso10_L2_0.dat"]
filenames = ["Caso9_L2_0.dat"]
for filename in filenames:
    h = list()
    m2err = dict()
    with open(filename, "r") as f:
        for line in f:
            if "h" in line:
                for w in line.split(','):
                    if "h" in w:
                        i=w.find("=")
                        h.append(float(w[i+1:]))
                    else:
                        h.append(float(w))
            else:
                first = True
                current = list()
                for w in line.split():
                    if "#" in w:
                        pass
                    elif first:
                        first = False
                        m=int(w[0:-1])
                    else:
                        current.append(float(w[0:-1]))
                if not first:
                    m2err[m]=current;
    counter = 0
    for hc in h:
        new_line_plot = list()
        for item in m2err:
            elem = (item,m2err[item][counter])
            new_line_plot.append(elem)
            print elem
        new_line_plot.sort(key = lambda el: el[0])
        aus =  zip(*new_line_plot)
        plt.loglog(aus[0],aus[1], label = filename)
        counter = counter + 1
    if filename is filenames[1]:
        for pippo in range(0,4):
            new_line_plot = list()
            for item in m2err:
                elem = (item,m2err[item][counter])
                new_line_plot.append(elem)
            counter =counter +1
            new_line_plot.sort(key = lambda el: el[0])
            aus =  zip(*new_line_plot)
            plt.loglog(aus[0],aus[1])
    f.close()
plt.legend()
plt.show()
    

        
