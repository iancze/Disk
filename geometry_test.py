import geometry as g
import numpy as np
import matplotlib.pyplot as plt


def test_delta_0():
    coord0 = g.StartingCoords(0.,0.0,0.0)
    coord1 = g.StartingCoords(1.0,0.0,0.0)
    coord2 = g.StartingCoords(1.11,0.0,0.0)
    coord3 = g.StartingCoords(1.5,0.0,0.0)
    coord4 = g.StartingCoords(np.pi/2,0.0,0.0)
    print(coord0)
    print(coord1)
    print(coord2)
    print(coord3)
    print(coord4)


def test_delta():
    coord0 = g.StartingCoords(0.,1.0,0.0)
    coord1 = g.StartingCoords(1.0,1.0,0.0)
    coord2 = g.StartingCoords(1.11,1.0,0.0)
    coord3 = g.StartingCoords(1.5,1.0,0.0)
    coord4 = g.StartingCoords(np.pi/2,1.0,0.0)
    print(coord0)
    print(coord1)
    print(coord2)
    print(coord3)
    print(coord4)

def test_plot():
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111)
    coord1 = g.StartingCoords(1.0,1.0,0.0)
    print(coord1)
    coord1.graph_on_ax(ax)
    fig.savefig("test_plot.png")


def main():
    #test_delta_0()
    #test_delta()
    test_plot()
    pass

if __name__=="__main__":
    main()
