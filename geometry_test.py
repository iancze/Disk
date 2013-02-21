import geometry as g
import numpy as np
import matplotlib.pyplot as plt


def test_delta_0():
    orientation0 = g.Orientation(0.0)
    orientation1 = g.Orientation(1.0)
    orientation2 = g.Orientation(1.11)
    orientation3 = g.Orientation(1.5)
    orientation4 = g.Orientation(np.pi/2)
    los0 = g.LineOfSight(orientation0,0.0,0.0)
    los1 = g.LineOfSight(orientation1,0.0,0.0)
    los2 = g.LineOfSight(orientation2,0.0,0.0)
    los3 = g.LineOfSight(orientation3,0.0,0.0)
    los4 = g.LineOfSight(orientation4,0.0,0.0)

    print(orientation0)
    print(los0)

    print(orientation1)
    print(los1)

    print(orientation2)
    print(los2)

    print(orientation3)
    print(los3)

    print(orientation4)
    print(los4)


def test_plot():
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111)
    orn1 = g.Orientation(1.0)
    los1 = g.LineOfSight(orn1,1.0,0.0)
    print(orn1)
    print(los1)
    los1.graph_on_ax(ax)
    fig.savefig("test_plot.png")


def main():
    #test_delta_0()
    #test_delta()
    test_plot()
    pass

if __name__=="__main__":
    main()
