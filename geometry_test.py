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


def test_line_of_sight_w_theta():
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111)
    deg = 89.999
    orn1 = g.Orientation(deg * np.pi/180.)
    delta_array = np.linspace(-1*orn1.delta_limit, orn1.delta_limit,num=20)
    los_array = [g.LineOfSight(orn1, i, 0.0) for i in delta_array]
    orn1.graph_on_ax(ax)
    print(orn1)
    for i in los_array:
        i.graph_on_ax(ax)
        #print(i)
    fig.savefig("Plots/Test_Theta/test_plot{theta:.0f}.png".format(theta=(orn1.theta*180./np.pi)))

def test_path_length(deg):
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111)
    #deg = 50.0
    orn1 = g.Orientation(deg * np.pi/180.)
    orn1.graph_on_ax(ax)
    delta_array = np.linspace(-1*orn1.delta_limit, orn1.delta_limit,num=20)
    los_array = [g.LineOfSight(orn1, i, 0.0) for i in delta_array]
    for i in los_array:
        i.graph_on_ax(ax)
        i.plot_path(ax)
        print(i)
        i.check_total_path_length()
    fig.savefig("Plots/Test_Path_Length/path_length{theta:.0f}.png".format(theta=(orn1.theta*180./np.pi)))

def test_polar(deg):
    orn1 = g.Orientation(deg * np.pi/180.)
    los1 = g.LineOfSight(orn1, 1.0, 0.0)
    s_test = los1.s_finish*0.6
    print(s_test)
    print(los1.cartesian(s_test))
    print(los1.polar(s_test))
    return

def test_coords(deg):
    orn1 = g.Orientation(deg * np.pi/180.)
    los1 = g.LineOfSight(orn1, 1.0, 0.0)
    s_test = los1.s_finish*0.6
    print(s_test)
    print(los1.dI(s_test))
    return


def main():
    #test_delta_0()
    #test_delta()
    #test_plot()
    #for deg in np.linspace(0.01,89.99,num=10):
    #    test_path_length(deg)
    #test_path_length(50.)
    #test_polar(50.)
    test_coords(50.)
    pass

if __name__=="__main__":
    main()
