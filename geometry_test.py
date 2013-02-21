import geometry as g
import numpy as np


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


def main():
    test_delta_0()

if __name__=="__main__":
    main()
