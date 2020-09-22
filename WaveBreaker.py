from mshr import Polygon, Rectangle
from fenics import Point
import numpy as np

def build_fin(loc, l, d):
    """
    Build a fin as a domain to be meshed

    Input: 
        loc - location of outer most corner coordinates [x, y]
        l - length of fin
        d - thickness of fin
    Output:
        fin - fin domain
    """

    # define the fin as a polygon
    angle = np.radians(45)
    # vertices, proceeding counter-clockwise
    pts = np.zeros((2, 6))
    # leave pts[:, 0] alone
    pts[:, 1] = [l * np.cos(angle), l * np.sin(angle)]
    pts[:, 2] = [(l + d) * np.cos(angle), (l - d) * np.sin(angle)]
    pts[:, 3] = [d / np.cos(angle), 0]
    
    pts[0, 4] = pts[0, 2]  # flip point along x-axis
    pts[1, 4] = -pts[1, 2]
    pts[0, 5] = pts[0, 1]
    pts[1, 5] = -pts[1, 1]

    # shift over to loc
    # convert to mshr Points
    points = []
    for j in range(6):
        points.append(Point(pts[:, -j] + loc))

    # build polygon
    poly = Polygon(points)
    return poly


def fin_sequence(start_loc, length, n_fins, l, d):
    """
    Build a sequence of fins

    Input:
        start_loc - starting location for placing fins
        length - length of array of fins, eg. 2in
        n_fins - number of fins to place in that length
        l - length of fin wing
        d - width of fin wing
    Output:
        fins - array of fins
    """

    # starting locations
    locs = np.zeros((2, n_fins))
    locs[1, :] = start_loc[1]
    locs[0, :] = np.linspace(start_loc[0], start_loc[0] + length, n_fins)

    fins = []
    for j in range(n_fins):
        fins.append(build_fin(locs[:, j], l, d))

    return fins

def build_wave_breaker(start_loc, length, n_fins, l, d, rod_start, rod_stop):
    """
    Build the wave breaker

    Input:
        start_loc - starting location
        length - length of breaker
        n_fins - number of fins
        l - length of fin wing
        d - width of fin wing
        rod_start - starting position for support rod
        rot_stop - stopping position for support rod
    Output:
        breaker - domain of wave breaker
    """

    # build support rod
    rod_pt1 = Point(rod_start)
    rod_pt2 = Point(rod_stop)

    # get fins
    fins = fin_sequence(start_loc, length, n_fins, l, d)

    # assemble together
    domain = Rectangle(rod_pt1, rod_pt2)

    for fin in fins:
        # += may have odd behavior
        domain = domain + fin

    return domain


