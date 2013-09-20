from entities import *
from mesh     import *

class triangle:
    def __init__(self, p0, p1, p2):
        #self.id     = 0
        self.points = [p0, p1, p2]
        #self.edges  = [[p0, p1], [p1,p2],[p2,p0]]
        self.adjacents  = [None, None, None]

    def set_adjacents(self, a0, a1, a2):
        self.adjacents  = [a0, a1, a2]

def locate_cell(cell_ini, x, y):
    #print "finding cell, cell_ini:",
    #print_cell(cell_ini)
    for i in range(3):
        pa = cell_ini.points[i]
        pb = cell_ini.points[(i+1)%3]
        #print "right:",i, (pb.x-pa.x)*(y-pa.y) - (pb.y-pa.y)*(x-pa.x)
        if (pb.x-pa.x)*(y-pa.y) - (pb.y-pa.y)*(x-pa.x) < 0: # point lies on right side
            return locate_cell(cell_ini.adjacents[i], x, y)

    return cell_ini


def print_cell(cell, spc=""):
    c = cell
    if c:
        print spc, "id:", c.id, (c.points[0].id, c.points[1].id, c.points[2].id), "coords:",(c.points[0].x, c.points[0].y), (c.points[1].x, c.points[1].y), (c.points[2].x, c.points[2].y)
    else:
        print spc, c

def print_cells(cells):
    for c in cells:
        print_cell(c)
        print "    neigh:"
        print_cell(c.adjacents[0], "     ")
        print_cell(c.adjacents[1], "     ")
        print_cell(c.adjacents[2], "     ")
        print

def check_adjacents(cell):
    #return
    ok = True
    for i in range(3):
        adjac = cell.adjacents[i]
        if adjac==None: continue

        pa = cell.points[i]
        pb = cell.points[(i+1)%3]
        if pa not in adjac.points: ok = False
        if pb not in adjac.points: ok = False

    if not ok:
        raise Exception("Adjacent fail for cell")

def swap_cells(cell0, cell1, p0_idx=None, check_adjacency=False):

    # Search p0
    if not p0_idx:
        for idx_p0 in range(3):
            p0 = cell0.points[idx_p0]
            if p0 not in cell1.points:
                break
    else:
        p0 = cell0.points[idx_p0]

    # Get v0 and v1
    v0 = cell0.points[idx_p0-1]
    v1 = cell0.points[idx_p0-2]

    # Get p1
    idx_p1 = cell1.points.index(v1)-2
    p1     = cell1.points[idx_p1]


    # Calculate if p0 from cell 0 is in circuncircle of cell1
    x01 = v0.x - p1.x  ;  y01 = v0.y - p1.y
    x11 = v1.x - p1.x  ;  y11 = v1.y - p1.y
    x00 = v0.x - p0.x  ;  y00 = v0.y - p0.y
    x10 = v1.x - p0.x  ;  y10 = v1.y - p0.y

    if (x01*x11 + y01*y11)*(x10*y00 - x00*y10) < (y01*x11 - x01*y11)*(x10*x00 + y00*y10):

        # related adjacent cells
        a = cell1.adjacents[idx_p1-1]
        b = cell1.adjacents[idx_p1]
        c = cell0.adjacents[idx_p0-1]
        d = cell0.adjacents[idx_p0]

        cell0.points = [p0, v1, p1]
        cell1.points = [p0, p1, v0]
        cell0.adjacents = [d, a, cell1]
        cell1.adjacents = [cell0, b, c]

        if a: a.adjacents[ a.points.index(p1) ] = cell0
        if c: c.adjacents[ c.points.index(p0) ] = cell1

        if check_adjacency:
            check_adjacents(cell0)
            check_adjacents(cell1)
            if a: check_adjacents(a)
            if b: check_adjacents(c)

        return True

    return False


def triangulate_polygon(pts_list, points, cells, check_adjacency=False):
    local_cells  = CollectionCell()

    minx = min([p[0] for p in pts_list])
    miny = min([p[1] for p in pts_list])
    maxx = max([p[0] for p in pts_list])
    maxy = max([p[1] for p in pts_list])

    dx = maxx - minx
    dy = maxy - miny
    de = max(dx, dy, 1.0)

    # Supertriangle
    pp0 = Point(minx-de  , miny-de)
    pp1 = Point(minx+4*de, miny-de)
    pp2 = Point(minx-de  , miny+4*de)

    pp0.id = -1
    pp1.id = -2
    pp2.id = -3

    cell_ini = local_cells.add_new(TRI3, [pp0, pp1, pp2])
    cell_ini.adjacents = [None, None, None]

    for pt in pts_list:
        P = points.add_new(pt)
        #print ">>>Introducing point " , P

        # Locate container cell
        curr_cell = locate_cell(local_cells[-1], P.x, P.y)

        # Get cell info
        p0, p1, p2 = curr_cell.points
        a0, a1, a2 = curr_cell.adjacents

        # Define new cells
        curr_cell.points = [P, p0, p1]
        t0 = curr_cell
        t1 = local_cells.add_new(TRI3, [P, p1, p2])
        t2 = local_cells.add_new(TRI3, [P, p2, p0])
        t0.adjacents = [t2, a0, t1]
        t1.adjacents = [t0, a1, t2]
        t2.adjacents = [t1, a2, t0]

        if a0: a0.adjacents[ a0.points.index(p1) ] = t0
        if a1: a1.adjacents[ a1.points.index(p2) ] = t1
        if a2: a2.adjacents[ a2.points.index(p0) ] = t2

        if check_adjacency:
            check_adjacents(t0)
            check_adjacents(t1)
            check_adjacents(t2)
            if a0: check_adjacents(a0)
            if a1: check_adjacents(a1)
            if a2: check_adjacents(a2)

        # Check circumcircle (swap)
        cells_for_swap = [t0, t1, t2]

        while cells_for_swap:
            cell0 = cells_for_swap.pop()
            cell1 = cell0.adjacents[1] # adjacent cell opposito to point P in cell0 cell

            if not cell1:
                continue

            if swap_cells(cell0, cell1, 0):
                # add cells to check
                cells_for_swap.append(cell0)
                cells_for_swap.append(cell1)

    #points.append(pp0)
    #points.append(pp1)
    #points.append(pp2)

    # Remove extra related to supertriangle

    # Find a triangle that contains one supertriangle point
    curr_cell = None
    for cell in local_cells:
        if pp0 in cell.points:
            curr_cell = cell
            break

    # Iteratively remove neighbors to curr_cell
    spoints =  [pp0, pp1, pp2]

    while curr_cell:
        # remove adjacency from neighbors
        adjs = curr_cell.adjacents
        pts  = curr_cell.points
        for p, adj in zip(pts, adjs):
            if adj:
                adj.adjacents[adj.points.index(p)-1] = None

        # remove curr_cell
        local_cells[curr_cell.id] = None

        # find next cell
        for c in adjs:
            if c:
                p0, p1, p2 = c.points
                if any([p0 in spoints, p1 in spoints, p2 in spoints]):
                    curr_cell = c
                    break
        else:
            break



    #while True:
    #    # remove adjacency from neighbors
    #    adjs = curr_cell.adjacents
    #    for p, adj in zip(curr_cell.points, adjs):
    #        if adj is None: continue
    #        adj.adjacents[adj.points.index(p)-1] = None

    #    # remove curr_cell
    #    local_cells[curr_cell.id] = None
    #    for c in adjs:

    #        if not c:
    #            continue

    #        p0, p1, p2 = c.points
    #        if any([p0 in spoints, p1 in spoints, p2 in spoints]):
    #            curr_cell = c
    #            break
    #    else:
    #        break

    local_cells = [c for c in local_cells if c is not None]


    cells.extend(local_cells)
    local_cells = None


if __name__=="__main__":

    points_lst = []
    points_lst.append([0,0])
    points_lst.append([1,0])
    points_lst.append([1,1])
    points_lst.append([0,1])


    from random import *
    for k in range(1000):
        points_lst.append( [random(), random()] )


    #points.sort(key=lambda p:p.__hash__())
    #for p in points[:10]:
        #print p

    points = CollectionPoint()
    cells  = CollectionCell()

    triangulate_polygon(points_lst, points, cells)


    Ps = CollectionPoint()
    for p in points:
        Ps.add_new([p.x, p.y])


    Cs = CollectionCell()
    for c in cells:
        Cs.add_new(TRI3, c.points)



    mesh = Mesh()
    mesh.points = Ps
    mesh.cells  = Cs

    mesh.write_file("delaunay0.vtk")

