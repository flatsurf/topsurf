r"""
This file provides the operations that appear in the gasket decomposition of 3-colored planar maps
"""

from topsurf import *

def atomic_3colmap(color=0):
    return ColoredOrientedMap(vp=[], vcolors={-1:color}, mutable = True)

color_corresp = {0: 'white', 1: 'black', 2:'red'}

def plot_3col(m):
    return m.plot(vertex_colors=color_corresp)
        

def bridge_3col(m1, m2, r1=None, r2=None, check=True):
    r"""
    Add a copy of m2 in m1 and joins them with an edge from the root corner of m1 to the co-root corner of m2. If the root half edge is not provided it is assumed to be 2n-2. 

    INPUT:
        - ``m1``, ``m2`` : gray maps
        - ``r1``, ``r2`` : their root half edges, if ``None`` then assumed to be 2ni-2

    EXAMPLES::

            sage: M = atomic_3colmap()
            sage: bridge_3col(M, atomic_3colmap())
            ColoredOrientedMap("(0)(~0)", "(0,~0)", edge colors: [None], vertex colors: {0: 0, 1: 1})
            
            sage: M1 = ColoredOrientedMap(vp="(0,1)(~0)(~1)", vcolors={0: 0, 1: 1, 3: 1}, mutable = True)
            sage: bridge_3col(M1, M)
             ColoredOrientedMap("(0,1,3)(~0)(~1)(2)(~2,~3)", "(0,~0,3,~2,2,~3,1,~1)", 
             edge colors: [None, None, None, None], vertex colors: {0: 0, 1: 1, 3: 1, 4: 0, 5: 1})
    """

    if check:
        m1._assert_mutable()
        m1._check()
        m2._check()
    
    if r1 is None :
        r1 = len(m1._vp)-2
    if r2 is None:
        r2 = len(m2._vp)-2  
    n1 = len(m1._vp)
    n2 = len(m2._vp)
    if n2 == 0:
        m1.add_edge(r1, -1, v1_color=1)
    else:
        m1.disjoint_union(m2, check=False)
        m1.add_edge(r1, m1._fp[n1+r2])

    return m1
    


def peninsula_3col(m1, m2, r1=None, r2=None, check=True):
    r"""
    Add an edge to m2 so that its outer face has degree 2 and merge its root corner to the one of m1. If the root half edge is not provided it is assumed to be 2n-2. If m2 is a single vertex then m2 is replaced by an edge.

    INPUT:
        - ``m1``, ``m2``: red maps
        - ``r1``, ``r2``: their root half edges, if ``None`` then assumed to be 2ni-2
    
    EXAMPLES::

            sage: M = atomic_3colmap(0)
            sage: peninsula_3col(M, atomic_3colmap())
             ColoredOrientedMap("(0)(~0)", "(0,~0)", edge colors: [None], vertex colors: {0: 0, 1: 2})
            sage: M1 = ColoredOrientedMap("(0,1)(~0,~1)", "(0,~1)(~0,1)", vcolors: {-1: 0, 0: 0, 1: 2}, mutable = True)
            sage: M2 = ColoredOrientedMap("(0)(~0)", "(0,~0)", vcolors: {0: 0, 1: 2})
            sage: peninsula_3col(M1, M2)
             ColoredOrientedMap("(0,1,2,3)(~0,~1)(~2,~3)", "(0,~1)(~0,3,~2,1)(2,~3)", 
             edge colors: [None, None, None, None], vertex colors: {-1: 0, 1: 2, 5: 2, 0: 0})
    """

    if check:
        m1._assert_mutable()
        m1._check()
        m2._check()
    
    if r1 is None :
        r1 = len(m1._vp)-2
    if r2 is None:
        r2 = len(m2._vp)-2  
    n1 = len(m1._vp)
    n2 = len(m2._vp)
    if n2 == 0:
        m1.add_edge(r1, -1, v1_color=2)
    else:
        m1.disjoint_union(m2, check=False)
        if n1 > 0:
            h = m1.next_in_face(n1+r2)
            m1.merge_vertices(r1, n1+r2)
            m1.add_edge(n1+r2, h)
        else:
            m1.add_edge(r2, m1.next_in_face(r2))

    return m1


def sew_face_3col(mgray, mred, cw, rgray=None, rred = None, check=True):
    r"""
    Add an edge from the corner ``rgray`` to close a face on ``mgray`` and sew ``mred`` inside as described in ``cw``.

    INPUT:
        - ``mgray``: a gray map
        - ``rgray``: its root half edge, if ``None`` then assumed to be 2n-2
        - ``mred``: a red map
        - ``rred``: its root half edge, if ``None`` then assumed to be 2n-2
        - ``cw``: a list of integers encoding how to glue mred. ``cw[i]`` is the number of corners to skip between le ith and the i+1th edge to glue.

    EXAMPLES::

            sage: M = ColoredOrientedMap(vp="(5, ~0)(0, ~1)(1,~2)(2,~3)(3,~4)(4,~5)", vcolors={0:1, 2:0, 4:1, 6:0, 8:1, 10:0}, mutable =True)
            sage: H = ColoredOrientedMap(vp="(0,1,2)(~0,~2,~1)", vcolors={0:0, 1:2})
            sage: sew_face_3col(M, H, [0, 2, 1, 0])
             ColoredOrientedMap("(0,~1)(~0,5,8,9)(1,~2)(2,~9,6,~3)(3,7,~4)(4,~5)(~6,~8,~7)", "(0,9,2,1)(~0,~1,~2,~3,~4,~5)(3,6,~7)(4,7,~8,5)(~6,~9,8)", 
             edge colors: [None, None, None, None, None, None, None, None, None, None], 
             vertex colors: {0: 1, 1: 0, 2: 0, 4: 1, 6: 0, 8: 1, 13: 2})
    """

    if check:
        mgray._assert_mutable()
        mgray._check()
        mred._check()
    
    ngray = len(mgray._vp)
    nred = len(mred._vp)
    if rgray is None:
        rgray = ngray -2
    if rred is None:
        rred = nred -2
        
    if check:    
        if mgray.face_degree(rgray) < sum(cw)+1:
            raise ValueError("The face to create is bigger than the outer face.")
        if mred.vertex_degree(rred) != len(cw)-1:
            raise ValueError("The length of the connection word doesn't match the degree of the root vertex.")
        if sum(cw) % 2 == 0:
            raise ValueError("The sum of the moves in cw should be odd.")
    
    mgray.disjoint_union(mred, check=False)
    if cw[0] > 0:
        new_root = rgray
    else:
        new_root = ngray+rred
    cur_corner = rgray
    cur_half_edge = ngray+rred
    for l in cw[:-1]:
        for i in range(l):
            cur_corner = mgray.next_in_face(cur_corner)
        newh = mgray.previous_at_vertex(cur_half_edge)
        mgray.move_half_edge(cur_half_edge, cur_corner)
        cur_half_edge = newh
    for i in range(cw[-1]):
        cur_corner = mgray.next_in_face(cur_corner)
    mgray.add_edge(new_root, cur_corner)

    return mgray




def expand_vertex_3col(mred, mgray, cw,  rred=None, rgray=None, swap=False, shift=False, check=True):
    r"""
    Splits the root vertex of ``mred`` in two and expands one half into ``mgray``, positionning the edge as described in ``cw``

    INPUT:
        - ``mred``: a red map
        - ``rred``: its root half edge, if ``None`` then assumed to be 2n-2
        - ``mgray``: a gray map, if ``None`` then assumed to be 2n-2
        - ``rgray``: its root half edge
        - ``cw``: a list of integers encoding how to expand the new vertex into mgray. cw[i] is the number of edges to glue in the ith corner.
        - ``swap``: if True then the colors 0 and 1 will be swapped in mgray.
        - ``shift``: if True then the root corner of ``mgray`` will be changed for the next corner in the outer face. If mgray is atomic then the new vertex will have color ``1``instead of ``0``. (faster alternative to swap for random sampling)

    EXAMPLES::

            sage: H = ColoredOrientedMap(vp="(0,1,2,3)(~0,~1)(~3,~2)", vcolors={0:0,1:2,5:2}, mutable=True)
            sage: expand_vertex_3col(H, atomic_3colmap(), [2])
             ColoredOrientedMap("(0,4)(~0,~1,~4)(1,2,3)(~2,~3)", "(0,~4)(~0,4,~1,3,~2,1)(2,~3)", 
             edge colors: [None, None, None, None, None], vertex colors: {0: 0, 1: 2, 5: 2, 2: 0})

            sage: M = ColoredOrientedMap(vp="(~0,1)(~1,2)(~2,3,4)(~3,0)(~4}", vcolors={0:1,2:0,4:1,6:0,9:1}, mutable=True)
            sage: H = ColoredOrientedMap(vp="(0,1,2,3)(~0,~1)(~3,~2)", vcolors={0:0,1:2,5:2}, mutable=True)
            sage: expand_vertex_3col(H.copy(), M, [0,1,0,1,0,0,0])
             ColoredOrientedMap("(0,9)(~0,~1,~9)(1,~6,7,2,8)(~2,~3)(3,~8)(4,~7)(~4,5)(~5,6)", "(0,~9)(~0,9,~1,8,3,~2,7,4,5,6,1)(2,~3,~8)(~4,~7,~6,~5)", 
             edge colors: [None, None, None, None, None, None, None, None, None, None], 
             vertex colors: {0: 0, 1: 2, 5: 2, 8: 1, 9: 0, 11: 1, 6: 1, 2: 0})
            sage: expand_vertex_3col(H.copy(), M, [0,1,0,1,0,0,0], shift=True)
             ColoredOrientedMap("(0,9)(~0,~1,~9)(1,~8)(2,~7,4)(~2,~3)(3,8,~6,7)(~4,5)(~5,6)", "(0,~9)(~0,9,~1,~8,3,~2,4,5,6,8,1)(2,~3,7)(~4,~7,~6,~5)", 
             edge colors: [None, None, None, None, None, None, None, None, None, None], 
             vertex colors: {0: 0, 1: 2, 5: 2, 9: 0, 11: 1, 6: 0, 4: 1, 2: 1})
    """

    if check:
        mred._assert_mutable()
        mgray._check()
        mred._check()

    ngray = len(mgray._vp)
    nred = len(mred._vp)
    if rgray is None:
        rgray = ngray -2
    if rred is None:
        rred = nred -2
    if swap:
        mgray.swap_color(check=check)
    
    if check:
        if mred.vertex_degree(rred) < sum(cw)+1:
            raise ValueError("The vertex to create is bigger than the root vertex.")
        if mgray.face_degree(rgray) != len(cw)-1:
            raise ValueError("The length of the connection word doesn't match the degree of the outer face.")

    empty_deg = mred.vertex_degree(rred) == sum(cw)+1
    
    cur_half_edge = rred
    if ngray >0:
        mred.disjoint_union(mgray, check=False)
        cur_corner = rgray+nred
        if shift:
            cur_corner = mred.next_in_face(cur_corner)
        for l in cw[:-1]:
            for i in range(l):
                nexth = mred.previous_at_vertex(cur_half_edge)
                mred.move_half_edge(cur_half_edge, cur_corner)
                cur_half_edge = nexth
            cur_corner = mred.next_in_face(cur_corner)
        for i in range(cw[-1]):
            nexth = mred.previous_at_vertex(cur_half_edge)
            mred.move_half_edge(cur_half_edge, cur_corner)
            cur_half_edge = nexth
        nexth = mred.previous_at_vertex(cur_half_edge)
        mred.move_half_edge(cur_half_edge, cur_corner)
        if empty_deg:
            mred.add_edge(-1, mred.previous_in_face(cur_corner), v0_color=0)
        else:
            mred.add_edge(nexth, mred.previous_in_face(cur_corner))
        
    else:
        cur_corner = cur_half_edge
        nexth = mred.previous_at_vertex(cur_half_edge)
        if swap or shift:
            mred.move_half_edge(cur_half_edge, -1, v_color = 1)
        else:
            mred.move_half_edge(cur_half_edge, -1, v_color = 0)    
        for i in range(cw[0]):
            cur_half_edge = nexth
            nexth = mred.previous_at_vertex(cur_half_edge)
            mred.move_half_edge(cur_half_edge, cur_corner)
        if empty_deg:
            mred.add_edge(-1, mred._ep(cur_half_edge), v0_color=0)
        else:
            mred.add_edge(nexth, mred._ep(cur_half_edge))    

    return mred
 