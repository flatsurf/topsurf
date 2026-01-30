r"""
This file provides the stadart operations that appear in the recursive decomposition of planar maps.
"""

from topsurf import *

def atomic():
    OrientedMap(vp=[], mutable = True)

def bridge(m1, m2, r1=None, r2=None, check=True):
    r"""
    Add a copy of ``m2`` in ``m1`` and joins them with an edge from the root corner of ``m1`` to the co-root corner of ``m2``. 

    INPUT:
        - ``m1`, ``m2``: maps
        - ``r1``, ``r2``: their root half edges, if ``None`` then assumed to be 2ni-2

    EXAMPLES::

            sage: from topsurf import OrientedMap
            sage: M1 = OrientedMap(vp = [3, 4, 5, 0, 1, 2], mutable=True)
            sage: M2 = OrientedMap(vp = [2, 3, 0, 1])
            sage: bridge(M1, M2)
            sage: M1
            OrientedMap("(0,~1)(~0,2,5)(1,~2)(3,4)(~3,~4,~5)", "(0,5,~4,3,~5,2,1)(~0,~1,~2)(~3,4)")
            sage: M0 = OrientedMap(vp = 0)
            sage: bridge(M1, M0)
            sage: 0rientedMap("(0,~1)(~0,2,5,6)(1,~2)(3,4)(~3,~4,~5,~6)", "(0,6,~5,2,1)(~0,~1,~2)(3,~6,5,~4)(~3,4)")
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
        m1.add_edge(r1, -1)
    else:
        m1.disjoint_union(m2, check=False)
        m1.add_edge(r1, m1._fp[n1+r2])


def peninsula(m1, m2, r1=None, r2=None, check=True):
    r"""
    Add an edge to ``m2`` so that its outer face has degree 2 and merge its root corner to the one of ``m1``. If the root half edge is not provided it is assumed to be 2n-2. If ``m2`` is a single vertex then ``m2`` is replaced by an edge.

    INPUT:
        - ``m1``, ``m2``: maps
        - ``r1``, ``r2``: their root half edges, if ``None`` then assumed to be 2ni-2
    
    EXAMPLES::

            sage: from topsurf import OrientedMap
            sage: M1 = OrientedMap(vp = [3, 4, 5, 0, 1, 2], mutable=True)
            sage: M2 = OrientedMap(vp = [2, 3, 0, 1])
            sage: peninsula(M1, M2)
            sage: M1
            OrientedMap("(0,~1)(~0,2,5,3,4)(1,~2)(~3,~5,~4)", "(0,4,~5,2,1)(~0,~1,~2)(3,~4)(~3,5)")
            sage: M0 = OrientedMap(vp = [])
            sage: peninsula(M1, M0)
            sage: M1
            OrientedMap("(0,~1)(~0,2,5,6,3,4)(1,~2)(~3,~5,~4)(~6)", "(0,4,~5,2,1)(~0,~1,~2)(3,~4)(~3,6,~6,5)")
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
        m1.add_edge(r1, -1)
    else:
        m1.disjoint_union(m2, check=False)
        if n1 > 0:
            h = m1.next_in_face(n1+r2)
            m1.merge_vertices(r1, n1+r2)
            m1.add_edge(n1+r2, h)
        else:
            m1.add_edge(r2, m1.next_in_face(r2))


def close_face(m, k, r=None, check=True):
    r"""
    Add a edge between ``r`` and the (k-1)th next corner in its face, creating a face of degree k. Raises an error is the face has degree less than k-1.

    INPUT:
        - ``m``: a map
        - ``k``: the degree of the face to create
        - ``r``: the root half edge, if ``None`` then assumed to be 2n-2

    EXAMPLES::

            sage: M = OrientedMap(vp=[0, 2, 1, 4, 3, 5], mutable=True)
            sage: close_face(M, 4)
            sage: M
            OrientedMap("(0)(~0,~3,1)(~1,2,3)(~2)", "(0,1,3,~0)(~1,~3,2,~2)")
    """

    if check:
        m._assert_mutable()
        m._check()
        if k < 1:
            raise ValueError("k < 1")

    n = len(m._vp)
    if r is None:
        r = n -2
    if n == 0:
        if k>1:
            raise ValueError("k is larger than the degree of the face")
        else:
            m.add_edge(-1,-1)
    if k==1:
        m.add_edge(r, r)
        m.reverse_orientation(n//2)
    else:
        c = r
        for i in range(k-1):
            c = m.next_in_face(c)
            if c == r and i != k-2:
                raise ValueError("k is larger than the degree of the face")
        m.add_edge(r, c)



def split_vertex(m, k, r=None, check=True):
    r"""
    Splits the root vertex and the kth edge around it from ``r`` in two, creating a vertex of degree ``k`` and a new edge.

    INPUT:
        - ``m``: a map
        - ``k``: the degree of the certex to create
        - ``r``: the root half edge, if ``None`` then assumed to be 2n-2

    EXAMPLES::

            sage: M = OrientedMap(vp= [2, 1, 4, 3, 0, 5], mutable=True)
            sage: split_vertex(M, 2)
            sage: M
            OrientedMap("(0,3)(~0)(1,2)(~1,~3)(~2)", "(0,~0,3,~1,2,~2,1,~3)")
    """
    
    if check:
        m._assert_mutable()
        m._check()
        if k < 1:
            raise ValueError("k < 1")

    n = len(m._vp)
    if r is None:
        r = n -2
    if n == 0:
        raise ValueError("Cannot split a single vertex")
    else:
        h0 = m.previous_in_face(r)
        c = r
        k_max = False
        for i in range(k):
            c = m.previous_at_vertex(c)
            if c == r:
                if i == k-1:
                    k_max = True
                else:
                    raise ValueError("k is larger than the degree of the vertex")
        if k_max:
            m.add_edge(-1, h0)
        else:
            h1 = m.previous_in_face(c)
            m.insert_edge(h0, h1)
            m.move_half_edge(n+1, h1)
            
        
        






    