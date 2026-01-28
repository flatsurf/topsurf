


from topsurf import OrientedMap



def labels(G):
    r"""
    Assign a label (an integer) to each vertex and each face of G.
    Compute for each half-edge the vertex and the face it belongs. 
    """
    vertices = G.vertices()
    faces = G.faces()
    hedge_to_vertex = [None for _ in G.half_edges()]
    hedge_to_face = [None for _ in G.half_edges()]
    for i in range(len(vertices)):
        for e in vertices[i]:
            hedge_to_vertex[e] = i
    for i in range(len(faces)):
        for e in faces[i]:
            hedge_to_face[e] = i
    return hedge_to_vertex, hedge_to_face


def tree_co_tree(G):
    r"""
    Compute a tree/co-tree decomposition of G.
    Return a list res of length the number of edges of G.
    res[e] is 0 if the edge e is in the tree, 1 if e is in the co-tree and 2 otherwise.
    """
    hedge_to_vertex, hedge_to_face = labels(G)
    vp = G.vertex_permutation(copy=False)
    T_found = [False for _ in G.vertices()]
    T_found[0] = True
    res = [2 for _ in G.edges()]
    e = 0
    path = [e]

    if not T_found[hedge_to_vertex[1]]:
        T_found[hedge_to_vertex[1]]=True
        res[0]=0
        path.append(1)
        e=vp[1]
    else:
        e=vp[0]
    
    while len(path) > 0:
        e1 = G._ep(e)
        if e == path[-1]:
            path.pop()
            e = vp[e1]
        elif not T_found[hedge_to_vertex[e1]]:
            T_found[hedge_to_vertex[e1]] = True
            res[e // 2] = 0
            path.append(e1)
            e = vp[e1]
        else:
            e = vp[e]

    fp = G.face_permutation(copy=False)
    C_found = [False for _ in G.faces()]
    C_found[0] = True
    e = fp[0]
    path = [0]
    
    while len(path) > 0:
        e1 = G._ep(e)
        if e == path[-1]:
            path.pop()
            e = fp[e1]
        elif res[e//2] == 0:
            e = fp[e]
        elif not C_found[hedge_to_face[e1]]:
            C_found[hedge_to_face[e1]] = True
            res[e//2] = 1
            path.append(e1)
            e = fp[e1]
        else:
            e = fp[e]
        
    return res


def tree_contraction(G, treecotree):
    r"""
    Compute the OrientedMap F obtained from G by:
        - contracting all the edges in the tree of treecotree
        - removing all the edges in the co-tree of treecotree
    Compute a list correspondence of length the number of half-edges of G such that correspondence[e] is:
        - None if e is an edge of the tree of treecotree
        - d such that fp[d]=e in G after contracting the edge of the tree if e is an edge of the cotree
        - the corresponding edge of F otherwise
    Compute a list recor such that recor[f] for f an half-edge in F is the corresponding half_edge in G
    Compute a list rank such that rank[e] is:
        - None if e is an edge of the tree of treecotree
        - the index of e around the vertex if e is an edge of the cotree
        - the total number of edge of the cotree between e and the next edge in F

    INPUT:

    treecotree should be a list of lenght the number of edges in G where G[e] is:
        - 0 if e is in the tree
        - 1 if e is in the co-tree
        - 2 otherwise
    """
    fp = G.face_permutation(copy=False)
    correspondence = [None for _ in G.half_edges()]
    rank = [None for _ in G.half_edges()]
    nfp = [None for _ in range(4*G.genus())]
    recor = [None for _ in range(4*G.genus())]
    i = 0
    while treecotree[i] != 2:
        i += 1
    last=0
    correspondence[2 * i] = 0
    correspondence[2 * i + 1] = 1
    recor[0] = 2 * i
    recor[1] = 2 * i + 1
    e = fp[2 * i]
    free_index = 2
    r = 0
    
    while e != 2 * i:
        e1 = G._ep(e)
        if treecotree[e // 2] == 1:
            correspondence[e] = last
            rank[e] = r
            r += 1
            e = fp[e1]
            
        elif treecotree[e//2] == 0:
            e = fp[e]
        else:
            if correspondence[e] == None:
                correspondence[e] = free_index
                correspondence[e1] = free_index+1
                recor[free_index] = e
                recor[free_index + 1] = e1
                free_index += 2
            rank[e] = r
            r = 0
            nfp[last] = correspondence[e]
            last = correspondence[e]
            e = fp[e]
    rank[2 * i] = r
    nfp[last] = 0
    F = OrientedMap(fp=nfp)
    return F, correspondence, recor, rank



class QuadSystem:

    #TODO : Documentation

    def __init__(self, G, treecotree=None, check=True):
        r"""
        Methods :
            _origin_map : the original OrientedMap
            _genus : the genus of the underlying surface
            _quad : the quad system
            _proj : the projection from _origin_map half-edges to path of length 2 of _quad
        """

        if check:
            G._check()
            G._assert_connected()
            if G.has_folded_edge():
                raise NotImplementedError
            if G.genus() == 0:
                raise ValueError("Cannot look for homotopy in genus 0")
            elif G.genus() == 1:
                raise ValueError("Quad systems and homotopy are not effective in genus 1")

        if G.is_mutable():
            G = G.copy(mutable=False)

        self._origin_map = G
        self._genus = G.genus()

        if treecotree == None:
            treecotree = tree_co_tree(G)
        
        F,cor, _, _ = tree_contraction(G,treecotree)
        vp = F.vertex_permutation(copy=False)
        fp = F.face_permutation(copy=False)
        qvp = [None for _ in range(8 * self._genus)]
        remember = [None for _ in range(4 * self._genus)]
        remember2 = [None for _ in range(4 * self._genus)]
        
        e = fp[0]
        last_edge = 0
        next_edge = 1
        remember[0] = vp[0]
        remember2[0] = 0
        while e != 0:
            qvp[2 * last_edge] = 2 * next_edge
            remember[next_edge] = vp[e]
            remember2[e] = next_edge
            e = fp[e]
            last_edge = next_edge
            next_edge += 1
        qvp[-2] = 0
        for ne in range(4*self._genus):
            qvp[2 * ne + 1] = 2 * remember2[remember[ne]] + 1
        self._quad = OrientedMap(vp=qvp)
        
        qfp = self._quad.face_permutation(copy=False)
        cor2 = []
        for e in G.half_edges():
            if treecotree[e//2] == 0:
                cor2.append([])
            elif treecotree[e//2] == 1:
                start = fp[cor[e]]
                e1 = G._ep(e)
                end = fp[cor[e1]]
                das = 2 * remember2[start] + 1
                dae = 2 * remember2[end]
                cor2.append([das,dae])
            else:
                edge = cor[e]
                da = 2 * remember2[edge] + 1
                cor2.append([da, qvp[da - 1]])
        self._proj = cor2





