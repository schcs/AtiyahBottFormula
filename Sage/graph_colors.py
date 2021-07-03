from sympy.combinatorics.prufer import Prufer 
from sage.misc.search import search 

def write_graph_list_for_julia( grlist, filename, varname ):

    fn = open( filename, 'w' )
    fn.write( varname + ' =  [\n' )
    for gr in grlist:
        fn.write( write_graph_string_for_julia(gr)+ '\n' )

    fn.write( '];\n' )

#from sage.misc.search import search

def color_dict_to_color_list( cdict ):

    return [ cdict[k] for k in [0..max( cdict.keys())]]


def extend_graph( gr, col, nr_cols ):

    gr_new = copy( gr )
    col_new = copy( col )
    nv = gr_new.num_verts()-1

    for i in range( nr_cols+1 ):
        if not i in col_new.keys():
            gr_new.add_vertex()
            nv += 1 
            col_new[i] = [nv]

    return gr_new, col_new

def add_loops( gr, col ):

    new_g = Graph( gr, weighted = True )
    new_g.allow_loops( true )

    for e in new_g.edges():
        new_g.set_edge_label( e[0], e[1], -1 )

    for k in col.keys():
            new_g.add_edge( k, k, col[k] )
    
    return new_g

def canonical_coloring_label( gr, col ):
    return add_loops( gr, col ).canonical_label( edge_labels = true )

def canonical_coloring_label_1(G,c):
    """Given a coloring dictionary,
    {color1 : [u1, u2, ...], color2 : [v1, v2, ... ], ... }
    return a string which uniquely identifies the isomorphism
    class of the coloring."""
    
    H = G.copy()
    #H.allow_loops( true )

    for i in c:
        print( i )
        H.add_edges([(i,j) for j in c[i]])

    P = [G.vertices(), c.keys()]
    return H.canonical_label(partition=P)

def is_isomorphic_colorings( gr, col1, col2 ):

    # if the two colorings have different number of colors, then return false
    if len(col1.keys()) != len(col2.keys()):
        return False
    
    # if the color classes have different lengths then return false
    if sorted([ len(col1[i]) for i in col1.keys()]) != sorted([ len(col1[i]) for i in col1.keys()]):
                return False

    max_col = max( list( col1.keys()) + list( col2.keys()))

    gr_ext1, col1_new = extend_graph( gr, col1, max_col )
    gr_ext2, col2_new = extend_graph( gr, col2, max_col )

    keys = list( col1_new.keys())
    keys.sort()

    part1 = [ col1_new[i] for i in keys ]
    part2 = [ col2_new[i] for i in keys ]

    return gr_ext1.canonical_label( partition = part1 ) == gr_ext2.canonical_label( partition = part2 )


def is_isomorphic_colorings_2( gr, col1, col2 ):
        c_gr = gr.canonical_label()
        
        a = gr.automorphism_group( edge_labels = False )
        p1 = tuple( set( col1[k] ) for k in col1.keys())
        p2 = tuple( tuple( col2[k] ) for k in col2.keys())
        orb = a.orbit( p2, action = "OnTuplesSets" )
        return p1 in orb

def is_isomorphic_colorings_3( gr, col1, col2 ):

    # if the two colorings have different number of colors, then return false
    if set(col1.values()) != set(col1.values()):
        return False
    
    # if the color classes have different lengths then return false
    #if sorted([ len(col1[i]) for i in col1.keys()]) != sorted([ len(col1[i]) for i in #col1.keys()]):
    #            return False
        
    gr1 = add_loops( gr, col1 )
    gr2 = add_loops( gr, col2 )

    return gr1.is_isomorphic( gr2, edge_labels = true )

def is_isomorphic_colorings_4( gr, col1, col2 ):

    # if the two colorings have different number of colors, then return false
    if len(col1.keys()) != len(col2.keys()):
        return False
    
    # if the color classes have different lengths then return false
    if sorted([ len(col1[i]) for i in col1.keys()]) != sorted([ len(col1[i]) for i in col1.keys()]):
                return False
        
    c1 = canonical_coloring_label( gr, col1 )
    c2 = canonical_coloring_label( gr, col2 )

    return c1 == c2

def is_isomorphic_coloring_in_list( gr, col, col_list ):

    for c in col_list:
        if is_isomorphic_colorings_3( gr, col, c ):
            return True

    return False 

def colorings_up_to_isomorphism( gr, num_col ):

    can_list = []
    cols_list = []
    auts = []
    cols = sage.graphs.graph_coloring.all_graph_colorings( gr, num_col, 
                    vertex_color_dict = True )
    for col in cols:
        gr_with_loops = add_loops( gr, col )
        cr_with_loops = gr_with_loops.canonical_label( edge_labels = True )        
        cr_edge_set = [ list( x ) for x in cr_with_loops.edges()]
        cr_edge_set.sort()
        v, pos = search( can_list, cr_edge_set )

        if not v:
            a = gr_with_loops.automorphism_group( edge_labels = True ).order()
            can_list.insert( pos, cr_edge_set )
            cols_list.append( col )
            auts.append( a )
    
    return cols_list, auts

def write_graph_colorings( gr, num_col, fname ):
    letters = [ '1','2','3','4','5','6','7','8','9',
                'a','b','c','d','e','f','g','h','i','j','k',
                'l','m','n','o','p','q','r','s','t','u'];
    
    cols, auts = colorings_up_to_isomorphism( gr, num_col )
    f = open( fname, "w" )
    pr = [ letters[x] for x in 
                    Prufer( gr.edges( labels = false )).prufer_repr ]
    # write Prufer sequence as a string 
    f.write( "".join( pr ) + "\n")

    for i in range( len( cols )):
        col_st = [ letters[cols[i][x]] for x in range( len( gr ))]
        f.write( "".join( col_st ))
        f.write( "," + str( auts[i]) + "\n" )
    
    f.write( "STOP\n" )
    f.close()
