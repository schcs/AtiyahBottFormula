load( "graph_colors.py" )
from sympy.combinatorics.prufer import Prufer 
from sage.misc.search import search 

def generate( n , min_deg, max_deg):

    letters = [ '1','2','3','4','5','6','7','8','9',
                'a','b','c','d','e','f','g','h','i','j','k',
                'l','m','n','o','p','q','r','s','t','u'];

    for v in range((min_deg+1),(max_deg+2)):
        tr = graphs.trees(v)
        for t in tr:
            #compute prufer sequence
            pr = [ letters[x] for x in 
                    Prufer( t.edges( labels = false )).prufer_repr ]
            filename = "../Data/Dimension"+str(n)+"/"+"".join(pr)+"0.clr"
            write_graph_colorings( t, n+1, filename )

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

    for i in range( len( cols )):
        col_st = [ letters[cols[i][x]] for x in range( len( gr ))]
        f.write( "".join( col_st ))
        f.write( "," + str( auts[i]) + "\n" )
    
    f.write( "STOP" )
    f.close()
