from sympy.combinatorics.prufer import Prufer 
from sage.misc.search import search 

def generate( n , min_deg, max_deg):
    
    for v in range((min_deg+1),(max_deg+2)):
        tr = graphs.trees(v)
        counter = 1
        for t in tr:
            write_graph_colorings( t, n+1, "name_"+str(n)+"_"+str(v)+"_"+str(counter) )
            counter += 1

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
                    'a','b','c','d','e','f','g','h','i','j','k'];
    
    cols, auts = colorings_up_to_isomorphism( gr, num_col )
    f = open( fname, "w" )
    pr = [ letters[x] for x in 
                    Prufer( gr.edges( labels = false )).prufer_repr ]
    # write Prufer sequence as a string 
    f.write( "".join( pr ) +","+ "\n")

    for i in range( len( cols )):
        col_st = [ letters[cols[i][x]] for x in range( len( gr ))]
        f.write( "".join( col_st ))
        f.write( "," + str( auts[i]) + "\n" )
    
    f.write( "STOP" )
    f.close()
