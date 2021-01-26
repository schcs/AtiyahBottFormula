class ColoredEdgeWeightedGraph(sage.graphs.graph.Graph):
    """a new class for vertex colored and edge weighted graphs"""
    
    def __init__(self, *args, **kwargs ):
        
        if len(args) == 1 and type(args[0]) == ColoredEdgeWeightedGraph:
            super().__init__(args[0])
            ver_colors, edge_weights = args[0].vertex_colors, args[0].edge_weights
        else:
            gr, ver_colors, edge_weights = args[0], args[1], args[2]
            super().__init__( gr )
        
        self.edge_weights = edge_weights
        self.vertex_colors = ver_colors
        for w in edge_weights:
            self.set_edge_label( w[0], w[1], w[2] )

        num_vars = max( ver_colors.keys())+1
        self.fraction_field = PolynomialRing( QQ, 'x', 
               num_vars ).fraction_field()       

    def vertex_color( self, v ):
        for k in self.vertex_colors:
            if v in self.vertex_colors[k]:
                return k
                
    
    def automorphism_group( self ):
        parts = [ self.vertex_colors[k] for k in self.vertex_colors ]
        g = super().automorphism_group( edge_labels = True, partition = parts )
 
        return g

    def canonical_label( self ):
        g, a = super().canonical_label( edge_labels = True, certificate = True )
        max_cols = max( self.vertex_colors.keys())
        cols = { y: [ a[x] for x in self.vertex_colors[y]] for y in [0..max_cols]}
        return  ColoredEdgeWeightedGraph( g, cols, g.edges()), a

    def is_isomorphic( self, gr ):
        c_self, ac = self.canonical_label() 
        c_gr, agr = gr.canonical_label()
        
        if c_self != c_gr:
            return False
        a = super(ColoredEdgeWeightedGraph,c_self).automorphism_group( 
                            edge_labels = True )
        p1 = tuple( set(c_self.vertex_colors[k]) for k in c_self.vertex_colors.keys())
        p2 = tuple( tuple(c_gr.vertex_colors[k]) for k in c_gr.vertex_colors.keys())

        orb = a.orbit( p2, action = "OnTuplesSets" )
        return p1 in orb
    def function_V( self ): 

        k = self.fraction_field
        vcols = [ self.vertex_color( x ) for x in self.vertices()]
        maxcols = max( self.vertex_colors.keys()) 
        gens = k.gens()
        lam = lambda v:  gens[vcols[v]]

        
        return  product([ product( [ lam(v)-gens[j] for j in [0..maxcols] if 
                        j != vcols[v]])^(len(self.neighbors(v))-1)*
                sum( [ self.edge_label( v, w )/(lam(v)-lam(w)) for 
                        w in self.neighbors(v)])^(len(self.neighbors(v))-3)*
                product( [self.edge_label( v, w )/(lam(v)-lam(w)) for 
                            w in self.neighbors(v)]) for v in self. vertices()])        
            
        
    def function_E( self ):
            
        k = self.fraction_field
        I = k.one()
        Z = k.zero()
        vcols = [ self.vertex_color( x ) for x in self.vertices()]
        gens = k.gens()
        lam = lambda v:  gens[vcols[v]]
        
        nrcols = max( self.vertex_colors.keys() )
        prod = I
        for e in self.edges():
            v, w, de = e[0], e[1], e[2]
            # first factor
            fac1 = (-1)^de*(de/(lam(v)-lam(w)))^(2*de)/factorial(de)^2
        
            # second factor (double product)
            fac2 = I 
            for col in [0..nrcols]:
                if col != vcols[v] and col != vcols[w]:
                    lam0 = gens[col]
                    for alph in range( de+1 ):
                        fac2 *= ((alph*lam(v)+(de-alph)*lam(w))/de-lam0)^-1
            prod *= fac1*fac2

        return prod
    
    def function_Theorem43( self, r ):
        
        k = self.fraction_field
        vcols = [ self.vertex_color( x ) for x in self.vertices()]
        gens = k.gens()
        lam = lambda v: gens[vcols[v]]

        vcols = [ self.vertex_color( x ) for x in self.vertices()]

        return sum([ e[2]*sum([ lam(e[0])^t*lam(e[1])^(r-t) for t in [0..r]]) 
                    for e in self.edges()])


def is_isom_in_list( list, gr ):
    for g in list:
        if gr.is_isomorphic( g ):
            return True
    
    return False

def isom_types_graphs( num_v, #number of vertices 
                       num_col,  #number of colors
                       sum_weight ): # sum of edge weights

    tr = graphs.trees(num_v)
    grs = []

    for t in tr:    
        pos = len(grs)
        cols = sage.graphs.graph_coloring.all_graph_colorings(t,num_col)
        ed = t.edges()
        num_ed = len(ed)
        for col in cols:
            for c in [0..num_col-1]:
                if not c in col:
                    col[c] = []
        
            for p in Partitions( sum_weight, length = num_ed ):
                for w in Arrangements( p, len(p)):                    
                    e_w = [[ed[i][0],ed[i][1],w[i]] for i in range( len( w ))]
                    gr = ColoredEdgeWeightedGraph( t, col, e_w )
                    if not is_isom_in_list( grs[pos:len(grs)], gr ):
                        grs.append( ColoredEdgeWeightedGraph( t, col, e_w))
                

    return grs




# the examples of Giosuè
#gr1 = ColoredEdgeWeightedGraph( Graph( { 0:[1]}), {0:[0],1:[1], 2:[], 3:[]}, 
#            [[0,1,3]] )
#
#gr2 = ColoredEdgeWeightedGraph( Graph( { 0:[1], 1:[2] }), {0:[0],1:[1], 2:[2],3:[]}, 
#            [[0,1,2],[1,2,1]] )
#
#gr3 = ColoredEdgeWeightedGraph( Graph( { 0:[1], 1:[2] }), {0:[0,2],1:[1], 2:[],3:[]}, 
#            [[0,1,2],[1,2,1]] )
#
#gr4 = ColoredEdgeWeightedGraph( Graph( { 0:[1], 1:[2], 2:[3] }), 
#            {0:[0,2],1:[1,3], 2:[],3:[]}, [[0,1,1],[1,2,1],[2,3,1]] )
#
#gr5 = ColoredEdgeWeightedGraph( Graph( { 0:[1], 1:[2], 2:[3] }), 
#            {0:[0,2],1:[1], 2:[3],3:[]}, [[0,1,1],[1,2,1],[2,3,1]] )

#gr8 = ColoredEdgeWeightedGraph( Graph({1:[0,2,3]}),
#    {0:[0,2,3],1:[1],2:[],3:[]},[[0,1,1],[1,2,1],[1,3,1]])

#gr10 = ColoredEdgeWeightedGraph(g,{1:[0],2:[1],3:[2],0:[3]},[[0,1,1],[1,2,1],[1,3,1]])