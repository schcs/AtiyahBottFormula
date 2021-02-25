load( "treesGiosue.sage" )
class ModuliStableCurve(SageObject):
    
    def __init__(self, n, d):
        self._dim = (n+1)*d+n-3
        self._d = d
        self._n = n
        self._num_col = n+1
        pass
    
    def _repr_(self): 
        #Returns a string representation of this moduli space.
        return "The moduli spaces of degree d stable maps from a genus 0 curve to a projective space of dimension n."   
    def dim(self):
        #Returns the dimension of this moduli space.
        return self._dim
    def degree(self):
        #Returns the degree of the maps parameterized by this moduli space.
        return self._d
    def dim_proj_space(self):
        #Returns the dimension of the projective space.
        return self._n
    def fixed_points(self):
        #Returns the number of fixed points of the torus action on the projective space.
        return self._num_col
    
    
    def fixed_locus(self):
        #Returns the fixed point loci, **with duplicades**, of this moduli space under the natural torus action.
        grs = []
        for num_v in list(IntegerRange(2,self._d+2)):
            tr = graphs.trees(num_v)
            for t in tr:
                cols = sage.graphs.graph_coloring.all_graph_colorings(t,self._num_col)
                ed = t.edges()
                num_ed = len(ed)
                for col in cols:
                    for c in list(IntegerRange(0,self._num_col)):
                        if not c in col:
                            col[c] = []
                    for p in Partitions( self._d, length = num_ed ):
                        for w in Arrangements( p, len(p)):
                            e_w = [[ed[i][0],ed[i][1],w[i]] for i in range(len(w))]
                            gr = ColoredEdgeWeightedGraph( t, col, e_w)
                            grs.append(gr)
        return grs

    
    def fixed_locus_pure(self):
        #Returns the fixed point loci, **without duplicades**, of this moduli space under the natural torus action.
        grs = []
        for num_v in list(IntegerRange(2,self._d+2)):
            tr = graphs.trees(num_v)
            for t in tr:
                cols = sage.graphs.graph_coloring.all_graph_colorings(t,self._num_col)
                ed = t.edges()
                num_ed = len(ed)
                for col in cols:
                    for c in list(IntegerRange(0,self._num_col)):
                        if not c in col:
                            col[c] = []
                    for p in Partitions( self._d, length = num_ed ):
                        for w in Arrangements( p, len(p)):
                            e_w = [[ed[i][0],ed[i][1],w[i]] for i in range(len(w))]
                            gr = ColoredEdgeWeightedGraph( t, col, e_w)
                            t=0
                            for temp in grs:
                                #print("a")
                                if gr.is_isomorphic( temp ):
                                    t=1
                                    break
                            if t==0:
                                grs.append(gr)
        return grs
    
    def fixed_locus_pure2(self):
        #Returns the fixed point loci, **without duplicades**, of this moduli space under the natural torus action.
        grs = []
        for gr in self.fixed_locus():
            t=0
            for temp in grs:
                if gr.is_isomorphic( temp ):
                    t=1
                    break
            if t==0:
                grs.append(gr)
        return grs
            
    
            
                    
            
            
            