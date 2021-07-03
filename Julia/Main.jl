#module  AtiyahBott

using ProgressMeter
using LightGraphs
using Combinatorics

<<<<<<< HEAD
include("Marked.jl")
=======
include( "Marked.jl" )
>>>>>>> 7a72e15ec329d7dc7850bea9b1fecf423f43e442
include("GraphFunctions.jl")
include("EquivariantClasses.jl")
include("Checks.jl")

#export AtiyahBottFormula, AtiyahBottFormulaForGraph

number_trees = [1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]
#the number of non-isomorphic graphs with given number of vertices (starting from 2)

# ProgressData contains the data necessary to keep the progress bar up-to-date
mutable struct ProgressData
    progress_bar::Progress
    current_graph::Int64
    top_aut::Int64
    threshold::Int64
end

# create trivial progress data in case the user does not want progress bar
function EmptyProgressData()::ProgressData
    P::Progress = Progress( 0 )
    P.enabled = false
    return ProgressData( P, 0, 0, 0 )
end

# the following function performs the computation of the Atiyah-Bott formula
# for a particular graph.
function AtiyahBottFormulaForGraph( g::SimpleGraph, pruf_str::String, 
    aut::Int64, n::Int64, deg::Int64, n_marks::Int64, P, s::Vector{Rational{BigInt}},
    progress_data::ProgressData)::Vector{Rational{BigInt}}

    local weights::Vector{Vector{Int64}} = get_weights(nv(g)-1, deg) #the array of array of weights

    local n_results::Int64 = length(P)   #this store the number of final results of our computation
    local result::Vector{Rational{BigInt}} = [Rational{BigInt}(0) for _ in 1:n_results] #the array of final results
    
    
    try_to_find_color_file = true #IS THIS USEFUL???
    if try_to_find_color_file   #we always try to read the colorations from the a file where they are unique modulo isomorphisms
        from_file, file_name = exists_file_with_colorings( pruf_str, n )

        if from_file 
            cols = graph_coloring_from_file( file_name )        
        else 
<<<<<<< HEAD
            #println( "graph not found!!!" )
=======
            println( "graph not found!!!" )
>>>>>>> 7a72e15ec329d7dc7850bea9b1fecf423f43e442
            cols = graph_coloring( g, UInt8(n+1) )
        end 
    else 
        cols = graph_coloring( g, UInt8(n+1) )
        from_file = false
    end

    for c in cols   #we run among all colorations of g
        if from_file 
            aut = cols.current_aut  #we are read the number of automorphisms of the colorated graph from the file
        end
        
        for m in marks(nv(g),n_marks)    #we run among all marks of g, if n_marks==0 we have only the empty mark
            for w in weights          #we run among all weights of g
                try
                    local Euler::Rational{BigInt} = Euler_inv(g,c,w,s,m)//(aut*prod(w)) #the contribuition of the Euler class in the Atiyah-Bott formula
                    for res in 1:n_results      #compute each term of the array P
                        result[res] += P[res](g,c,w,s,m)*Euler    #apply Atiyah-Bott
                    end
                catch err 
                    if isa(err, DivideError) 
                        error("Some division by zero occurred. Try again")
                    end
                    println(err)
                    error("Some error occurred")
                    return [Rational{BigInt}(0)]
                end
            end

            if progress_data.progress_bar.enabled
            
                progress_data.current_graph += progress_data.top_aut÷aut   
<<<<<<< HEAD
                #update the progress bar
=======
                #upgrade the progress bar
>>>>>>> 7a72e15ec329d7dc7850bea9b1fecf423f43e442
                update!(progress_data.progress_bar, 
                        progress_data.current_graph,
                        showvalues = [(:"Total number of graphs",progress_data.threshold),
                        (:"Current graph",progress_data.current_graph)])
            end
        end
    end
    return result
end 


"""
    AtiyahBottFormula(n, d, m, P; do_check)

Apply the Atiyah-Bott residue formula to the class `P`, in the moduli space of rational marked stable maps to the projective space of dimension `n` of degree `d` with `m` marks.
# Arguments
- `n::Int64`: the dimension of the projective space.
- `d::Int64`: the degree of the stable maps.
- `m::Int64`: the number of marks.
- `P`: the equivariant class.
- `do_check::Bool`: if `true`, checks if `P` is a well defined zero cycle, and stops the computation if this is not true. If `false`, the computation may have an unexpected behaviour. By default is `true`.
- `try_to_find_color_file::Bool`: if `true` will read the coloration from the folder Data. If `false` or the file are not present, it will generate the colorations internally. By default is `true`. 

The general construction of `P` is the following:
```jldoctest
julia> P = (g,c,w,s,m) ->
```
After `->`, one has to write an expression in the equivariant classes. All such equivariant classes are functions starting with `(g,c,w,s)` or `(g,c,w,s,m)`. At the end, they can have more arguments. The expression is a polynomial combination of the equivariant classes. We compute the degree of `P` by

```jldoctest
julia> AtiyahBottFormula(n,d,m,P);
```

# Example
```jldoctest
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);
julia> AtiyahBottFormula(3,1,0,P);
Warning: the class is not a 0-cycle.
julia> AtiyahBottFormula(4,1,0,P);
Result: 2875//1
```

The function returns an array of the same dimension of `P` (non-vectorized classes are assumed as 1-dimensional arrays). The Julia notation for accessing to array is `name_of_array[i]` where `i` is an index starting from 1.

# Example
```jldoctest
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)*Hypersurface(g,c,w,s,3);
julia> x = AtiyahBottFormula(3,2,0,P)[1];
Result: 81//1
julia> x
81//1
```

The class `P` supports parameters.
```jldoctest
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,3)*(Incidency(g,c,w,s,2)//3)^(d-1);
julia> d = 2;
julia> AtiyahBottFormula(3,d,0,P);
Result: 27//1
julia> d = 3;
julia> AtiyahBottFormula(3,d,0,P);
Result: 84//1
```

More examples are available in the support of the equivariant classes. It is enough to type `?` and then the name of the class. Currently, the supported classes are:

* `O1_i`
* `O1`
* `Incidency`
* `Hypersurface`
* `Contact`
* `R1`
* `Psi`
* `Jet`
"""
function AtiyahBottFormula(n::Int64, deg::Int64, n_marks::Int64, P; do_check::Bool = true)::Vector{Rational{BigInt}}
    
    try_to_find_color_file::Bool = true
    if n < 1 || deg < 1
        printstyled("ERROR: ", bold=true, color=:red)
        println("n and d must be positive, correct ", n, " and ", deg)
        return [Rational{BigInt}(0)]
    end
    if n_marks < 0
        printstyled("ERROR: ", bold=true, color=:red)
        println("m must be non negative, correct ", n_marks)
        return [Rational{BigInt}(0)]
    end
    if !isa(P, Array)  #we want that P is an array, possibly with only one element
        P = [P]
    end

    if do_check
        if !is_zero_cycle(n, deg, n_marks, P)
            return [Rational{BigInt}(0)]
        end
    end
    
    local n_results::Int64 = length(P)   #this store the number of final results of our computation
    local result::Vector{Rational{BigInt}} = [Rational{BigInt}(0) for _ in 1:n_results] #the array of final results
    local max_col::Int64 = n+1   #the colors are number from 1 to n+1
     
    list_g::IOStream = open( "list_trees.txt" ) 
    #open the file containing the list of Prufer sequences of graphs
    
    s::Vector{Rational{BigInt}} = [ Rational{BigInt}(rand(-1000*max_col*deg:1000*max_col*deg)) 
                                        for i in 1:max_col] 

    #set up progress data                                    
    threshold::Int64 = sum([number_trees[v-1]*max_col*(n^(v-1))*(v^n_marks) 
                for v in 2:deg+1])


    progress_data = ProgressData( 
        Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green), #progress_bar
        0, #current_graph
        0, #top_aut
        threshold ) #threshold

    for v in 2:(deg+1) #run the computation among all graphs with fixed number of vertices
                
        n_trees_nv = number_trees[v - 1]  #we known how many graphs there are with fixed number of vertices
        
        for n_g in 1:n_trees_nv  #run the computation for a fixed graph

            str = readline(list_g) #read a new line, we expect a Prufer seq plus the number of automorphisms
            local (g, aut) = get_graph(str)  #g is the graph, aut is the number of automorphisms
            progress_data.top_aut = aut
            #local top_aut::Int64 = copy(aut) #the number of automorphisms of g without any decoration

            pruf_str = string(split( str, ',')[1])
            result += AtiyahBottFormulaForGraph( g, pruf_str, aut, n, deg, 
                    n_marks, P, s, 
                    progress_data)

        end
    end
    
    close(list_g)     #close the file with Prufer sequences
    if n_results == 1
        println("Result: ",result[1])
    else 
        for res in 1:n_results
            println("Result number ",res,": ",result[res])
        end
    end
    return result
end
<<<<<<< HEAD
=======

function check_Data(data_dir = "..")
    
    data_dir = data_dir*"/Data/"
    
    if !isdir(data_dir)
        println("Folder ", data_dir, " not found.")
        return
    end
    
    Dimension_dirs = [x for x in readdir(data_dir) if startswith(x,"Dimension")]
    
    if length(Dimension_dirs) == 0
        println("""Folder "Data" is empty.""")
        return
    end
    
    println("""Folder "Data" found.""")
    
    for current_dir in Dimension_dirs
        println(current_dir," contains:")
        files = [x for x in readdir(data_dir*current_dir)]
        for v in 2:14
            num = count(x->length(x)==v+3 && endswith(x,".clr"), files)
            if num == 0
                println( "No colored graphs with ", v, " vertices.")
                continue
            elseif num < number_trees[v-1]
                println("Some colored graph with ",v," vertices is missing.")
            else
                println("All colored graph with ",v," vertices.")
            end
            
        end
    end
    
    return
end
#end
>>>>>>> 7a72e15ec329d7dc7850bea9b1fecf423f43e442
