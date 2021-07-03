"""
The structure `graph_coloring` define the colorations of a graph.
# Arguments
- `graph::SimpleGraph`: the graph.
- `num_cols::UInt8`: the number of colors.
- `current_color::Array{UInt8,1}`: the array of colors.
"""
mutable struct graph_coloring
    graph::SimpleGraph
    num_cols::UInt8
    current_color::Array{UInt8,1}
end

function graph_coloring( G::SimpleGraph, num_cols::UInt8 )
    return graph_coloring( G, num_cols, [] )
end

mutable struct graph_coloring_from_file
    file_name::String
    color_file::IOStream
    current_color::Array{UInt8,1}
    current_aut::Int64
end

function graph_coloring_from_file( file_name::String )
    color_file = open( file_name )
    return graph_coloring_from_file( file_name, color_file, [], 1 )
end

function Base.iterate( GC::graph_coloring_from_file, c=0 )

    st = readline( GC.color_file )
    if st == "STOP" 
        close( GC.color_file )
        return nothing
    end 

    st = split( st, "," )
    GC.current_color = [ parse( UInt8, x ) for x in split( st[1], "" )]
    GC.current_aut = parse( Int64, st[2] )
    return GC.current_color, 0
end

function exists_file_with_colorings( pruf_str::String, 
        dim::Int64; data_dir = "../Data/" )

    dir = data_dir*"Dimension"*string( dim )*"/"

    if isdir( dir ) && pruf_str*"0.clr" in readdir( dir )
        return true, dir*pruf_str*"0.clr"
    else
        return false, nothing
    end
end

function is_coloring( G::SimpleGraph, cols::Array{UInt8,1} )

    num_v = nv( G )
    for i in 1:num_v
        for j in (i+1):num_v
            if cols[i] == cols[j] && has_edge( G, i, j ) 
                return false
            end
        end
    end

    return true
end

function smallest_coloring( G::SimpleGraph, num_cols::UInt8 )

    cols = [ 0x1 for _ in 1:nv(G)]
    while true
        if is_coloring( G, cols )
            return cols
        end

        next_tuple!( cols, num_cols )
    end
end


function Base.iterate( GC::graph_coloring, c=0 )

    if GC.current_color == []
        GC.current_color = smallest_coloring( GC.graph, GC.num_cols )
        return GC.current_color, 0
    end

    while true
        v = next_tuple!( GC.current_color, GC.num_cols )
        if v == false
            return nothing
        end

        if is_coloring( GC.graph, GC.current_color )
            return GC.current_color, 0
        end

    end
end
    

function next_tuple!( tuple::Array{UInt8,1}, max_entry::UInt8 )

    k = length( tuple )
    for k in k:-1:0

        if k == 0  
            return false        
        elseif tuple[k] < max_entry 
            tuple[k] += 1
            break
        else
            tuple[k] = 1
        end
    end

    return nothing
end


"""
Compute the graph with a specific Prufer sequence.
"""
function PruferToGraph(prufer::Vector{UInt8})::SimpleGraph
    
    num_v = length(prufer)+2 #number of vertices
    g = Graph(num_v) #initialize a graph with specific number of vertices
    
    degree = [count(i->(i==j),prufer) for j in 1:num_v ] #new version of degree

    for i in 1:(num_v-2)
        for j in 1:num_v
            if degree[j] == 0 #If j is not present in prufer set 
                add_edge!(g, j, prufer[i]) #add edge
                degree[j] = -1 #Remove from Prufer
                degree[prufer[i]] -= 1 #low the degree of the vertex i
                break
            end
        end
    end

    # For the last element, we look at those vertices still with degree 0
    for i in 1:num_v
        if degree[i] == 0
            for j in (i+1):num_v
                if degree[j] == 0
                    add_edge!(g, i, j)
                    break
                end
            end
            break
        end
    end

    return g
end

function GraphToPrufer(g::SimpleGraph)::Vector{UInt8}
    
    local G = copy(g) #this is not strictly necessary, but it is useful because we are going to modify g
    local prufer = Vector{UInt8}(undef,(nv(G)-2))
        i::Int64 = 1
    #push!(prufer,convert(UInt8, nv(G))) #comment off this if we want that the first element of the sequence is the number of vertices
    
    while (nv(G)>2) #when G is a single edge, we stop
        for v in vertices(G) #we look for the smallest leaf
            local neighb_v = all_neighbors(G, v)
            if length(neighb_v) == 1 #if this is true, then we found the leaf with smallest label
                other_v = neighb_v[1] #the other vertex
                prufer[i] = convert(UInt8, other_v)
                i += 1
                rem_vertex!(G, v) #remove the vertex
                break #back to while
            end
        end
    end
    
    return sort(prufer) #I like to have my sequence sorted
end


"""
Extract a Prufer sequence and a number from a string. The number is the number of automorphisms of the graph with that Prufer sequence.
"""
function get_graph(str::String)::Tuple{SimpleGraph{Int64}, Int64}

    s = split(str, ',')
    g = PruferToGraph([parse(UInt8,s[1][i]) for i in 1:length(s[1])]) #get the graph
    a = parse(Int64,s[2])       
    return (g, a)
end

"""
Return all the arrays of length `l` such that the sum of all elements of the array is `d`.
"""
function get_weights(l::Int64, d::Int64)::Vector{Vector{Int64}}
    return vcat(unique.(permutations.(partitions(d,l)))...)
end