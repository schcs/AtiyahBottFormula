"""
The structure `marks` define the marked vertex of a graph.
# Arguments
- `n_vert::Int64`: the number of vertices.
- `m::Int64`: the number of marks.
- `get_vertex::Dict{Int64,Int64}`: a dictionary Dict{Int64,Int64} which associate to each mark a vertex.
"""
#if !@isdefined marks 
    mutable struct marks
        n_vert::Int64
        m::Int64
        get_vertex::Dict{Int64,Int64}
    end
#end


# function marks(g::SimpleGraph,m::Int64)::marks
#     return marks( nv(g), m , Dict{Int64, Int64}() )
# end

####???????????????????
function empty_mark()::marks
    return marks( 0, 0 , Dict{Int64, Int64}() )
end
# empty_mark = marks( 0, 0 , Dict{Int64, Int64}() )
####???????????????????

function marks(v::Int64,m::Int64)::marks
    return marks( v, m , Dict{Int64, Int64}() )
end

function Base.iterate( marks::marks, c=0)
    
    if marks.m == 0 && c == 1 
        return nothing
    end

    if marks.get_vertex == Dict{Int64, Int64}()
        marks.get_vertex = Dict([i for i in 1:marks.m].=>[1 for i in 1:marks.m])
        return marks, 1
    end

    index = findfirst(x -> x<marks.n_vert, [marks.get_vertex[x] for x in 1:marks.m])
    if index == nothing
        return nothing
    end

    next_mark!(marks, index)
    return marks, 1
end


function next_mark!(marks::marks, index::Int64)
    
    marks.get_vertex[index] += 1
    for i in 1:(index-1)
        marks.get_vertex[i] = 1
    end
    
    return nothing
end            

"""
Take the marks of a graph and return a dictionary Dict{Int64,Vector{Int64}} which associate to each vertex the set of marks on the vertex.
"""
function invert_marks(mark::marks)::Dict{Int64,Vector{Int64}}
    new_D=Dict{Int64,Vector{Int64}}()
    D = mark.get_vertex
    for k in keys(D)
         if D[k] in keys(new_D)
           push!(new_D[D[k]],k)
         else
           new_D[D[k]] = [k]
         end
    end
    for v in 1:mark.n_vert
        if !(v in keys(new_D))
            new_D[v]=Int64[]
        end
    end
    return new_D
end

"""
Return the number of marks assigned to some vertex v.
# Arguments
- `mark::marks`: the marks assigned to a graph.
- `v::Int64`: the vertex.
"""
function num_marks(mark::marks, v::Int64)::Int64
    return sum([mark.get_vertex[w] == v for w in 1:mark.m])
end