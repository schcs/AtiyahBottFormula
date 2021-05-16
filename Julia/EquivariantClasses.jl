using LightGraphs
include("Marked.jl")

"""
Equivariant class of curves meeting a linear subspace of codimension r.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `r::Int64`: the codimension of the subvariety. It can be either a number or a vector.
"""
function Incidency(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, r::Int64)::Rational{BigInt}
    
    local p1 = Rational{BigInt}(0);
    r -= 1
    
    col = Dict(vertices(g).=> coloration) #assing colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges

    for e in edges(g)
        for t in (0:r)
            p1 += d[e]*(scalars[col[src(e)]]^(t))*(scalars[col[dst(e)]]^(r-t))
        end
    end

    return p1

end

function Incidency(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, r::Vector{Int64})::Rational{BigInt}
    
    local p1 = Rational{BigInt}(1);
    
    for j in unique(r)
        p1 *= Incidency(g, coloration, weights, scalars, j)^count(x->x==j, r)
    end
    return p1

end


"""
Equivariant class of curves contained in a hypersurface of degree b.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `b::Int64`: the degrees of the hypersurface. It can be either a number or a vector.
"""
function Hypersurface(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, b::Int64)::Rational{BigInt}

    local p1=Rational{BigInt}(1)
    local q1=Rational{BigInt}(1)
    
    col = Dict(vertices(g).=> coloration) #assing colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    for e in edges(g)
        for alph in (0:b*d[e])
            p1 *= (alph*scalars[col[src(e)]]+(b*d[e]-alph)*scalars[col[dst(e)]])//d[e]
        end
    end
    
    for v in vertices(g)
        q1 *= (b*scalars[col[v]])^(1-length(all_neighbors(g, v)))   
    end

    return p1*q1

end

function Hypersurface(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, b::Vector{Int64})::Rational{BigInt}

    local p1 = Rational{BigInt}(1)
    
    for j in unique(b)
        p1 *= Hypersurface(g, coloration, weights, scalars, j)^count(x->x==j, b)
    end
    return p1
end


"""
Equivariant class of contact curves.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
"""
function Contact(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}})::Rational{BigInt}

    local p1 = Rational{BigInt}(1)
    local q1 = Rational{BigInt}(1)
    
    col = Dict(vertices(g).=> coloration) #assing colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    for e in edges(g)
        for alph in (1:2*d[e]-1)
            p1 *= (alph*scalars[col[src(e)]]+(2*d[e]-alph)*scalars[col[dst(e)]])//d[e]
        end
    end
    
    for v in vertices(g)
        q1 *= (2*scalars[col[v]])^(length(all_neighbors(g, v))-1)
    end

    return p1*q1

end


"""
Equivariant class of the pull-back of O(1) with respect to the i-th evaluation map.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `m::Rational{BigInt}`: the marks.
- `i::Rational{BigInt}`: the evaluation map.
"""
function O1_i(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, marks::marks, i::Int64)::Rational{BigInt}
    
    col = Dict(vertices(g).=> coloration) #assing colors to vertices
    
    return scalars[col[marks.get_vertex[i]]]
end


"""
Equivariant class of the pull-back of O(1) with respect to the product of all evaluation maps.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `m::Rational{BigInt}`: the marks.
"""
function O1(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, marks::marks)::Rational{BigInt}
    
    local p1 = Rational{BigInt}(1)
    col = Dict(vertices(g).=> coloration)
    
    for t in 1:marks.m
        p1 *= scalars[col[marks.get_vertex[t]]]
    end
    
    return p1
end


"""
The inverse of the equivariant class of the Euler class of the normal bundle.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `m::Rational{BigInt}`: the marks.
"""
function Euler_inv(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, marks::marks)::Rational{BigInt}
    local p1::Rational{BigInt} = one(Rational{BigInt})
    local q1::Rational{BigInt} = one(Rational{BigInt})
    local V::Rational{BigInt} = one(Rational{BigInt})
    local E::Rational{BigInt} = one(Rational{BigInt})
    
    col = Dict(vertices(g).=> coloration) #assing colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    omega_inv = Dict(edges(g).=> [d[e]//(scalars[col[src(e)]]-scalars[col[dst(e)]]) for e in edges(g)]) 
    merge!(omega_inv,Dict(reverse.(edges(g)).=> [d[e]//(scalars[col[dst(e)]]-scalars[col[src(e)]]) for e in edges(g)]))
    
    max_col = length(scalars)
    
    for e in edges(g)
        q1 = one(Rational{BigInt})
        for var_col in 1:max_col
            if var_col != col[src(e)] && var_col != col[dst(e)]
                    for alph in 0:d[e]
                        q1 *= ((alph*scalars[col[src(e)]]+(d[e]-alph)*scalars[col[dst(e)]])*(1//d[e])-scalars[var_col])
                    end
            end
        end
        E *= ((omega_inv[e])^(2*d[e]))*((-1)^d[e])//(factorial(d[e])^2)*(1//q1)
    end
    
    for v in vertices(g)
        p1 = one(Rational{BigInt})
        for j in 1:max_col
            if j != col[v]
                p1 *= scalars[col[v]]-scalars[j]
            end
        end
        p1 = p1^(length(all_neighbors(g, v))-1)

        s = sum([omega_inv[LightGraphs.SimpleEdge(v,w)] for w in all_neighbors(g, v)])^(length(all_neighbors(g, v))-3+num_marks(marks,v))
        pr = prod([omega_inv[LightGraphs.SimpleEdge(v,w)] for w in all_neighbors(g, v)])

        V *= p1*s*pr
    end

    return V*E
end