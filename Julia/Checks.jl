import Base: *, //, /, ^, +, -, inv

struct Cycle #We define the structure Cycle. It keeps track of the codimension and of the fact that it contains psi-classes
    c::Int64 #The codimension. If negative, the class is not well defined
    is_psi::Int64 #this is a counter that measure how many times the class Psi is called in P. If this number is strictly greater than 1, then the program stops.
end

#Cycles satisfy the following arithmetic. They are not invertible (except in codimension 0), and cycles of different codimension are not summable.
#So the Cycle(-1, 0) is a meaningless cycle.

*(P1::Cycle, P2::Cycle)::Cycle = (P1.c < 0 || P2.c < 0 ) ? Cycle(-1, 0) : Cycle(P1.c + P2.c, P1.is_psi + P2.is_psi)
^(P1::Cycle, n::Int64) ::Cycle = Cycle(n*P1.c, n*P1.is_psi) #if n=0, we got a zero cycle
+(P1::Cycle, P2::Cycle)::Cycle = P1.c == P2.c ? Cycle(P1.c, max(P1.is_psi, P2.is_psi)) : Cycle(-1, 0)
-(P1::Cycle, P2::Cycle)::Cycle = +(P1, P2)
inv(P1::Cycle)         ::Cycle = P1.c == 0 ? Cycle(0, 0) : Cycle(-1, 0)
*(P1::Cycle, ::Number) ::Cycle = P1
*(::Number, P1::Cycle) ::Cycle = P1
//(P1::Cycle, ::Number)::Cycle = P1
//(::Number, P1::Cycle)::Cycle = P1
/(P1::Cycle, ::Number) ::Cycle = P1
/(::Number, P1::Cycle) ::Cycle = P1

"""
    dim_M(n, d, m)

The dimension of the moduli space of stable rational map to the projective space of dimension `n`, of degree `d` with `m` marks.
# Arguments
- `n::Int64`: the dimension of the projective space.
- `d::Int64`: the degree of the stable maps.
- `m::Int64`: the number of marks.

# Example
```jldoctest
julia> dim_M(2,2,5)
10
```
"""
function dim_M(n::Int64, deg::Int64, n_marks::Int64)::Int64
    
    return n + (n + 1)*deg + n_marks - 3
end


"""
    codim(n, d, m, P)

The codimension of the equivariant class `P`.
# Arguments
- `n::Int64`: the dimension of the projective space.
- `d::Int64`: the degree of the stable maps.
- `m::Int64`: the number of marks.
- `P`: the equivariant class.

# Example
```jldoctest
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);
julia> codim(4,1,0,P)
6
```
"""
function codim(n::Int64, deg::Int64, n_marks::Int64, P)::Int64
    
    return P(n, deg, n_marks, 0, 0).c
end

"""
    is_zero_cycle(n, d, m, P)

Return `true` if the equivariant class `P` is a 0-cycle in the moduli space, `false` otherwise.
# Arguments
- `n::Int64`: the dimension of the projective space.
- `deg::Int64`: the degree of the stable maps.
- `n_marks::Int64`: the number of marks.
- `P`: the equivariant class.

# Example
```jldoctest
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);
julia> is_zero_cycle(4,1,0,P)
true
```
"""
function is_zero_cycle(n::Int64, deg::Int64, n_marks::Int64, P)::Bool
    
    if !isa(P, Array)  #we want that P is an array, possibly with only one element
        P = [P]
    end

    for res in 1:length(P)
        local P_cycle = P[res](n, deg, n_marks, 0, 0)

        if P_cycle.c < 0
            return false
        end

        if P_cycle.is_psi > 1
            printstyled("Warning: ", bold=true, color=:light_yellow)
            println("more instances of Psi has been found. Type:")
            printstyled("julia> ", bold=true, color=:light_green)
            println("?Psi")
            println("for support.")
            return false
        end

        if P_cycle.c != dim_M(n, deg, n_marks)
            printstyled("Warning: ", bold=true, color=:light_yellow)
            length(P)==1 ? println("the class is not a zero cycle") : println("some classes are not zero cycles")
            return false
        end
    end
    return true
end

"""
    check_Data()

List of all files containing the colorations in the folder Data.
"""
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
