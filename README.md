# AtiyahBottFormula

This package contains an implementation of the Atiyah-Bott residue formula in the Julia language. The theory behind the package is described in the paper 
"Effective computations of the Atiyah-Bott formula" by Giosuè Muratore e Csaba Schneider (https://arxiv.org/pdf/2105.11183.pdf).

Run Julia on your computer from the Julia directory of the package. Once the Julia prompt appears, type 

julia> include( "Main.jl" )

You can check your current working directory inside Julia by typing 

julia> pwd() 

To use our code, the you should first define the equivariant classes to be calculated as "P = (g,c,w,s,m) ->...".
After the "->", one has to write an expression in the equivariant classes. After P is defined, one has to call the
Atiyah-Bott formula by the command AtiyahBottFormula(n,d,m,P). 

The full list of the currently supported equivariant classes is the following:

O1_i(g,c,w,s,m,i)<br>
O1(g,c,w,s,m)
Incidency(g,c,w,s,r)
Hypersurface(g,c,w,s,b)
Contact(g,c,w,s)
R1(g,c,w,s)

Brief descriptions on these functions can be obtained through the standard help functionality of Julia by typing "?" and then the name of the function.

For example, to compute the number of rational plane curves of degree d through 3d−1 general points, one may write

julia> d = 1 #for other values of d, change this line
julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2
julia> AtiyahBottFormula(2,d,3*d-1,P);


The number of degree d curves on a cubic surface passing through d-1 points:
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,3)*(Incidency(g,c,w,s,2)/3)^(d-1)
julia> AtiyahBottFormula(3,d,0,P);

The number of Degree d plane curves passing through 3d-1 points:
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^(3d-1)
julia> AtiyahBottFormula(2,d,0,P);

