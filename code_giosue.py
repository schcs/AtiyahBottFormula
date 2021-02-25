load("moduli.py")
import time
start = time.time()
n = 4 #the dimension of the projective space
d = 3 #the degree of the curves
M = ModuliStableCurve(n,d)
result1 = 0
for g in M.fixed_locus():
    result1 += ((1/g.altern())*(g.quintic())*(1/g.c_top())).subs(x0=0,x1=1,x2=5,x3=23,x4=31)
end = time.time()
print("Degree ", d)
print("Result ", result1)
print("Elapsed time in seconds ", (end-start))