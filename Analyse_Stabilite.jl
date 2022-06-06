#Analyse_Stabilité
using LinearAlgebra
using Plots

function maxi(Liste, Len)
    Lsol = []
    for i in 1:Len
        push!(Lsol, maximum(Liste[i]))
    end
    return(Lsol)
end


function Stabilite(n, Nb=100, inf=0, sup=0.1, p=2.0, q=0.05)
    u1 = p/( (p+q)^2)
    u2 = p+q 
    df = [ (-u2^2) (-2*u1*u2) ; (u2^2) (2*u1*u2 - 1)] # calcul de la jacobienne
    Val_d = range(inf,sup, Nb)
    Vp = []
    for d in Val_d
        D = [1. 0. ; 0. d]
        push!(Vp, real.(eigvals(df - (n^2)*D)) )
    end
    plot(Val_d, maxi(Vp, Nb), label="n="*string(n))
end

# Stabilite(1)
# Stabilite(2)
# Stabilite(3)
Stabilite(4)
# Stabilite(5)