# Code complet

# Romain Gérard
# Nicolas Sereyjol-Garros

using LinearAlgebra
using Plots


#_____________________________________________________________________________________________________________________
# 1) Analyse de stabilité
#________________________
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
# Stabilite(4)
# Stabilite(5)

#_____________________________________________________________________________________________________________________
# 2) Euleur explicite
#____________________
# Principe: On a un intervalle [inf, sup] = [0, 2pi] découpé en J points expacés de dx, pour une durée T de N espaces temporels espacés de dt
# Pour n <= N, on stocke les deux composantes de U(t=n*dt) dans un vecteur Un tq pour tout j<=J, U(t=n*dt, x=j*dx)_1 = Un[j] et  U(t=ndt)_2 = Un[j+J] 

function f(Un, J, p, q) # calcule le terme de gauche généré par la fonction f
    F = zeros(2*J,1)
    for j in 1:J
        F[j] = p - Un[j]*(Un[J+j]^2) # p - U1(n,j)* U2(n,j)^2
        F[J+j] = q - Un[J+j] + Un[j]* Un[J+j]^2
    end
    return(F)
end

function Deriv_spat(U, J, dx, d) #calcule le terme de droite généré par le produit de la matrice D et de la dérivée seconde de U en x
    Deriv2 = zeros(2*J,1)
    dx2 = dx^2
    # Conditions périodiques pour U_1
    Deriv2[1] = (U[2] - 2*U[1] + U[J]) /dx2
    Deriv2[J] = (U[1] - 2*U[J] + U[J-1]) /dx2
    # Conditions périodiques pour U_2
    Deriv2[J+1] = d* (U[J+2] - 2*U[J+1] + U[2*J])/dx2
    Deriv2[2J] = d* (U[J+1] - 2*U[2*J] + U[2*J-1])/dx2
    # Pas de pb pour les autres indices
    for j in 2:(J-1)
        Deriv2[j ] = (U[j+1] - 2* U[j] + U[j-1])/dx2
        Deriv2[j + J] = (U[(j+1)+ J] - 2* U[j+ J] + U[(j-1)+ J])*(d/dx2)
    end
    return(Deriv2)
end

function next_Un(Un, dt, J, dx, p, q, d) # calcule Un+1 en fonction de Un
    U_next = Un + dt*( f(Un, J, p, q) + Deriv_spat(Un, J, dx, d) )
    return(U_next)
end


function generer_initial(inf, sup, J, p, q) # cosinus
    U0=zeros(2*J,1)
    X=range(inf, sup, J)
    U0[1:J]= cos.(X)
    U0[(J+1):(2*J)]= cos.(X)
    U0
end

function generer_U0(inf, sup, J, p, q) # cosinus + pt d'équilibre
    U0=zeros(2*J,1)
    X=range(inf, sup, J)
    U0[1:J]=p/(p+q)^2 .+ cos.(X)
    U0[(J+1):(2*J)]=p+q .+ cos.(X)
    U0
end


function Euleur_explicite( # Cette fonction retourne la liste des Un pour chaque pas de temps n dans [1:N]
    T, N, # durée, nb de pas de temps, écart temporel
    J, # Nb de pas spatiaux
    d, # Paramètre de la matrice D
    inf=0., sup=(2*pi), # intervalle = [inf, sup] = [0, 2pi]
    p=2., q=0.05) # (p,q)=(2. , 0.05)

    dx = (sup-inf)/(J-1)

    dt = T/N

    Un = generer_U0(inf, sup, J, p, q)
    Utot = [Un] # condition initiale
    for n in 1:N
        Un = next_Un(Un, dt, J, dx, p, q, d)
        push!(Utot, Un)
    end
    return(Utot)
end

Eul_ex = Euleur_explicite(25, 5000, 50, 0.0037)

function plotter(Interval, U, t, J)
    plot(Interval, [U[t][1:J] U[t][(J+1):2*J] ]  )
end

function animer_explicite(T, N, J, d, # T: temps total, N: nb de pas de temps, J: nb de pas spatiaux, d: cst de la matrice D
    inf=0.0, sup=(2*pi), # intervalle spatial (ici [0,2pi])
    p=2.0,q=0.05) # cstes du problème

    Intrvl = range(inf, sup, J)
    Usol = Euleur_explicite(T, N, J, d)
    
    anim = @animate for t in 1:N
        plotter(Intrvl, LU, t, J)
    end
    gif(anim,fps=30)
end

# animer_explicite(10,1000,40,0.003)

#_____________________________________________________________________________________________________________________
# 3) Implicite Euler
#____________________
# Note: on a utilisé ici un schéma implicte pour la dérivée seconde mais on a gardé les termes en fonction de Un dans la partie du f pour pouvoir utiliser la méthode de Richardson

Time = 0.5
dt = 0.0009
dx = 0.09
t0 = 0.1
It = Int(floor(Time/dt))
dt = Time/It
Ix = Int(floor(1.0/dx))
dx = 1.0/Ix

u=(a=[0.0 for x in 1:Ix];[a for t in 1:It])
v=(a=[0.0 for x in 1:Ix];[a for t in 1:It])
p=2.0
q=0.05
d1=1
d2=0.1
alpha=1



function matrix_D(d) #matrice de diffusion
    B = zeros(Int(Ix),Int(Ix))
    # could use diagflat to do this in a less pedestrian way
    for i in 1:Ix
        B[i,i] = 2
    end
    for i in 2:Ix
        B[i,i-1] = -1
    end
    for i in 1:Ix-1
        B[i,i+1] = -1
    end 
    B[1,Ix] = -1
    B[Ix,1] = -1
    return (I+(dt/dx^2)*d*B)
end

function vector_b(n,u,v,p,q;equation)
    b=[0. for j in 1:Ix]
    if equation==:U
        for j in 1:Ix
            b[j]=u[n][j]*(1-dt*v[n][j]^2)-dt*p
        end
    elseif equation==:V
        for j in 1:Ix
            b[j]=u[n][j]*(1+dt*v[n][j]^2)-dt*q-v[n][j]
        end
    end
    return b
end

function richardson_it(A,b,x,alpha)
    while norm(A*x-b)>10^-8
        x=x-alpha*(A*x-b)
    end
    return x
end

function eul_imp(u,v,p,q,d1,d2,alpha)
    A1=matrix_D(d1)
    A2=matrix_D(d2)
    for n in 1:It-1
        b1=vector_b(n,u,v,p,q;equation=:U)
        u[n+1]=richardson_it(A1,b1,u[n],alpha)
        b2=vector_b(n,u,v,p,q;equation=:V)
        v[n+1]=richardson_it(A2,b2,v[n],alpha)
    end
    return u,v
end

eul_imp(u,v,p,q,d1,d2,alpha)

#_____________________________________________________________________________________________________________________
# 4) Black box ODE
#____________________

using DifferentialEquations

function fonc(du,Un,paramtr,t)
    (J, d, p, q) = paramtr
    dx = 2*pi/(J-1)
    res=f(Un, J, p, q) + Deriv_spat(Un, J, dx, d)
    for i in 1:(2*J)
        du[i]=res[i]
    end
end

function black_box_ode(T, J, d, p, q)
    paramtr = (J, d, p, q)
    U0 = generer_U0(0., 2*pi, J, p, q)
    prob = ODEProblem(fonc,U0,(0.,T),paramtr)
    sol = solve(prob)
end

# black_box_ode(10, 50, 0.037, 2., 0.05)