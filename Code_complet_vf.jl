using LinearAlgebra
using Plots
using DifferentialEquations

#_____________________________________________________________________________________________________________________
# Analyse_Stabilité
#_____________________________________________________________________________________________________________________

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
    for i in 1:n
        plot!(Val_d, maxi(Vp, Nb), label="n="*string(n))
    end
end
# Stabilite(5)

#_____________________________________________________________________________________________________________________
# Fonction f (caractéristique du problème)
#_____________________________________________________________________________________________________________________

function f(Un, p = 2.0, q = 0.05) # calcule le terme de gauche généré par la fonction f
    J = Int(size(Un)[1]/2)
    F = zeros(2*J,1)
    for j in 1:J
        F[j] = p - Un[j]*(Un[J+j]^2) # p - U1(n,j)* U2(n,j)^2
        F[J+j] = q - Un[J+j] + Un[j]* Un[J+j]^2
    end
    return F
end

#_____________________________________________________________________________________________________________________
# Schéma Explicite
#_____________________________________________________________________________________________________________________

function Deriv_spat(U, d) #calcule le terme de droite généré par le produit de la matrice D et de la dérivée seconde de U en x
    J = Int(size(U)[1]/2)
    dx = (2*pi)/(J-1)
    dx2 = dx^2

    Deriv2 = zeros(2*J,1)
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
    return Deriv2
end

function next_Un_explicit(Un, dt, d) # calcule Un+1 en fonction de Un
    U_next = Un + dt*( f(Un) + Deriv_spat(Un, d) )
    return U_next
end

#_____________________________________________________________________________________________________________________
# Schéma Implicite
#_____________________________________________________________________________________________________________________

function Jacobienne(U,J) # calcule la matrice de la Jacobienne de f pour chaque pas de temps d'un vecteur U donné
    Jcb=zeros(2*J,2*J)
    for i in 1:J
        Jcb[i,i]= -U[J+i]^2
        Jcb[J+i,i]= -2*U[J+i]*U[i]
        Jcb[i,J+i]= U[J+i]^2 
        Jcb[J+i,J+i]= 2*U[J+i]*U[i] - 1
    end

    return Jcb
end

function MatDeriv(J,d)
    d2D=zeros(2*J,2*J)

    for i in 1:J # diagonale
        d2D[i,i] = -2.0
        d2D[J+i,J+i] = -2.0*d
    end
    for i in 1:(J-1) 
        # surdiagonale
        d2D[i,i+1] = 1.0
        d2D[J+i,J+i+1] = d
        # sousdiagonale
        d2D[i+1,i] = 1.0
        d2D[J+i+1,J+i] = d
    end
    # termes de coins
    d2D[1,J] = 1.0
    d2D[J,1] = 1.0
    d2D[J+1, 2*J] = d
    d2D[2*J, J+1] = d
    # on divise par dx^2
    dx = (2*pi)/(J-1)
    d2D = d2D/(dx^2)

    return d2D
end

function newton(U, dt, DX2, max_ite = 4, tol = 0.1)
    J = Int(size(U)[1]/2)
    Uit = U
    dU=zeros(2*J,1)
    i = 0
    while (i == 0) || ( ( norm( dt*( f(Uit) + DX2*Uit ) + U - Uit ) > tol) && (i < max_ite) )
        G = dt * (f(Uit) + DX2*Uit)
        jcbn = Jacobienne(Uit, J)
        dU = (I - dt*(jcbn+DX2))\ G
        Uit = Uit + dU
        i+=1
    end
    return Uit
end

function next_Un_implicit(Un, dt, d)
    J = Int(size(Un)[1]/2)
    MatDrv = MatDeriv(J,d)
    Unext = newton(Un, dt, MatDrv)

    return Unext
end

#_____________________________________________________________________________________________________________________
# Fonction liant le tout
#_____________________________________________________________________________________________________________________

function schema_euler(type, u0, T, dt, d = 0.037)
    N = Int(floor(T/dt))
    Un = u0
    Uevol = [Un]
    if type == "explicit"
        fct = next_Un_explicit
    else
        fct = next_Un_implicit
    end
    for i in 1:N
        Un = fct(Un, dt, d)
        push!(Uevol, Un)
    end
    return Uevol
end

#_____________________________________________________________________________________________________________________
# Fonctions supplémentaires pour afficher facilement
#_____________________________________________________________________________________________________________________

function generer_U0(J, p=2.0, q=0.05) # cosinus + pt d'équilibre
    U0=zeros(2*J,1)
    X=range(0.0, 2.0*pi, J)
    U0[1:J] = p/(p+q)^2 .+ 0.2*cos.(X)
    U0[(J+1):(2*J)] = p+q .+ 0.2*cos.(X)
    return U0
end

function plotter(Interval, U, t, J)
    plot(Interval, [U[t][1:J] U[t][(J+1):2*J] ]  )
end

function animer_schema_euler(T, N, J, type = "explicit", d=0.037)
    Intrvl = range(0.0, 2.0*pi, J)
    dt = T/N
    u0  = generer_U0(J)
    Usol = schema_euler(type, u0, T, dt, d)

    anim = @animate for t in 1:N
        plotter(Intrvl, Usol, t, J)
    end
    gif(anim,fps=30)
end

#animer_schema_euler(10,4000,80, "implicit")
#animer_schema_euler(10,4000,80, "explicit")

#_____________________________________________________________________________________________________________________
# Black Box ODE
#_____________________________________________________________________________________________________________________

function fctn_evolution(dU,U,prmtr,t)
    (J,d) = prmtr
    Unxt = f(U) + Deriv_spat(U,d)
    for i in 1:(2*J)
        dU[i]= Unxt[i]
    end
end

function BlackBox(U0,T,d=0.037)
    J = Int(size(U0)[1]/2)
    prmtr = (J,d)
    time_span = (0., T)
    pblm = ODEProblem(fctn_evolution, U0, time_span, prmtr)
    sol = solve(pblm)
    return sol
end

#BlackBox(generer_U0(80), 10)
