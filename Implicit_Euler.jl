#Implicit_Euler

using LinearAlgebra
using Plots

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

function f(Un, p = 2.0, q = 0.05) # calcule le terme de gauche généré par la fonction f
    J = Int(size(Un)[1]/2)
    F = zeros(2*J,1)
    for j in 1:J
        F[j] = p - Un[j]*(Un[J+j]^2) # p - U1(n,j)* U2(n,j)^2
        F[J+j] = q - Un[J+j] + Un[j]* Un[J+j]^2
    end
    return F
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

function implicit_euler(fct, u0, T, dt, d = 0.037)
    N = Int(floor(T/dt))
    Un = u0
    Uevol = [Un]
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


function animer_implicite(T, N, J, d=0.037)
    Intrvl = range(0.0, 2.0*pi, J)
    dt = T/N
    u0  = generer_U0(J)
    Usol = implicit_euler(next_Un_implicit, u0, T, dt, d)
    
    anim = @animate for t in 1:N
        plotter(Intrvl, Usol, t, J)
    end
    gif(anim,fps=30)
end

@time animer_implicite(10,4000,80, 1.0)