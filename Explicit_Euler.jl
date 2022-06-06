using LinearAlgebra
using Plots
#_____________________________________________________________________________________________________________________
#Explicit Euler
#_____________________________________________________________________________________________________________________

function f(Un, J, p = 2.0, q = 0.05) # calcule le terme de gauche généré par la fonction f
    F = zeros(2*J,1)
    for j in 1:J
        F[j] = p - Un[j]*(Un[J+j]^2) # p - U1(n,j)* U2(n,j)^2
        F[J+j] = q - Un[J+j] + Un[j]* Un[J+j]^2
    end
    return F
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
    return Deriv2
end

function next_Un(Un, J, dx, dt, d) # calcule Un+1 en fonction de Un
    #U_next = Un + dt*( f(Un, J) + Deriv_spat(Un, J, dx, d) )
    U_next = Un + dt*Deriv_spat(Un, J, dx, d) 
    return U_next
end

function explicit_euler(f, u0, T, dt, d = 0.037)
    J = Int(size(u0)[1] /2)
    dx = (2.0*pi)/J
    N = Int(floor(T/dt))
    Un = u0
    Uevol = [Un]
    for i in 1:N
        Un = f(Un, J, dx, dt, d)
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


function animer_explicite(T, N, J, d=0.037)
    Intrvl = range(0.0, 2.0*pi, J)
    dt = T/N
    u0  = generer_U0(J)
    Usol = explicit_euler(next_Un, u0, T, dt, d)
    
    anim = @animate for t in 1:N
        plotter(Intrvl, Usol, t, J)
    end
    gif(anim,fps=30)
end

@time animer_explicite(10,4000,80, 1.0)



