#implicite

#Crée la jacobienne de F() en U
function DF(U,J)
    res=zeros(2*J,2*J)
    for i in 1:Nx
        res[i,i]=-U[Nx+i]^2
        res[Nx+i,i]= -2U[Nx+i]*U[i]
        res[i,Nx+i]= U[Nx+i]^2 
        res[Nx+i,Nx+i]=-1+2U[Nx+i]*U[i] 
    end
    res
end

#Crée la matrice DX2 tel que DX2*U= D*d²U/dx²
function D(Nx,d)
    res=zeros(2*Nx,2*Nx)

    for i in 1:Nx
        if i<Nx
            res[i,i+1]=1.
            res[i+1,i]=1.
        end
        res[i,i]=-2.
    end

    for i in Nx+1:2Nx
        if i<2Nx
            res[i,i+1]=d
            res[i+1,i]=d
        end
        res[i,i]=-2. * d
    end

    for K in 0:1
        res[K*Nx+1,K*Nx+Nx]=1. * d^K
        res[K*Nx+Nx,K*Nx+1]=1. * d^K
    end
    dx=2π/(Nx-1)
    res/=dx^2
end

#Calcule le terme suivant en utilisant la méthode de Newton avec
#une seule iterations
function implicite(p,q,d,dt,U,Nx,DX2) 
    M=I-dt*DX2
    #G=dt*F(p,q,U,Nx) -M*U
    G=dt*(F(p,q,U,Nx) +DX2*U)
    df=DF(U,Nx)
    res=(-dt*df +M)\ G
    return U+res
end

#Calcule le terme suivant en utilisant la méthode de Newton avec
# 1 à 3 iterations
function implicite1(p,q,d,dt,U,Nx,DX2)
    tol=1e-1
    nU=U
    dU=zeros(2Nx,1)
    ite=0
    while (norm(dt*(F(p,q,nU,Nx) +DX2*nU)-(nU-U))>tol && ite<=3) || ite<=1
        G=dt*(F(p,q,nU,Nx) +DX2*nU)
        df=DF(nU,Nx)
        dU=(I-dt*(df+DX2 ))\ G
        nU+=dU
        ite+=1
    end
    return nU
end