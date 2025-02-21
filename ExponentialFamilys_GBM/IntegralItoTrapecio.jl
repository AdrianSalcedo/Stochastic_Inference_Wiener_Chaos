#integral de ito por trapecio acumulado
function Int_Trap_Acum_2(senos,cosenos,j,t,W)
    npoints=length(senos)
    integral= zeros(npoints)

    delta=t[2]-t[1]
    TF=last(t)
    for i in 2:npoints
        integral[i]=(senos[i]*W[i]-senos[i-1]*W[i-1])-(delta*j*pi)*(W[i-1]*cosenos[i-1]+W[i-1]*cosenos[i-1])/(2*TF)
    end

    integrals=sum(integral)
    return(integrals)
end
