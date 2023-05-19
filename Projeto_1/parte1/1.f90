program proj1
    implicit real*8(a-h,o-z)
    write(*,*)"Escreva a quantidade de iterações"
    read(*,*)num
    pi=acos(0.0)*2.0
    do n=0,10
        flag=1
        ener=0.1
        do while(flag==1)
        !calculando os extremos
            soma=0
            xo=sqrt(ener*2)
            xi=-sqrt(ener*2)
            dex=(xo-xi)/(1.0*num)
            !fazendo a integral trapezio
            do j=0,num
                soma=soma+sqrt(ener-(xi+j*dex)**2*0.5)*dex
            end do
            !vendo a condição de parada
            if(((n+0.5)*pi-sqrt(2.0)*soma)<0.001)then
                flag=0
                do j=1,num
                    !saida de dados, espaço de fase
                    write(10+n,*)xi+j*dex,sqrt(ener-(xi+j*dex)**2*0.5)
                end do
            else
                ener=ener+0.01
            endif
        end do
        write(*,*)ener,n
    enddo 
end program