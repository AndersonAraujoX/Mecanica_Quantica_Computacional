program proj1
    implicit real*8(a-h,o-z)
    write(*,*)"Escreva a quantidade de iterações"
    read(*,*)num
    write(*,*)"Digite o valor de omega"
    read(*,*)omega
    pi=acos(0.0)*2.0
    do n=0,10
        flag=1
        ener=-1+0.0001
        do while(flag==1 .and. ener<=0.0)
            !calculando os extremos
            soma=0
            xo=((-2.0-2.0*sqrt(1.0+ener))/ener)**(1.0/6.0)
            xi=((-2.0+2.0*sqrt(1.0+ener))/ener)**(1.0/6.0)
            dex=(xo-xi)/(1.0*num)
            do j=1,num-1!trapezio
                soma=soma+sqrt(ener-4*((xi+j*dex)**(-12)-(xi+j*dex)**(-6)))*dex
            end do
            !vendo se a integral está no alcance do erro
            if(((n+0.5)*pi-omega*soma)<0.0001)then
                flag=0
                do j=0,num!saida de dados
                    write(10+n,*)xi+j*dex,sqrt(ener-4*((xi+j*dex)**(-12)-(xi+j*dex)**(-6)))
                end do
            else
                ener=ener+0.01
            endif
        end do
        if (ener<=0.0)then
            write(*,*)ener,n
        endif
    enddo 
end program