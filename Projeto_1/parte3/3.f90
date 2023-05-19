program proj1
    implicit real*8(a-h,o-z)
    real*8, dimension(0:9)::ener1,ener2
    write(*,*)"Escreva a quantidade de iterações"
    read(*,*)num
    write(*,*)"Digite o valor de omega"
    read(*,*)omega
    pi=acos(0.0)*2.0
    do n=0,9
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
        ener1(n)=ener
    enddo
    write(*,*)'metade' 
    do n=0,9
        flag=1
        ener=-1+0.01
        do while(flag==1)
            !calculando os extremos
                soma=0
                xo=-(-36*2**(5.0/6.0) - sqrt(72*2**(2.0/3.0)*(ener - 35) + 2592*2**(2.0/3.0)))/(36.0*2**(2.0/3.0))
                xi=-(-36*2**(5.0/6.0) + sqrt(72*2**(2.0/3.0)*(ener - 35) + 2592*2**(2.0/3.0)))/(36.0*2**(2.0/3.0))
                dex=(xo-xi)/(1.0*num)
                
                do j=0,num!trapezio
                    
                    soma=soma+sqrt(ener-func(xi+j*dex))*dex
                end do
                !vendo se a integral está no alcance do erro
                if(((n+0.5)*pi-omega*soma)<0.001)then
                    flag=0
                    do j=1,num
                        !saida de dados
                        write(10+n,*)xi+j*dex,func(xi+j*dex)
                    end do
                else
                    ener=ener+0.01
                endif
        end do
        ener2(n)=ener
    enddo
    do n=0,9
        write(*,'(f12.3,A,f12.3,A,1I5,A)')ener1(n),' &',ener2(n)," & ",n ,'\\'
        write(*,*)'\hline'
    enddo
end program

function func(x) result(f)
    real*8, intent (in) :: x ! input

    real*8              :: f ! output
    f = -1 + 18*2**(2.0/3.0)*(x - 2**(1.0/6.0))**2 
end function