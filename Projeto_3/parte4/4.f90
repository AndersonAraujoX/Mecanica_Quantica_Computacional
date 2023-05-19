program onda
    implicit real*8(a-h,o-z)
    real*8 ::norm1,norm2,l
    integer::xm
    real*8, dimension(0:10000)::yl,yr
    !condições iniciais
    write(*,*)'escreva xr:'
    read(*,*)xr
    write(*,*)'escreva o ponto de encontro:'
    read(*,*)enc
    write(*,*)'escreva a nivel de energia:'
    read(*,*)l
    xl=0
    detx=(xr-xl)/10000.0
    yz=0
    yzl=0.0001
    j=0
    E=-0.6
    do i=0,10000
        write(30,*)(i*detx),-1/(i*detx)
    enddo
    do while(j<3)
        !resetando a normalização
        E=E+0.000001
        norm2=0
        norm1=0
        do i=1,9999
            !condições iniciais para o lado esquerdo
            yl(0)=yz
            yl(1)=yzl*detx+yl(0) 
            !passo seguinte
            yl(i+1)=2*yl(i)-yl(i-1)+detx*detx*g(i*detx,E,l)*yl(i) 
            !condições iniciais para direita
            yr(10000)=yz
            yr(9999)=yzl*detx+yr(10000)     
            !sentido contrário
            yr(10000-(i+1))=2*yr(10000-(i))-yr(10000-(i-1))+detx*detx*g(xr-i*detx,E,l)*yr(10000-i) 
        enddo
        !ponto de encontro
        xm=int((enc-xl)*10000/(xr-xl))
        !derivada
        const=yr(xm)/yl(xm)
        derl=const*(yl(xm+1)-yl(xm-1))/(2*detx)
        derr=(yr(xm+1)-yr(xm-1))/(2*detx)
        !write(*,*)derr,derl,derr,E 
        !plotar as seguintes energias   
        if((1+0.01)>=derr/derl .and. derr/derl>=(1-0.01))then!or. abs(E+0.055)<0.001
            j=j+1
            write(*,*)"energias permitidas"
            write(*,*)E
            E=E+0.001 
            !normalização
            do i=0,xm
                norm2=const*const*yl(i)*yl(i)*detx+norm2
            enddo
            do i=xm,10000
                norm1=yr(i)*yr(i)*detx+norm1
            enddo
            do i=0,xm
                write(10+j,*)i*detx+xl,const*yl(i)/(sqrt(norm1+norm2)*(i*detx))
            enddo
            do i=xm,10000
                write(10+j,*)i*detx+xl,yr(i)/(sqrt(norm1+norm2)*(i*detx))
            enddo
        endif
    enddo
    end program
    function g(r,E,l)
        real*8, intent (in) :: r,E,l
        real*8             :: g
        g = (2.0*(-1.0/r-E)+l*(l+1.0)/(r*r)) 
    end function