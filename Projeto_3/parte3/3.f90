program onda
implicit real*8(a-h,o-z)
real*8 ::norm1,norm2
integer::xm
real*8, dimension(0:10000)::yl,yr
!condições iniciais
write(*,*)'escreva xl:'
read(*,*)xl
write(*,*)'escreva xr:'
read(*,*)xr
write(*,*)'escreva o ponto de encontro:'
read(*,*)enc
write(*,*)'escreva a energia do lennard-jones:'
read(*,*)eps
write(*,*)'escreva o sigma do potencial:'
read(*,*)sig
detx=(xr-xl)/10000.0
yz=0
yzl=0.001
j=0
E=-eps
do while(E<0.0)
    !resetando a normalização
    E=E+0.01
    norm2=0
    norm1=0
    do i=1,9999
        !condições iniciais para o lado esquerdo
        yl(0)=yz
        yl(1)=yzl*detx+yl(0) 
        !passo seguinte 
        yl(i+1)=2*yl(i)-yl(i-1)-2*detx*detx*g((i+1)*detx+xl,E,eps,sig)*yl(i) 
        !condições iniciais para direita
        yr(10000)=yz
        yr(9999)=yzl*detx+yr(10000)     
        !sentido contrário
        yr(10000-(i+1))=2*yr(10000-(i))-yr(10000-(i-1))-2*detx*detx*g(xr-i*detx,E,eps,sig)*yr(10000-(i)) 
    enddo
    !ponto de encontro
    xm=int((enc-xl)*10000/(xr-xl))
    !derivada
    const=yr(xm)/yl(xm)
    derl=const*(yl(xm+1)-yl(xm-1))/(2*detx)
    derr=(yr(xm+1)-yr(xm-1))/(2*detx)
    !write(*,*)derr/derl,E 
    !plotar as seguintes energias   
    if((1+0.01)>=derr/derl .and. derr/derl>=(1-0.01) )then
        j=j+1
        write(*,*)"energias permitidas"
        write(*,*)E
        E=E+1 
        !normalização
        do i=0,xm
            norm2=const*const*yl(i)*yl(i)*detx+norm2
        enddo
        do i=xm,10000
            norm1=yr(i)*yr(i)*detx+norm1
        enddo
        do i=0,xm
            write(10+j,*)i*detx+xl,const*yl(i)/sqrt(norm1+norm2)
        enddo
        do i=xm,10000
            write(10+j,*)i*detx+xl,yr(i)/sqrt(norm1+norm2)
        enddo
    endif
enddo
end program
function g(x,E,eps,sig)
    real*8, intent (in) :: x,E,eps,sig
    real*8             :: g
    g = (E-4*eps*((sig/x)**(12.0)-(sig/x)**(6.0))) 
end function