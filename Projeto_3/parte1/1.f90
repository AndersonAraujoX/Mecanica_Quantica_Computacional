program onda
implicit real*8(a-h,o-z)
real*8 ::norm1,norm2,k,L
integer::xm
real*8, dimension(0:1000000)::yl,yr
!condições iniciais
write(*,*)'escreva  da distância que vamos começar a onda L e -L:'
read(*,*)L
write(*,*)'escreva a constante da mola:'
read(*,*)k
N=10000
detx=(2.0*L)/N
yz=0
yzl=0.0001
j=0
E=0
write(*,*)"energias permitidas"
do while(j<10)
    !resetando a normalização
    E=E+0.0001
    norm2=0
    norm1=0
    do i=1,N-1
        !condições iniciais para o lado esquerdo
        yl(0)=yz
        yl(1)=yzl*detx+yl(0) 
        !passo seguinte 
        yl(i+1)=2*yl(i)-yl(i-1)-2*detx*detx*g(i*detx-L,E,k)*yl(i) 
        !condições iniciais para direita
        yr(N)=yz
        yr(N-1)=yzl*detx+yr(N)     
        !sentido contrário
        yr(N-(i+1))=2*yr(N-(i))-yr(N-(i-1))-2*detx*detx*g(L-i*detx,E,k)*yr(N-(i)) 
    enddo
    !ponto de encontro
    xm=5700
    !derivada
    const=yr(xm)/yl(xm)
    derl=(yl(xm+1)-yl(xm-1))/(2*detx)
    derr=(yr(xm+1)-yr(xm-1))/(2*detx*const)
    !plotar as seguintes energias
    !write(*,*)derr,derl,abs(derr-derl),E 
    if((1+0.001)>=derr/derl .and. derr/derl>=(1-0.001))then
        
        write(*,*)E,'&',(j+0.5)*sqrt(k)
        write(1,*)E,j
        j=j+1
        E=E+0.5  
        !normalização
        do i=0,xm
            norm2=const*const*yl(i)*yl(i)*detx+norm2
        enddo
        do i=xm,N
            norm1=yr(i)*yr(i)*detx+norm1
        enddo
        do i=0,xm
            write(10+j,*)i*detx-L,const*yl(i)/sqrt(norm1+norm2)
        enddo
        do i=xm,N
            write(10+j,*)i*detx-L,yr(i)/sqrt(norm1+norm2)
        enddo
    endif
enddo
end program
function g(x,E,k)
    real*8, intent (in) :: x,E,k
    real*8 :: g
    g = (E-k*x*x/2.0) 
end function