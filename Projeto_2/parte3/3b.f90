program onda
    implicit real*8(a-h,o-z)
    real*8 ::L
    real*8, dimension(0:10000)::y
    !condições iniciais
    write(*,*)'escreva o tamanho da caixa  [-L,L]:'
    read(*,*)L
    write(*,*)'escreve a condição inicial y(0)'
    read(*,*)yz
    write(*,*)'escreve a condição inicial dy(0)/dx'
    read(*,*)yzl
    write(*,*)'escreve a energia da particula'
    read(*,*)E
    detx=(L*2.0)/10000
    pi=acos(-1.0)
    !metodo numerov
    !2*E=g(x)
    soma_err=0
    do i=1,9999
        !condições iniciais
        y(1)=yz
        y(2)=yzl*detx+y(1) 
        !passo seguinte  
        y(i+1)=(y(i)*(2.0-5.0*(detx*detx*2*E)/6.0)-y(i-1)*(1.0+(detx*detx*2*E)/12.0))/(1.0+(detx*detx*2*E)/12.0)    
        write(1,*)i*detx-L,y(i)/sqrt(2.0*L)
        !um erro simplificado
        soma_err=soma_err+abs(y(i)-cos(pi*i*detx/L))
    enddo
    write(*,*)"erro do metodo numerov",soma_err/10000!media do erro
end program