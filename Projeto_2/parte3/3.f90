program onda
    implicit real*8(a-h,o-z)
    real*8 ::L
    real*8, dimension(0:10000)::y
    !condições iniciais
    write(*,*)'escreva o tamanho da caixa indo [-L,L]:'
    read(*,*)L
    write(*,*)'escreve a condição inicial y(0)'
    read(*,*)yz
    write(*,*)'escreve a condição inicial dy(0)/dx'
    read(*,*)yzl
    write(*,*)'escreve a energia da particula'
    read(*,*)E
    !metodo da discretização
    detx=(L*2.0)/10000
    pi=acos(-1.0)
    soma_err=0.0
    do i=1,9999
        !condições iniciais
        y(0)=yz
        y(1)=yzl*detx+y(0)   
        !passo seguinte  
        y(i+1)=2*y(i)-y(i-1)-2*(detx*detx)*E*y(i)   
        write(1,*)i*detx-L,y(i)/sqrt(2.0*L)
        !um erro simplificado
        soma_err=soma_err+abs(y(i)-cos(pi*i*detx/L))
    enddo
    write(*,*)"erro do metodo numerov",soma_err/10000!media do erro
end program