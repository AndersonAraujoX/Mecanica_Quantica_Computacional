program onda
implicit real*8(a-h,o-z)
!condições iniciais
write(*,*)'escreva xl:'
read(*,*)xl
write(*,*)'escreva xr:'
read(*,*)xr
detx=(xr-xl)/10000.0
do i=0,10000
    !resetando a normalização
    write(2,*)i*detx+xl,g(xl+i*detx)
enddo
end program
function g(x)
    real*8, intent (in) :: x
    real*8             :: g
    g = 4*(x**(-12.0)-x**(-6.0))
end function