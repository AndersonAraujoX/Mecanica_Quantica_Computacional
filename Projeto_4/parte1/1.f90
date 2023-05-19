program onda
  implicit real*8(a-h,o-z)
  real*8, dimension(0:100000)::El
  !metodo da discretização
  write(*,*)"escreva delta"
  read(*,*)delta
  x=0
  i=0
  N=100000
  total=0
  !temos que energia local é igual 2 na parte do numerador
  do while(i<=N)
    !razão entre as probabilidades, como estamos trabalhando com numeros reais
    total=total+1
    xn=x+(rand()-0.5)*delta
    !write(*,*)xn,f(xn)
    if((f(xn)*f(xn)/(f(x)*f(x)))>=rand())then
      x=xn
      i=i+1
      El(i)=2/(x*(1-x))
    endif
  enddo
  Et=0
  do i=0,N
    Et=Et+El(i)
  enddo
  write(*,*)delta,'&',i/total,'&',Et/N
  !write(*,*)
end program
function f(r)
  real*8::r,f,g
  if (1.0>r .and. r>0.0 )then
    g=r*(1-r)
  else
    g=0
  endif
  f=g
end function