program onda
  implicit real*8(a-h,o-z)
  real*8, dimension(0:10000)::El
  integer,parameter :: seed = 15000!gerando uma semente
  call srand(seed)
  nao=rand(seed)!ignorar essa linha 
  !entrada de dados
  write(*,*)"escreva delta"
  read(*,*)delta
  write(*,*)"escreva a constante da exponencial "
  read(*,*)ax
  write(*,*)"escreva o modulo"
  read(*,*)a
  write(*,*)"escreva o inicio"
  read(*,*)x0

  x=0
  i=0
  total=0
  do while(i<=2000)
    !razão entre as probabilidades, como estamos trabalhando com numeros reais
    xn=x+(rand()-0.5)*delta
    !write(*,*)xn,f(xn,ax)
    pxn=f(xn,ax,a,x0)
    px=f(x,ax,a,x0)
    if((pxn*pxn/(px*px))>=rand())then
      x=xn
      i=i+1
    endif
  enddo
  i=0
  !temos que a energia local é igual 2 na parte do numerador
  do while(i<=10000)
    !razão entre as probabilidades, como estamos trabalhando com numeros reais
    total=total+1
    xn=x+(rand()-0.5)*delta

    !write(*,*)xn,f(xn,ax)
    pxn=f(xn,ax,a,x0)
    px=f(x,ax,a,x0)
    if((pxn*pxn/(px*px))>=rand())then
      x=xn
      i=i+1
      El(i)=-(ax*(2*ax*(x-x0)**2-1))+800*((x)**(-12.0)-(x)**(-6.0))
    endif

  enddo

  Et=0
  do i=0,10000
    Et=Et+El(i)
  enddo


  write(*,*)i/total
  write(*,*)Et/10000

end program


function f(x,ax,a,x0)

  real*8::x,ax,f,a,x0
  f=exp(-((x-x0)**2*ax))*a

end function