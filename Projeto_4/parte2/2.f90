program onda
  implicit real*8(a-h,o-z)
  real*8, dimension(0:10000)::El
  integer,parameter :: seed = 15000!gerando uma semente
  call srand(seed)
  nao=rand(seed)!ignorar essa linha 

  !entrada de dados
  write(*,*)"escreva delta"
  read(*,*)delta

  write(*,*)"escreva frequencia x"
  read(*,*)wx

  write(*,*)"escreva frequencia y"
  read(*,*)wy

  write(*,*)"escreva constante x"
  read(*,*)ax

  write(*,*)"escreva constante y"
  read(*,*)ay

  x=0
  y=0
  i=0
  total=0

  !temos que a energia local é igual 2 na parte do numerador
  do while(i<=10000)
    !razão entre as probabilidades, como estamos trabalhando com numeros reais
    total=total+1
    xn=x+(rand()-0.5)*delta
    yn=y+(rand()-0.5)*delta

    !write(*,*)xn,f(xn)
    pxn=f(xn,yn,ax,ay)
    px=f(x,y,ax,ay)
    
    if((pxn*pxn/(px*px))>=rand())then
      x=xn
      y=yn
      i=i+1
      El(i)=-2*(ax*ax*x*x+ay*ay*y*y)+ax+ay+wx*wx*x*x*0.5+wy*wy*y*y*0.5
    endif

  enddo

  Et=0
  do i=0,10000
    Et=Et+El(i)
  enddo


  write(*,*)'fracao de acerto',i/total
write(*,"(1F5.2,1A,1F5.2,1A,1F5.2,1A,1F5.2,1A,1F5.2,1A,1F5.2,1A,1F5.2)")wx,' & ',wy,' & ',ax,' & ',ay,' & ',0.5*wx+0.5*wy,&
&' & ',Et/10000

end program


function f(x,y,ax,ay)

  real*8::x,y,ax,ay,f
  f=exp(-(x*x*ax+ay*y*y))

end function