program onda
  implicit real*8(a-h,o-z)
  complex*8 :: psi
  integer,parameter :: seed = 15000!gerando uma semente
  !metodo da discretização
  N=500
  deltat=0.1
  pi=acos(-1.0d0)
  !temos que energia local é igual 2 na parte do numerador
  do i=0,N
    !gerando diversas ondas
    t=i*deltat
    !write(*,*)t
    do j=0,1000
      x=(j-500)*0.01
      psi=(0.3568*zexp(complex(0,0.5d0*t))/sqrt(pi*pi)-0.151388*H(x,2)*zexp(complex(0,2.5d0*t))/sqrt(pi*pi*8)&
      &+0.0786635*H(x,4)*zexp(complex(0,4.5d0*t))/sqrt(pi*pi*384))*exp(-0.5*x*x)
      r=psi*conjg(psi)
      write(i+15,*)x,r

    enddo
  enddo
end program
RECURSIVE function H(x,n) RESULT (Her)
  real*8::x,Her
  integer::n
  if(n>1)then
    Her=2*x*H(x,n-1)-2*(n-1)*H(x,n-2)
  elseif (n==1)then
    Her=2*x
  elseif(n==0)then
    Her=1
  endif
end function