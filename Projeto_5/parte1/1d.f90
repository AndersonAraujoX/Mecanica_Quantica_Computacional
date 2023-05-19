program onda
  implicit real*8(a-h,o-z)
  complex*8 :: psi
  integer,parameter :: seed = 15000!gerando uma semente
  !metodo da discretização
  write(*,*)"digite o valor esperado da energia"
  read(*,*)E0

  deltat=0.1

  pi=acos(-1.0d0)
  !temos que energia local é igual 2 na parte do numerador
  do i=0,500
    !gerando diversas ondas
    t=i*deltat
    !write(*,*)t
    do j=0,1000
      x=(j-500)*0.01
      psi=0.0
      do n=0,10
        psi=psi+sqrt(E0**(n)/gamma(n+1.0)*exp(-E0))*H(x*1.0d0,n)*exp(-0.5*x*x)*zexp(complex(0.0,-(n+0.5d0)*t))&
        &/sqrt(pi*2*2**n*gamma(n+1.0))
      enddo
      !r=psi*conjg(psi)
      r=abs(psi)*abs(psi)
      write(i+15,*)x,r,abs(psi)*abs(psi)
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
