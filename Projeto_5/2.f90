program onda
  implicit real*8(a-h,o-z)
  real*8,dimension(0:10000,0:10000) ::re,im,po
  write(*,*)"digite o valor do potencial"
  read(*,*)v0
  N=5000
  deltat=0.01
  deltax=0.01
  x0=15
  pi=acos(-1.0d0)
  do j=0,N
    x=(j-N*0.5)*deltax
    im(j,0)=exp(-(x-15)*(x-15)*0.25)*sin((x-x0))/(2*pi)**(0.25)
    re(j,0)=exp(-(x-15)*(x-15)*0.25)*cos((x-x0))/(2*pi)**(0.25)
  enddo
  !temos que energia local é igual 2 na parte do numerador
  do i=0,int(N*0.2)!tempo
    !gerando diversas
    do j=1,N-1!posição
      !condição para o potencial

      if(abs(x)<1)then
        re(j,i+1)=re(j,i)-v0*im(j,i)+0.5*deltat*(im(j+1,i)-2*im(j,i)+im(j-1,i))/(deltax*deltax)        
      else
        re(j,i+1)=re(j,i)+0.5*deltat*(im(j+1,i)-2*im(j,i)+im(j-1,i))/(deltax*deltax)
      endif

    enddo
    do j=1,N-1
      !condição para o potencial

      if(abs(x)<1)then
        im(j,i+1)=im(j,i)+v0*re(j,i)-0.5*deltat*(re(j+1,i+1)-2*re(j,i+1)+re(j-1,i+1))/(deltax*deltax)
      else
        im(j,i+1)=im(j,i)-0.5*deltat*(re(j+1,i+1)-2*re(j,i+1)+re(j-1,i+1))/(deltax*deltax)
      endif

    enddo
  enddo
  do k=0,int(N*0.5)!tempo

    do i=1,N-1!posição

      po(i,2*k)=re(i,k)*re(i,k)+im(i,k+1)*im(i,k)
      po(i,2*k+1)=re(i,k+1)*re(i,k)+im(i,k)*im(i,k)

      write(15+2*k,*)po(i,2*k),i*deltax
      write(16+2*k,*)po(i,2*k+1),i*deltax

    enddo

  enddo
end program