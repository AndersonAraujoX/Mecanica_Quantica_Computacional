program onda
  implicit real*8(a-h,o-z)
  real*8,dimension(-1:10001,-1:10001) ::re,im
  write(*,*)"digite o valor do potencial"
  read(*,*)v0
  N=200
  deltat=0.025
  deltax=0.2
  !write(*,*)1/deltat-1/deltax**2
  !call sleep(10)
  x0=-15.0d0
  pi=acos(-1.0d0)
  !write(*,*)pi
  soma=0.0d0
  do j=0,N
    x=(j-N*0.5)*deltax
    im(j,0)=exp(-(x-x0)*(x-x0)*0.25d0)*sin((x-x0))*(2*pi)**(-0.25d0)
    re(j,0)=exp(-(x-x0)*(x-x0)*0.25d0)*cos((x-x0))*(2*pi)**(-0.25d0)
    soma=soma+(re(j,0)*re(j,0)+im(j,0)*im(j,0))*deltax
  enddo
  !write(*,*)soma
  !meio caminho do primeiro tempo para o imaginario
  do j=0,N
    !condição de contorno    
    re(-1,i)=0.0
    re(N+1,i)=0.0
    x=(j-N*0.5)*deltax
    if(abs(x)<0.5d0)then!condição para o potencial
      im(j,0)=im(j,0)-v0*re(j,0)*deltat*0.5d0+0.25*deltat*(re(j+1,0)-2*re(j,0)+re(j-1,0))/(deltax*deltax)
    else
      im(j,0)=im(j,0)+0.25*deltat*(re(j+1,0)-2*re(j,0)+re(j-1,0))/(deltax*deltax)
    endif

  enddo
  
  do i=0,1000!tempo para não estourar a quntidade de arquivos
    !fazendo a iteração da parte real e imaginaria
    !write(*,*)i*deltat
    re(-1,i)=0.0
    re(N+1,i)=0.0
    im(-1,i)=0.0
    im(N+1,i)=0.0
    do j=0,N!p1osição
      !condição para o potencial
      x=(j-N*0.5)*deltax
      if(abs(x)<0.5d0)then
        re(j,i+1)=re(j,i)+v0*im(j,i)*deltat-0.5*deltat*(im(j+1,i)-2*im(j,i)+im(j-1,i))/(deltax*deltax)        
      else
        re(j,i+1)=re(j,i)-0.5*deltat*(im(j+1,i)-2*im(j,i)+im(j-1,i))/(deltax*deltax)
      endif
    enddo
    
    do j=0,N
      !condição para o potencial
      x=(j-N*0.5)*deltax
      if(abs(x)<0.5d0)then
        im(j,i+1)=im(j,i)-v0*re(j,i)*deltat+0.5*deltat*(re(j+1,i+1)-2*re(j,i+1)+re(j-1,i+1))/(deltax*deltax)
      else
        im(j,i+1)=im(j,i)+0.5*deltat*(re(j+1,i+1)-2*re(j,i+1)+re(j-1,i+1))/(deltax*deltax)
      endif
    enddo
    !write(*,*)re(N,i),im(N,i)
    !call sleep(1)
  enddo
  do k=0,1000!tempo
    !write(*,*)k
    do i=0,N!posição

      !po(i,2*k)=re(i,k)*re(i,k)+im(i,k+1)*im(i,k)
      !po(i,2*k+1)=re(i,k+1)*re(i,k)+im(i,k)*im(i,k)
      write(15+2*k,*)(i-N*0.5)*deltax,re(i,k)*re(i,k)+im(i,k+1)*im(i,k)!primeiro tempo
      write(16+2*k,*)(i-N*0.5)*deltax,re(i,k+1)*re(i,k)+im(i,k)*im(i,k)!segundo tempo
      !write(15+2*k,*)(i-500)*deltax,re(i,k)
      !write(16+2*k,*)(i-500)*deltax,im(i,k)

    enddo

  enddo
end program