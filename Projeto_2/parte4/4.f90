program onda
    implicit real*8(a-h,o-z)
    real*8 ::norm1,norm2,L
    integer::xm
    real*8, dimension(0:10000)::yl,yr
    !condições iniciais
    write(*,*)'escreva o tamanho da caixa  [-L,L]:'
    read(*,*)L
    write(*,*)'escreve a condição inicial y(0)'
    read(*,*)yz
    write(*,*)'escreve a condição inicial dy(0)/dx'
    read(*,*)yzl
    detx=(2*L)/10000.0
    j=0
    E=0
    write(*,*)"energias"
    do while(j<10)
        !resetando a normalização
        E=E+0.001
        norm2=0
        norm1=0
        do i=1,9999
            !condições iniciais para o lado esquerdo
            yl(0)=yz
            yl(1)=yzl*detx+yl(0) 
            !passo seguinte 
            yl(i+1)=2*yl(i)-yl(i-1)-2*detx*detx*E*yl(i) 
            !condições iniciais para direita
            yr(10000)=yz
            yr(9999)=yzl*detx+yr(10000)     
            !sentido contrário
            yr(10000-(i+1))=2*yr(10000-(i))-yr(10000-(i-1))-2*detx*detx*E*yr(10000-(i)) 
        enddo
        !ponto de encontro
        xm=5000
        !derivada
        const=yr(xm)/yl(xm)
        derl=const*(yl(xm+1)-yl(xm-1))/(2*detx)
        derr=(yr(xm+1)-yr(xm-1))/(2*detx)
        !plotar as seguintes energias e valores onde derivada não é pequena 
        if(abs(derr-derl)<0.0001 .and. E>0.01)then
            j=j+1
            
            write(*,*)E
            E=E+0.1  
            !normalização
            do i=0,xm
                norm2=const*const*yl(i)*yl(i)*detx+norm2
            enddo
            do i=xm,10000
                norm1=yr(i)*yr(i)*detx+norm1
            enddo
            do i=0,xm
                write(10+j,*)i*detx-L,const*yl(i)/sqrt(norm1+norm2)
            enddo
            do i=xm,10000
                write(10+j,*)i*detx-L,yr(i)/sqrt(norm1+norm2)
            enddo
        endif
        
    enddo
    end program