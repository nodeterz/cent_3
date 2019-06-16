program main
    !use F95_LAPACK
    implicit none
    integer :: iat,nat,info
    real(8),allocatable :: rat(:,:),gw(:),qat(:),hardness(:),chi(:)
    real(8),allocatable :: a(:,:)
    character(LEN=3),allocatable :: sat(:)
    integer, allocatable :: ipiv(:)
    real(8),allocatable :: qq(:)
    open(2,file='posinp.rzx')

    read(2,*) nat 
    allocate(sat(nat),rat(3,nat),gw(nat),qat(nat+1),hardness(nat),chi(nat))
    qat = 0.d0
    do iat = 1 , nat
        read(2,*) sat(iat),rat(1,iat),rat(2,iat),rat(3,iat),qat(iat),gw(iat),hardness(iat)
    end do
    allocate(a(nat+1,nat+1))
    call get_amat_cent1(nat,rat,gw,hardness,a)
    do iat = 1 , nat
        chi(iat)=-1.d0*dot_product(a(iat,1:nat),qat(:))
    end do
    !do iat = 1 , nat
    !    write(*,'(a,es14.6)') 'chi from "A*q = -chi" ',chi(iat)
    !end do
    allocate(ipiv(nat+1))
    call DGETRF(nat+1,nat+1,a,nat+1,ipiv,info)
    if(info/=0) then
        write(*,*) 'ERROR: DGETRF info=',info
        stop
    endif
    allocate(qq(1:nat+1))
    qq(1:nat) = -chi(1:nat)
    qq(1+nat) = sum(qat)
    !do iat = 1 , nat+1
    !    write(*,'(a,es14.6)') 'input -chi for CEP (no_shift) ',qq(iat)
    !end do
    call DGETRS('N',nat+1,1,a,nat+1,ipiv,qq,nat+1,info)
    write(*,'(a,es14.6)') 'Lagrangian multiplier before shift : ' ,qq(nat+1) 
    !write(*,*) 'qat from CEP(col.1) and qat form mulliken(col.2) (no_shift in chi) :'
    !do iat = 1 , nat
    !   write(*,'(2es14.6)') qq(iat) , qat(iat)
    !end do 

    qq(1:nat) = -chi(1:nat)-1
    qq(1+nat) = sum(qat)
    !do iat = 1 , nat+1
    !    write(*,'(a,es14.6)') 'input -chi for CEP (shift = +1) ',qq(iat)
    !end do
    call DGETRS('N',nat+1,1,a,nat+1,ipiv,qq,nat+1,info)
    write(*,'(a,es14.6)') 'Lagrangian multiplier after shift : ' ,qq(nat+1) 
    !write(*,*) 'qat from CEP(col.1) and qat form mulliken(col.2) (shift in chi) :'
    !do iat = 1 , nat
    !   write(*,'(2es14.6)') qq(iat) , qat(iat)
    !end do 
    
end program main
!*****************************************************************************************
!subroutine get_amat_cent1(atoms,ann_arr,a)
subroutine get_amat_cent1(nat,rat,gw,hardness,a)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: hardness(nat)
    real(8), intent(out):: a(nat+1,nat+1)
    !local variables
    integer:: iat, jat
    real(8):: dx, dy, dz, r, pi, beta_iat, beta_jat, gama
    pi=4.d0*atan(1.d0)
    do iat=1,nat
        a(iat,nat+1)=1.d0
        a(nat+1,iat)=1.d0
        beta_iat=gw(iat)
        gama=1.d0/sqrt(beta_iat**2+beta_iat**2)
        a(iat,iat)=gama*2.d0/sqrt(pi)+hardness(nat)
        do jat=iat+1,nat
            dx=rat(1,jat)-rat(1,iat)
            dy=rat(2,jat)-rat(2,iat)
            dz=rat(3,jat)-rat(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            beta_jat=gw(jat)
            gama=1.d0/sqrt(beta_iat**2+beta_jat**2)
            a(iat,jat)=erf(gama*r)/r
            a(jat,iat)=a(iat,jat)
        enddo
    enddo
    a(nat+1,nat+1)=0.d0
end subroutine get_amat_cent1
