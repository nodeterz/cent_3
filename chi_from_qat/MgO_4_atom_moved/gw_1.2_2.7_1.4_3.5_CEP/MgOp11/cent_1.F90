program main
    !use F95_LAPACK
    implicit none
    integer :: iat,nat,info
    real(8),allocatable :: rat(:,:),gw(:),qat(:),hardness(:),chi(:)
    real(8),allocatable :: a(:,:),ainv(:,:)
    character(LEN=3),allocatable :: sat(:)
    integer, allocatable :: ipiv(:)
    real(8),allocatable :: qq(:)
    real(8) :: tt_1, tt1, tt2
    real(8):: bohr2ang=0.529177210d0, qtot = 0.d0
    real(8):: epot, epot_es, eref
    character(8) :: tt_2
    open(2,file='posinp.rzx')

    read(2,*) nat 
    allocate(sat(nat),rat(3,nat),gw(nat),qat(nat+1),hardness(nat),chi(nat+1))
    qat = 0.d0
    !do iat = 1 , nat
    !    read(2,*) sat(iat),rat(1,iat),rat(2,iat),rat(3,iat),qat(iat),ttt,gw(iat),hardness(iat),ttt
    !end do
    !----------------------------------------------------------------------------------------------
    do iat = 1 , nat
        read(2,*) tt_2,rat(1,iat),rat(2,iat),rat(3,iat),sat(iat)
    end do
    rat = rat/bohr2ang
    do iat = 1 , nat
        read(2,*) qat(iat)
    end do
    !write(*,*)'qat_1:', qat
    read(2,*) tt_2,tt_2,tt_2,tt_2 
    do iat = 1 , nat
        read(2,*) tt_1,gw(iat),hardness(iat),tt_1
    end do
    allocate(a(nat+1,nat+1),ainv(nat+1,nat+1))
    call get_amat_cent1(nat,rat,gw,hardness,a)
    do iat = 1 , nat
        chi(iat)=-1.d0*dot_product(a(iat,1:nat),qat(:))
    end do
!#CEP part    
    chi(1:nat) = -1.d0*chi(1:nat)
    chi(nat+1) = qtot
    qat = 0.d0
    call inv(a,nat+1,ainv)
    call mat_mult(ainv,chi,nat+1,nat+1,1,qat)
    !write(*,*) 'qat_2 : ',qat
    epot = 0.d0
    eref = 0.d0
    tt1=0.d0
    tt2=0.d0
    do iat=1,nat
        tt1=tt1+chi(iat)*qat(iat)
        tt2=tt2+qat(iat)**2*0.5d0*hardness(iat)
    enddo
    call cal_electrostatic_ann(nat,rat,gw,qat,epot_es)
    epot=epot_es+tt1+tt2+eref
    epot = epot*27.211384500d0

    write(*,*) 'EPOT(CENT_1)(eV) = ', epot

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
!*****************************************************************************************
subroutine inv(a,a_dim,ainv)
    implicit none
    integer,intent(in) :: a_dim
    real(8),intent(in) :: a(a_dim,a_dim)
    real(8),intent(out):: ainv(a_dim,a_dim)
    real(8)            :: work(a_dim)         ! work array for LAPACK
    integer            :: n,info,ipiv(a_dim)  ! pivot indices

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    ainv = a
    call DGETRF(a_dim,a_dim,ainv,a_dim,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    call DGETRI(a_dim,ainv,a_dim,ipiv,work,a_dim,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
end subroutine inv
!*****************************************************************************************
subroutine mat_mult(a,b,dim_1,dim_2,dim_3,ab)
    implicit none
    integer,intent(in) :: dim_1,dim_2,dim_3
    real(8),intent(in) :: a(dim_1,dim_2)
    real(8),intent(in) :: b(dim_2,dim_3)         
    real(8),intent(out):: ab(dim_1,dim_3)
    !local variables
    integer :: i , j 
    do i = 1 , dim_1
        do j = 1 , dim_3
            ab(i,j)=dot_product(a(i,:),b(:,j))
        end do
    end do
end subroutine mat_mult
!*****************************************************************************************
subroutine cal_electrostatic_ann(nat,rat,gw,qat,epot)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(inout) :: epot
    ! local variables
    integer :: iat, jat
    real(8) :: pi 
    real(8) :: dx, dy, dz, r
    real(8) :: gama
    real(8) :: beta_iat, beta_jat
    real(8) :: tt2, tt3 
    pi=4.d0*atan(1.d0)
    tt2=0.d0
    tt3=0.d0
    do iat=1,nat
        beta_iat=gw(iat)
        gama=1.d0/sqrt(beta_iat**2+beta_iat**2)
        tt2=tt2+qat(iat)**2*gama/sqrt(pi)
        do jat=iat+1,nat
            dx=rat(1,jat)-rat(1,iat)
            dy=rat(2,jat)-rat(2,iat)
            dz=rat(3,jat)-rat(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            beta_jat=gw(jat)
            gama=1.d0/sqrt(beta_iat**2+beta_jat**2)
            tt3=tt3+qat(iat)*qat(jat)*erf(gama*r)/r
        enddo
    enddo
    epot=tt2+tt3
end subroutine cal_electrostatic_ann
