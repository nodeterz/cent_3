program main
    !use F95_LAPACK
    implicit none
    integer :: iat,nat,info, i, j
    integer, allocatable :: ipiv(:)
    character(LEN=3),allocatable :: sat(:)
    real(8),allocatable :: rat(:,:),qat(:),chi(:)
    real(8),allocatable :: chi_1(:),chi_2(:),chi_tot(:)
    real(8),allocatable :: gw_1(:),gw_2(:)
    real(8),allocatable :: hardness_1(:),hardness_2(:)
    real(8),allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
    real(8),allocatable :: a_tot(:,:), a_tot_inv(:,:)
    real(8) :: q_tot, lag_mult
    real(8) :: epot , eref
    character :: tt_1
    
    open(2,file='posinp.rzx')
    read(2,*) nat 
    allocate(sat(nat),rat(3,nat))
    allocate(chi_1(nat),chi_2(nat))
    allocate(gw_1(nat),gw_2(nat))
    allocate(hardness_1(nat),hardness_2(nat))
    !--------------------------------------------------------------------------------------
    hardness_2(:) = 0.d0
    q_tot = 0.d0
    do iat = 1 , nat
        read(2,*) tt_1,rat(1,iat),rat(2,iat),rat(3,iat),sat(iat)
    end do
    do iat = 1 , nat
        read(2,*) chi_2(iat)
    end do
    read(2,*) tt_1,tt_1,tt_1,tt_1
    do iat = 1 , nat
        read(2,*) gw_1(iat),gw_2(iat),hardness_1(iat),chi_1(iat)
    end do
    allocate(a(nat,nat) , b(nat,nat) , c(nat,nat) , d(nat,nat))
    allocate(a_tot(2*nat+1,2*nat+1),chi_tot(2*nat+1))
    !## calculating a b c d and a_tot matrices
    call get_mat_cent1(nat,rat,gw_1,gw_2,hardness_1,hardness_2,a,b,c,d)
    a_tot(1:nat,1:nat) = a(1:nat,1:nat)
    a_tot(nat+1:2*nat,1:nat) = d(1:nat,1:nat)
    a_tot(1:nat,nat+1:2*nat) = c(1:nat,1:nat)
    a_tot(nat+1:2*nat,nat+1:2*nat) = b(1:nat,1:nat)
    a_tot(2*nat+1,1:2*nat+1)=1
    a_tot(1:2*nat+1,2*nat+1)=1
    a_tot((2*nat)+1,(2*nat)+1)=0
    do i = 1 , 2*nat+1
        write(*,*)  a_tot(i,:)
    end do
    chi_tot(1:nat) = chi_1(1:nat)
    chi_tot(nat+1:2*nat) = chi_2(1:nat)
    chi_tot(2*nat+1) = q_tot 
    do i = 1 , 2*nat+1
        write(*,*)  'chi_tot : ', chi_tot(i)
    end do
    !## CEP part
   
    allocate(a_tot_inv(2*nat+1,2*nat+1),qat(2*nat+1))
    call inv(a_tot,2*nat+1,a_tot_inv)
    do i = 1 , 2*nat+1
        write(*,*)  a_tot_inv(i,:)
    end do
    call mat_mult(a_tot_inv,chi_tot,2*nat+1,2*nat+1,1,qat)
    do i = 1 , 2*nat+1
        write(*,*) 'qat : ', qat(i)
    end do
        write(*,*)  '***********************************************************************'
    do i = 1 , nat
        write(*,*)  'sum qat : ', qat(i)+qat(i+nat)
    end do
    !## ENERGY part
    epot = 0.d0
    eref = 0.d0
    call cal_electrostatic_ann(nat,rat,gw_1,gw_2,qat,epot)
    write(*,*) 'EPOT = ', epot
    do i = 1, nat
        epot = epot + eref + chi_tot(i)*qat(i) + chi_tot(i+nat)*qat(i+nat) + &
                      0.5d0*hardness_1(i)*(qat(i)**2) + 0.5d0*hardness_2(i)*(qat(i+nat)**2)

    end do
    write(*,*) 'EPOT = ', epot
!*****************************************************************************************

!deallocation part
    
end program main
!*****************************************************************************************
!subroutine get_amat_cent1(atoms,ann_arr,a)
subroutine get_mat_cent1(nat,rat,gw_1,gw_2,hardness_1,hardness_2,a,b,c,d)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat)
    real(8), intent(in):: gw_1(nat),gw_2(nat)
    real(8), intent(in):: hardness_1(nat),hardness_2(nat)
    real(8), intent(out):: a(nat,nat),b(nat,nat),c(nat,nat),d(nat,nat)
    !local variables
    integer:: iat, jat
    real(8):: dx, dy, dz, r, pi, beta_iat_1, beta_jat_1,beta_iat_2,beta_jat_2,gama_1,gama_2,gama_3,gama_4
    pi=4.d0*atan(1.d0)
    do iat=1,nat
        beta_iat_1=gw_1(iat)
        beta_iat_2=gw_2(iat)
        gama_1=1.d0/sqrt(beta_iat_1**2+beta_iat_1**2)
        gama_2=1.d0/sqrt(beta_iat_1**2+beta_iat_2**2)
        gama_3=1.d0/sqrt(beta_iat_2**2+beta_iat_1**2)
        gama_4=1.d0/sqrt(beta_iat_2**2+beta_iat_2**2)
        a(iat,iat)=gama_1*2.d0/sqrt(pi)+hardness_1(iat)
        b(iat,iat)=gama_4*2.d0/sqrt(pi)+hardness_2(iat)
        c(iat,iat)=(gama_2+gama_3)/sqrt(pi)
        d(iat,iat)=(gama_2+gama_3)/sqrt(pi)
        do jat=iat+1,nat
            dx=rat(1,jat)-rat(1,iat)
            dy=rat(2,jat)-rat(2,iat)
            dz=rat(3,jat)-rat(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            beta_jat_1=gw_1(jat)
            beta_jat_2=gw_2(jat)
            gama_1=1.d0/sqrt(beta_iat_1**2+beta_jat_1**2)
            gama_2=1.d0/sqrt(beta_iat_1**2+beta_jat_2**2)
            gama_3=1.d0/sqrt(beta_iat_2**2+beta_jat_1**2)
            gama_4=1.d0/sqrt(beta_iat_2**2+beta_jat_2**2)
            a(iat,jat)=erf(gama_1*r)/r
            a(jat,iat)=a(iat,jat)
            b(iat,jat)=erf(gama_4*r)/r
            b(jat,iat)=b(iat,jat)
            c(iat,jat)=erf(gama_2*r)/r
            c(jat,iat)=c(iat,jat)
            d(iat,jat)=erf(gama_3*r)/r
            d(jat,iat)=d(iat,jat)
        enddo
    enddo
end subroutine get_mat_cent1
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
subroutine cal_electrostatic_ann(nat,rat,gw_1,gw_2,qat,epot)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat)
    real(8), intent(in):: gw_1(nat),gw_2(nat)
    real(8), intent(in):: qat(2*nat)
    real(8), intent(out) :: epot
    ! local variables

    integer :: iat, jat
    real(8) :: pi 
    real(8) :: dx, dy, dz, r
    real(8) :: gama_1, gama_2, gama_3, gama_4
    real(8) :: beta_iat_1, beta_iat_2, beta_jat_1, beta_jat_2

    epot = 0.d0
    pi = 4.d0*atan(1.d0)
    do iat = 1 , nat
        beta_iat_1=gw_1(iat)
        beta_iat_2=gw_2(iat)
        gama_1=1.d0/sqrt(beta_iat_1**2+beta_iat_1**2)
        gama_2=1.d0/sqrt(beta_iat_1**2+beta_iat_2**2)
        gama_3=1.d0/sqrt(beta_iat_2**2+beta_iat_1**2)
        gama_4=1.d0/sqrt(beta_iat_2**2+beta_iat_2**2)
        epot = epot + (qat(iat)**2*gama_1 + qat(iat+nat)**2*gama_4 + qat(iat)*qat(iat+nat)*(gama_2+gama_3))/sqrt(pi)
        !write(*,*) (qat(iat)**2*gama_1 + qat(iat+nat)**2*gama_4 + qat(iat)*qat(iat+nat)*(gama_2+gama_3))/sqrt(pi)
        do jat = iat+1 , nat
            dx=rat(1,jat)-rat(1,iat)
            dy=rat(2,jat)-rat(2,iat)
            dz=rat(3,jat)-rat(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            beta_jat_1=gw_1(jat)
            beta_jat_2=gw_2(jat)
            gama_1=1.d0/sqrt(beta_iat_1**2+beta_jat_1**2)
            gama_2=1.d0/sqrt(beta_iat_1**2+beta_jat_2**2)
            gama_3=1.d0/sqrt(beta_iat_2**2+beta_jat_1**2)
            gama_4=1.d0/sqrt(beta_iat_2**2+beta_jat_2**2)
            epot = epot + ( qat(iat)*qat(jat)*erf(gama_1*r) + qat(iat)*qat(jat+nat)*erf(gama_2*r) &
                    + qat(iat+nat)*qat(jat)*erf(gama_3*r) + qat(iat+nat)*qat(jat+nat)*erf(gama_4*r))/r
            write(*,*)iat,jat, ( qat(iat)*qat(jat)*erf(gama_1*r) + qat(iat)*qat(jat+nat)*erf(gama_2*r) &
                    + qat(iat+nat)*qat(jat)*erf(gama_3*r) + qat(iat+nat)*qat(jat+nat)*erf(gama_4*r))/r
        end do
    end do
end subroutine 
