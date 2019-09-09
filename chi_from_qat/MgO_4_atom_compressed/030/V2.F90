program main
    !use F95_LAPACK
    implicit none
    integer :: iat,nat,info, i, j
    integer, allocatable :: ipiv(:)
    character(LEN=3),allocatable :: sat(:)
    real(8),allocatable :: rat(:,:),qat(:),chi(:)
    real(8),allocatable :: chi_1(:),chi_2(:)
    real(8),allocatable :: gw_1(:),gw_2(:)
    real(8),allocatable :: hardness_1(:),hardness_2(:)
    real(8),allocatable :: a(:,:),ainv(:,:)
    real(8),allocatable :: b(:,:),binv(:,:)
    real(8),allocatable :: c(:,:),cinv(:,:)
    real(8),allocatable :: d(:,:),dinv(:,:)
    real(8),allocatable :: temp_1_1(:,:),temp_2_1(:,:)
    real(8),allocatable :: temp_1_2(:,:),temp_2_2(:,:)
    real(8),allocatable :: temp_1_3(:,:),temp_2_3(:,:)
    real(8),allocatable :: temp_1_4(:,:),temp_2_4(:,:)
    real(8),allocatable :: chi_mat_1(:,:),chi_mat_2(:,:),chi_mat_inv(:,:)
    real(8),allocatable :: mat_chi_1(:),mat_chi_2(:)
    real(8),allocatable :: qq(:)
    real(8),allocatable :: a_test(:,:),chi_test(:),qat_test(:),a_test_inv(:,:)
    integer :: chi_loop_1,chi_loop_2,chi_loop_3,chi_loop_4
    integer :: i_test
    real(8) :: coeff,chi_mean_1, chi_mean_2, chi_var_1, chi_var_2
    character(len=12) :: format_string, str_nat

    open(2,file='posinp.rzx')
    read(2,*) nat 
    allocate(sat(nat),rat(3,nat),qat(nat+1))
    allocate(chi_1(nat),chi_2(nat))
    allocate(gw_1(nat),gw_2(nat))
    allocate(hardness_1(nat),hardness_2(nat))

    qat = 0.d0
    do iat = 1 , nat
        read(2,*) sat(iat),rat(1,iat),rat(2,iat),rat(3,iat),qat(iat),gw_1(iat),gw_2(iat)&
                  ,hardness_1(iat),chi_1(iat)
    end do
    allocate(a(nat,nat),ainv(nat,nat))
    allocate(b(nat,nat),binv(nat,nat))
    allocate(c(nat,nat),cinv(nat,nat))
    allocate(d(nat,nat),dinv(nat,nat))

    allocate(temp_1_1(nat,nat),temp_2_1(nat,nat))
    allocate(temp_1_2(nat,nat),temp_2_2(nat,nat))
    allocate(temp_1_3(nat,nat),temp_2_3(nat,nat))
    allocate(temp_1_4(nat,nat),temp_2_4(nat,nat))
    allocate(chi_mat_1(nat,nat),chi_mat_2(nat,nat))
    allocate(qq(1:nat),mat_chi_1(1:nat),mat_chi_2(1:nat),chi_mat_inv(1:nat,1:nat))
    write(*,*) 'chi_2 and hardness_2, chi_var_1, chi_var_2'

    hardness_2(:) = 0.d0

    !## calculating a b c d matrices and their inverse
    call get_mat_cent1(nat,rat,gw_1,gw_2,hardness_1,hardness_2,a,b,c,d)
    call inv(a,nat,ainv)
    call inv(b,nat,binv)
    call inv(c,nat,cinv)
    call inv(d,nat,dinv)

    !## calculating inv(-(b-d*ainv*c))
    call mat_mult(ainv,c,nat,nat,nat,temp_1_1)
    call mat_mult(d,temp_1_1,nat,nat,nat,temp_2_1)
    temp_1_1 = b-temp_2_1
    call  inv(temp_1_1,nat,temp_2_1)
    
    !## calculating inv(a-c*binv*d)*c*binv 
    call mat_mult(binv,d,nat,nat,nat,temp_1_2) 
    call mat_mult(c,temp_1_2,nat,nat,nat,temp_2_2) 
    temp_1_2 = a - temp_2_2 
    call  inv(temp_1_2,nat,temp_2_2)
    call mat_mult(temp_2_2,c,nat,nat,nat,temp_1_2) 
    call mat_mult(temp_1_2,binv,nat,nat,nat,temp_2_2)

    chi_mat_2(1:nat,1:nat) = temp_2_2(1:nat,1:nat) - temp_2_1(1:nat,1:nat)
    
    !## calculating inv(b-d*ainv*c)*d*ainv
    call mat_mult(ainv,c,nat,nat,nat,temp_1_3)
    call mat_mult(d,temp_1_3,nat,nat,nat,temp_2_3)
    temp_1_3 = b-temp_2_3
    call  inv(temp_1_3,nat,temp_2_3)
    call mat_mult(temp_2_3,d,nat,nat,nat,temp_1_3)
    call mat_mult(temp_1_3,ainv,nat,nat,nat,temp_2_3)
    
    !## calculating inv(-(a-c*binv*d))
    call mat_mult(binv,d,nat,nat,nat,temp_1_4) 
    call mat_mult(c,temp_1_4,nat,nat,nat,temp_2_4) 
    temp_1_4 = a - temp_2_4 
    call  inv(temp_1_4,nat,temp_2_4)
    
    chi_mat_1(1:nat,1:nat) = temp_2_3(1:nat,1:nat) - temp_2_4(1:nat,1:nat)


    write(str_nat,'(I2.2)') 2*nat+2
    format_string = '('//trim(str_nat)//'es14.6)'
    
    call mat_mult(chi_mat_1,chi_1,nat,nat,1,mat_chi_1) 
    call  inv(chi_mat_2,nat,chi_mat_inv)
    qq(1:nat) = qat(1:nat)-mat_chi_1(1:nat)
    call mat_mult(chi_mat_inv,qq,nat,nat,1,chi_2) 
    chi_mean_1 = 2.d0*sum(chi_2(1:nat/2))/nat
    chi_mean_2 = 2.d0*sum(chi_2(nat/2+1:nat))/nat
    chi_var_1 = sqrt((2.d0/nat)*sum((chi_2(1:nat/2)-chi_mean_1)**2))
    chi_var_2 = sqrt((2.d0/nat)*sum((chi_2(nat/2+1:nat)-chi_mean_2)**2))
    write(*,format_string)chi_2,hardness_2,chi_var_1,chi_var_2
    !call mat_mult(chi_mat_1,chi_1,nat,nat,1,mat_chi_1) 
    !call mat_mult(chi_mat_2,chi_2,nat,nat,1,mat_chi_2)
    !write(*,*) 'qat_cep :'
    !write(*,format_string) mat_chi_1+mat_chi_2 
    !write(*,*) 'qat_mulliken :'
    !write(*,format_string) qat(1:nat)

!*****************************************************************************************
!Test part
    allocate(a_test(1:2*nat,1:2*nat),chi_test(1:2*nat),qat_test(1:2*nat),a_test_inv(1:2*nat,1:2*nat))
    a_test(1:nat,1:nat) = a(1:nat,1:nat)
    a_test(nat+1:2*nat,1:nat) = d(1:nat,1:nat)
    a_test(1:nat,nat+1:2*nat) = c(1:nat,1:nat)
    a_test(nat+1:2*nat,nat+1:2*nat) = b(1:nat,1:nat)
    chi_test(1:nat) = chi_1(1:nat)
    chi_test(nat+1:2*nat) = chi_2(1:nat)
    call  inv(a_test,2*nat,a_test_inv)
    call mat_mult(a_test_inv,-1.d0*chi_test,2*nat,2*nat,1,qat_test) 
    do i_test = 1 , nat
        write(*,'(4es14.6)') qat_test(i_test),qat_test(i_test+nat),qat_test(i_test)+qat_test(i_test+nat),qat(i_test)
    end do
    deallocate(a_test,chi_test,qat_test,a_test_inv)
!*****************************************************************************************
!deallocation part
    deallocate(sat,rat,qat ,chi_1,chi_2 ,gw_1,gw_2 ,hardness_1,hardness_2)
    deallocate(a,ainv ,b,binv ,c,cinv ,d,dinv)
    deallocate(temp_1_1,temp_2_1 ,temp_1_2,temp_2_2)
    deallocate(temp_1_3,temp_2_3 ,temp_1_4,temp_2_4)
    deallocate(chi_mat_1,chi_mat_2)
    deallocate(qq,mat_chi_1,mat_chi_2)
 
    
    !!allocate(ipiv(nat+1))
    !call DGETRF(nat+1,nat+1,a,nat+1,ipiv,info)
    !if(info/=0) then
    !    write(*,*) 'ERROR: DGETRF info=',info
    !    stop
    !endif
    !allocate(qq(1:nat+1))
    !qq(1:nat) = -chi(1:nat)
    !qq(1+nat) = sum(qat)
    !do iat = 1 , nat+1
    !    write(*,'(a,es14.6)') 'input -chi for CEP (no_shift) ',qq(iat)
    !end do
    !call DGETRS('N',nat+1,1,a,nat+1,ipiv,qq,nat+1,info)
    !write(*,'(a,es14.6)') 'Lagrangian multiplier before shift : ' ,qq(nat+1) 
    !write(*,*) 'qat from CEP(col.1) and qat form mulliken(col.2) (no_shift in chi) :'
    !do iat = 1 , nat
    !   write(*,'(2es14.6)') qq(iat) , qat(iat)
    !end do 

    !qq(1:nat) = -chi(1:nat)-1
    !qq(1+nat) = sum(qat)
    !do iat = 1 , nat+1
    !    write(*,'(a,es14.6)') 'input -chi for CEP (shift = +1) ',qq(iat)
    !end do
    !call DGETRS('N',nat+1,1,a,nat+1,ipiv,qq,nat+1,info)
    !write(*,'(a,es14.6)') 'Lagrangian multiplier after shift : ' ,qq(nat+1) 
    !write(*,*) 'qat from CEP(col.1) and qat form mulliken(col.2) (shift in chi) :'
    !do iat = 1 , nat
    !   write(*,'(2es14.6)') qq(iat) , qat(iat)
    !end do 
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
