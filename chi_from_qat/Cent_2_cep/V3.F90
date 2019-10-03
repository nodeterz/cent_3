program main
    !use F95_LAPACK
    implicit none
    integer :: nconf, iconf, nat, iat, sd_loop, i
    real(8) :: sd_s, q_tot, tmp_sd, e_rmse, e_rmse_old
    real(8) :: bohr2ang=0.529177210d0
    character(len=40) :: file_name, tt
    character(len=40) :: format_string, str_nat
    character(LEN=3),allocatable :: sat(:,:)
    character :: tt_1
    real(8), allocatable :: gw_1(:), gw_2(:), hardness_1(:), chi_1(:), e_aims(:), e_ref(:), hardness_2(:), e_cent(:), e_err(:)
    real(8), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:), chi_tot(:,:), qat(:,:)
    real(8), allocatable :: rat(:,:,:), a_tot(:,:,:), a_tot_inv(:,:,:)
    
    open(3,file='input.rzx')
    read(3,*) !'number of configurations: '
    read(3,*) nconf 
    read(3,*) !'number of atoms '
    read(3,*) nat
    read(3,*) !Steepest_Descent_Step
    read(3,*) sd_s 

    allocate(gw_1(nat), gw_2(nat), hardness_1(nat), hardness_2(nat), chi_1(nat))
    allocate(e_aims(nconf), e_ref(nconf), e_cent(nconf), e_err(nconf))
    allocate(sat(nconf,nat), qat(nconf,2*nat+1), chi_tot(nconf,2*nat+1))
    allocate(a(nat,nat), b(nat,nat), c(nat,nat), d(nat,nat))
    allocate(rat(nconf,3,nat), a_tot(nconf,2*nat+1,2*nat+1), a_tot_inv(nconf,2*nat+1,2*nat+1))

    read(3,*) !gw_1(iat),gw_2(iat),hardness_1(iat),chi_1(iat)
    do iat = 1, nat
            read(3,*) gw_1(iat),gw_2(iat),hardness_1(iat),chi_1(iat)
    end do
    !--------------------------------------------------------------------------------------
    q_tot = 0.d0
    do iconf = 1 , nconf
        write(tt,'(I3.3)') iconf
        file_name='posinp_'//trim(tt)//'.rzx'
        open(2,file=file_name)
        read(2,*) !E_FOR_STRUCTURE_IN_EQB., !E_AIMS_FOR_THIS_STRUCTURE
        read(2,*) e_ref(iconf), e_aims(iconf)
        write(*,*) e_ref(iconf), e_aims(iconf)
        do iat = 1 , nat
            read(2,*) tt_1,rat(iconf,1,iat),rat(iconf,2,iat),rat(iconf,3,iat),sat(iconf,iat)
            write(*,*) tt_1,rat(iconf,1,iat),rat(iconf,2,iat),rat(iconf,3,iat),sat(iconf,iat)
        end do
    end do
    rat(1:nconf,1:3,1:nat) = rat(1:nconf,1:3,1:nat)/bohr2ang
    e_aims = e_aims/27.211384500d0
    e_ref = e_ref/27.211384500d0 

    hardness_2(:) = 0.d0
    !## calculating a b c d and a_tot matrices
    do iconf = 1, nconf
        a = 0.d0
        b = 0.d0
        c = 0.d0
        d = 0.d0
        call get_mat_cent1(nat,rat(iconf,:,:),gw_1,gw_2,hardness_1,hardness_2,a,b,c,d)
        a_tot(iconf,1:nat,1:nat) = a(1:nat,1:nat)
        a_tot(iconf,nat+1:2*nat,1:nat) = d(1:nat,1:nat)
        a_tot(iconf,1:nat,nat+1:2*nat) = c(1:nat,1:nat)
        a_tot(iconf,nat+1:2*nat,nat+1:2*nat) = b(1:nat,1:nat)
        a_tot(iconf,2*nat+1,1:2*nat+1)=1
        a_tot(iconf,1:2*nat+1,2*nat+1)=1
        a_tot(iconf,(2*nat)+1,(2*nat)+1)=0

        chi_tot(iconf,1:nat) = -1.d0*chi_1(1:nat)
        chi_tot(iconf,nat+1:2*nat) = -1.d0*chi_1(1:nat)
        chi_tot(iconf,2*nat+1) = q_tot 
        call inv(a_tot(iconf,:,:),2*nat+1,a_tot_inv(iconf,:,:))
    end do
    
    !## CEP part
    !## Steepest_Descent part
    do sd_loop = 1 , Huge(sd_loop)
        do iconf = 1 , nconf
            if (sd_loop==1) then
                e_cent(iconf) = 0.d0
                call mat_mult(a_tot_inv(iconf,:,:),chi_tot(iconf,:),2*nat+1,2*nat+1,1,qat(iconf,:))
                call cal_electrostatic_ann(nat,rat(iconf,:,:),gw_1,gw_2,qat(iconf,:),e_cent(iconf))
                do i = 1, nat
                    e_cent(iconf) = e_cent(iconf) + chi_tot(iconf,i)*qat(iconf,i) + chi_tot(iconf,i+nat)*qat(iconf,i+nat) + &
                                  0.5d0*hardness_1(i)*(qat(iconf,i)**2) + 0.5d0*hardness_2(i)*(qat(iconf,i+nat)**2)
                end do
                !e_cent = e_cent*27.211384500d0
                e_err(iconf) = sqrt((e_cent(iconf) + e_ref(iconf) - e_aims(iconf))**2)
                write(*,*) 'ENERGY cent(eV) = ', e_cent(iconf)
            else
                tmp_sd = sd_s*(e_cent(iconf)+e_ref(iconf)-e_aims(iconf))/(e_err(iconf)+1.d-16)
                chi_tot(iconf,nat+1:2*nat) = chi_tot(iconf,nat+1:2*nat) - tmp_sd*qat(iconf,nat+1:2*nat)
                !chi_tot(iconf,1:2*nat) = chi_tot(iconf,1:2*nat) - tmp_sd*qat(iconf,1:2*nat)
                call mat_mult(a_tot_inv(iconf,:,:),chi_tot(iconf,:),2*nat+1,2*nat+1,1,qat(iconf,:))
                e_cent(iconf) = 0.d0
                call cal_electrostatic_ann(nat,rat(iconf,:,:),gw_1,gw_2,qat(iconf,:),e_cent(iconf))
                do i = 1, nat
                    e_cent(iconf) = e_cent(iconf) + chi_tot(iconf,i)*qat(iconf,i) + chi_tot(iconf,i+nat)*qat(iconf,i+nat) + &
                                  0.5d0*hardness_1(i)*(qat(iconf,i)**2) + 0.5d0*hardness_2(i)*(qat(iconf,i+nat)**2)
                end do
                !e_cent = e_cent*27.211384500d0
                e_err(iconf) = sqrt((e_cent(iconf) + e_ref(iconf) - e_aims(iconf))**2)
            end if
        end do
        e_rmse_old = e_rmse
        e_rmse = sqrt(sum(e_err)/nconf)
        if ((e_rmse-e_rmse_old)<0.d0) then
            sd_s = sd_s*1.001d0
        else
            sd_s = sd_s*0.99d0
        end if
        write(*,*) 'ITER, e_rmse = ',sd_loop,e_rmse, sd_s, tmp_sd, e_err
        if (abs(e_rmse-e_rmse_old)<1.d-16) exit
        if (isnan(e_rmse)) exit
    end do
    write(str_nat,'(I2.2)') 2*nat+1
    format_string = '(a6,I6.6,'//trim(str_nat)//'es14.6)'
    !write(*,*) format_string
    do iconf = 1, nconf
        write(*,trim(format_string)) 'chi : ', iconf, chi_tot(iconf,:)
        write(*,trim(format_string)) 'qat : ', iconf, qat(iconf,:)
    end do
    !write(*,*) chi_tot
    !write(*,*) qat
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
    real(8), intent(inout) :: epot
    ! local variables

    integer :: iat, jat
    real(8) :: pi 
    real(8) :: dx, dy, dz, r
    real(8) :: gama_1, gama_2, gama_3, gama_4
    real(8) :: beta_iat_1, beta_iat_2, beta_jat_1, beta_jat_2

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
!            write(*,*)iat,jat, ( qat(iat)*qat(jat)*erf(gama_1*r) + qat(iat)*qat(jat+nat)*erf(gama_2*r) &
!                    + qat(iat+nat)*qat(jat)*erf(gama_3*r) + qat(iat+nat)*qat(jat+nat)*erf(gama_4*r))/r
        end do
    end do
end subroutine 
