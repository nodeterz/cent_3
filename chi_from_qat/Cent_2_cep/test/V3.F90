program main
    !use F95_LAPACK
    implicit none
    integer :: nconf, iconf, nat, iat, sd_loop, i
    real(8) :: sd_s, sd_s_inp, q_tot, tmp_sd, e_rmse, e_rmse_old, e_cent, e_diff
    real(8) :: bohr2ang=0.529177210d0
    character(len=40) :: file_name, tt
    character(len=40) :: format_string, str_nat
    character(LEN=3),allocatable :: sat(:,:)
    character :: tt_1
    real(8), allocatable :: gw_1(:), gw_2(:), hardness_1(:), e_aims(:), e_ref(:), hardness_2(:), e_err(:), chi_tot(:)
    real(8), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:), qat(:,:)
    real(8), allocatable :: rat(:,:,:), a_tot(:,:,:), a_tot_inv(:,:,:)
    integer :: gw_Mg_1_loop, gw_Mg_2_loop, gw_O_1_loop, gw_O_2_loop, hardness_O_2_loop, hardness_Mg_2_loop, range_loop, i_loop 
    real(8) :: gw_Mg_1, gw_Mg_2, gw_O_1, gw_O_2, hardness_O_1, hardness_Mg_1, hardness_O_2, hardness_Mg_2, gw_step, hardness_step 
    real(8) :: range_6, progress
    real(8) :: qat_1_mean_1,qat_1_std_1, qat_1_mean_2,qat_1_std_2
    real(8) :: qat_2_mean_1,qat_2_std_1, qat_2_mean_2,qat_2_std_2
    real(8), allocatable :: qat_Mg_1(:),qat_O_1(:),qat_Mg_2(:),qat_O_2(:)
    real(8) :: grad_Mg_1, grad_O_1 , grad_Mg_2, grad_O_2 , chi_Mg_1, chi_O_1 , chi_Mg_2, chi_O_2
    
    open(3,file='input.rzx')
    open(13660716,file='err.rzx',status='UNKNOWN')
!    open(13701108,file='chi.rzx',status='UNKNOWN')
!    open(13770514,file='qat.rzx',status='UNKNOWN')

    read(3,*) tt_1, nconf 
    read(3,*) tt_1, nat
    read(3,*) tt_1, sd_s_inp 
    read(3,*) tt_1, chi_Mg_1, chi_Mg_2
    read(3,*) tt_1, chi_O_1 , chi_O_2
    read(3,*) tt_1, hardness_Mg_1 , hardness_O_1

    allocate(gw_1(nat), gw_2(nat), hardness_1(nat), hardness_2(nat))
    allocate(e_aims(nconf), e_ref(nconf), e_err(nconf))
    allocate(sat(nconf,nat), qat(nconf,2*nat+1), chi_tot(2*nat+1))
    allocate(a(nat,nat), b(nat,nat), c(nat,nat), d(nat,nat))
    allocate(rat(nconf,3,nat), a_tot(nconf,2*nat+1,2*nat+1), a_tot_inv(nconf,2*nat+1,2*nat+1))
    allocate(qat_Mg_1(nconf*nat/2),qat_O_1(nconf*nat/2),qat_Mg_2(nconf*nat/2),qat_O_2(nconf*nat/2))

    !--------------------------------------------------------------------------------------
    q_tot = 0.d0
    do iconf = 1 , nconf
        write(tt,'(I3.3)') iconf
        file_name='pos_'//trim(tt)//'.rzx'
        open(2,file=file_name)
        read(2,*) !E_FOR_STRUCTURE_IN_EQB., !E_AIMS_FOR_THIS_STRUCTURE
        read(2,*) e_ref(iconf), e_aims(iconf)
        do iat = 1 , nat
            read(2,*) tt_1,rat(iconf,1,iat),rat(iconf,2,iat),rat(iconf,3,iat),sat(iconf,iat)
        end do
    end do
    rat(1:nconf,1:3,1:nat) = rat(1:nconf,1:3,1:nat)/bohr2ang
    e_aims = e_aims/27.211384500d0
    e_ref = e_ref/27.211384500d0 
    range_loop = 2
    range_6 = range_loop**6
    gw_step = 0.5d0
    hardness_step = 0.04
    i_loop = 0.d0

    write(13660716,*)'progress  i_loop sd_loop  e_rmse     sd_s    |         GW          &
        |                     CHI                     |         HARDNESS        &
        |                                       QAT \& STD                                        |'
    
    write(13660716,*)'                                             | Mg_1  O_1 Mg_2 O_2  &
        |     Mg_1       O_1        Mg_2       O_2    | Mg_1   O_1   Mg_2  O_2  &
        |  qat_Mg_1   std_Mg_1    qat_O_1    std_O_1   qat_Mg_2   std_Mg_2    qat_O_2    std_O_2  |'
    
    write(13660716,*)'=============================================|=====================&
        |=============================================|=========================&
        |=========================================================================================|'

        gw_Mg_1 = 1.d0-gw_step
    do gw_Mg_1_loop=1, range_loop
        gw_Mg_1=gw_Mg_1 +gw_step
        gw_O_1 = 1.d0-gw_step
    do gw_O_1_loop =1, range_loop
        gw_O_1=gw_O_1  +gw_step
        gw_Mg_2 = 2.d0-gw_step
    do gw_Mg_2_loop=1, range_loop
        gw_Mg_2=gw_Mg_2 +gw_step
        gw_O_2 = 2.d0-gw_step
    do gw_O_2_loop =1, range_loop
        gw_O_2 =gw_O_2  +gw_step
        hardness_O_2 = 0.d0-hardness_step
    do hardness_O_2_loop=1, range_loop
        hardness_O_2 = hardness_O_2 + hardness_step
        hardness_Mg_2 = 0.d0-hardness_step
    do hardness_Mg_2_loop=1, range_loop
        sd_s = sd_s_inp
        hardness_Mg_2 = hardness_Mg_2 + hardness_step
        i_loop = i_loop+1
        progress = (i_loop*100.d0)/(range_6)

        gw_1(1:nat/2)   =   gw_Mg_1
        gw_1(nat/2+1:nat) = gw_O_1
        gw_2(1:nat/2)   =   gw_Mg_2
        gw_2(nat/2+1:nat) = gw_O_2
        hardness_1(1:nat/2) = hardness_Mg_1 
        hardness_1(nat/2+1:nat) = hardness_O_1 
        hardness_2(1:nat/2) = hardness_Mg_2 
        hardness_2(nat/2+1:nat) = hardness_O_2 
        
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
            call inv(a_tot(iconf,:,:),2*nat+1,a_tot_inv(iconf,:,:))
        end do !iconf

        
        !## Steepest_Descent part
        e_rmse = 0.d0
        do sd_loop = 1 , Huge(sd_loop)
            chi_tot(1:nat/2) = -1.d0*chi_Mg_1
            chi_tot(nat/2+1:nat) = -1.d0*chi_O_1
            chi_tot(nat+1:nat+nat/2) = -1.d0*chi_Mg_2
            chi_tot(nat+nat/2+1:2*nat) = -1.d0*chi_O_2
            chi_tot(2*nat+1) = q_tot 
            e_err = 0.d0
            e_diff= 0.d0
            !## CEP part
            do iconf = 1, nconf
                call mat_mult(a_tot_inv(iconf,:,:),chi_tot,2*nat+1,2*nat+1,1,qat(iconf,:))
                e_cent = 0.d0
                call cal_electrostatic_ann(nat,rat(iconf,:,:),gw_1,gw_2,qat(iconf,:),e_cent)
                do i = 1, nat
                    e_cent = e_cent + chi_tot(i)*qat(iconf,i) + chi_tot(i+nat)*qat(iconf,i+nat) + &
                                  0.5d0*hardness_1(i)*(qat(iconf,i)**2) + 0.5d0*hardness_2(i)*(qat(iconf,i+nat)**2)
                end do
                e_diff = e_diff + e_cent + e_ref(iconf) - e_aims(iconf)
                e_err = e_err + (e_cent + e_ref(iconf) - e_aims(iconf))**2
            end do! iconf in sd_loop
            grad_Mg_1 = 0.d0 
            grad_O_1  = 0.d0 
            grad_Mg_2 = 0.d0 
            grad_O_2  = 0.d0 
            do iconf = 1 , nconf
                grad_Mg_1 = grad_Mg_1 + sum(qat(iconf,1:nat/2))
                grad_O_1  = grad_O_1  + sum(qat(iconf,nat/2+1:nat))
                grad_Mg_2 = grad_Mg_2 + sum(qat(iconf,nat+1:nat+nat/2))
                grad_O_2  = grad_O_2  + sum(qat(iconf,nat+nat/2+1:2*nat))
            end do
            !chi_Mg_1 = chi_Mg_1 + sd_s*sign(1.d0,e_diff)*grad_Mg_1
            !chi_O_1 = chi_O_1 + sd_s*sign(1.d0,e_diff)*grad_O_1
            chi_Mg_2 = chi_Mg_2 + sd_s*sign(1.d0,e_diff)*grad_Mg_2
            chi_O_2 = chi_O_2 + sd_s*sign(1.d0,e_diff)*grad_O_2

            e_rmse_old = e_rmse
            e_rmse = sqrt(sum(e_err)/nconf)*27.211384500d0
            if ((e_rmse-e_rmse_old)<0.d0) then
                    sd_s = sd_s*1.05d0
                else
                    sd_s = sd_s*0.50d0
            end if
            if ((abs(e_rmse-e_rmse_old)<1.d-16)) exit
            !if ((abs(e_rmse-e_rmse_old)<1.d-16) .and. (e_rmse<1.d0)) exit
            !if (abs(e_rmse)<1.d-3) exit
            !if (abs(e_rmse)<1.d-4) exit
            if (isnan(e_rmse)) then
                write(*,*) 'ERROR :',e_err, tmp_sd
                exit
            end if
        end do ! sd_loop
        do iconf = 1 , nconf
           qat_Mg_1((iconf-1)*nat/2+1:iconf*nat/2)=qat(iconf,1:nat/2)   
           qat_O_1((iconf-1)*nat/2+1:iconf*nat/2)=qat(iconf,nat/2+1:nat)   
           qat_Mg_2((iconf-1)*nat/2+1:iconf*nat/2)=qat(iconf,nat+1:nat+nat/2)   
           qat_O_2((iconf-1)*nat/2+1:iconf*nat/2)=qat(iconf,nat+nat/2+1:2*nat)   
        end do
        call std_mean(qat_Mg_1(1:nconf*nat/2),nconf*nat/2,qat_1_mean_1,qat_1_std_1)
        call std_mean(qat_O_1,nconf*nat/2,qat_1_mean_2,qat_1_std_2)
        call std_mean(qat_Mg_2,nconf*nat/2,qat_2_mean_1,qat_2_std_1)
        call std_mean(qat_O_2,nconf*nat/2,qat_2_mean_2,qat_2_std_2)

        write(13660716,'(f7.3,2i8,2es11.3,a2,4f5.2,a2,4es11.3,a2,4f6.3,a2,8es11.3,a2)')& 
            progress ,i_loop, sd_loop, e_rmse, sd_s,'|', &
            gw_Mg_1, gw_O_1, gw_Mg_2, gw_O_2,'|',&
            chi_Mg_1, chi_O_1, chi_Mg_2, chi_O_2,'|',&
            hardness_Mg_1,hardness_O_1,hardness_Mg_2,hardness_O_2,'|',&
            qat_1_mean_1,qat_1_std_1, qat_1_mean_2,qat_1_std_2,&
            qat_2_mean_1,qat_2_std_1, qat_2_mean_2,qat_2_std_2,'|'
        call flush(13660716)
                                                      
    end do !gw_Mg_1_loop
    end do !gw_O_1_loop 
    end do !gw_Mg_2_loop
    end do !gw_O_2_loop 
    end do !hardness_O_2_loop
    end do !hardness_Mg_2_loop
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
!*****************************************************************************************
subroutine std_mean(a,len_a,mean_a,std_a)
    implicit none
    integer, intent(in) :: len_a
    real(8), intent(in) :: a(len_a)
    real(8), intent(out):: mean_a, std_a
    mean_a = sum(a)/len_a
    std_a = sqrt(sum((a(:)-mean_a)**2)/(len_a-1.d0))
end subroutine
