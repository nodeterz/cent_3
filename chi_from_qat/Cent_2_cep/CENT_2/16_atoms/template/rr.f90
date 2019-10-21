program std_mean
    implicit none
    integer :: len_a=6
    real(8) :: a(6)
    real(8) :: mean_a, std_a
    a(1) = 727.7
    a(2) = 1086.5
    a(3) = 1091.0
    a(4) = 1361.3
    a(5) = 1490.5
    a(6) = 1956.1
    mean_a = sum(a)/len_a
    std_a = sqrt(sum((a(:)-mean_a)**2)/(len_a-1.d0))
    write(*,*) a
    write(*,*) '---------------------------------'
    write(*,*) mean_a, std_a
end program std_mean
