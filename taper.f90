subroutine  taper(aa, n)
integer i, m,n
real aa(n),tapers
real tt, pi1
tapers=0.05
m = tapers*n
pi1 = 3.1415926/m
do i=1,m
    tt = 0.5*(1.-cos(i*pi1))
    aa(i) = tt*aa(i)
    aa(n-i+1) = tt*aa(n-i+1)
enddo
end subroutine
