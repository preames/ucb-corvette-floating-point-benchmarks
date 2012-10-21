program main

!   This collects the PSLQ output from tpslqcnk.out, then performs
!   polynomial regression on each transposed row.
!   David H Bailey     2008-01-09

use ddmodule
implicit none
integer i, i1, i2, j, job, k, m, n, n1, ndp
parameter (n = 9, ndp = 400)
integer ipvt (n), info
double precision d1, dc(0:n), dc1(0:n), dv(n), dv1(n)
double precision c(0:n), r(n), x(n,n), v(n), v1(n), t1, t2, t3, t4

job = 0

!   Set data

c(0) = 5.d0
c(1) = 2304.d0
c(2) = 118101.d0
c(3) = 1838336.d0
c(4) = 14855109.d0
c(5) = 79514880.d0
c(6) = 321537749.d0
c(7) = 1062287616.d0
c(8) = 3014530821.d0

do i1 = 0, n - 1
  dc(i1) = c(i1)
enddo

write (6, *) 'input data ='
write (6, '(4f19.1)') dc

!   Calculate polynomial regression coefficients.

do i2 = 1, n
  do i1 = 1, n
    t1 = 0.d0

    do k = 0, n - 1
      t2 = k
      if (i1 == 1 .and. i2 == 1) then
        t1 = t1 + 1.d0
      else
        t1 = t1 + t2**(i1+i2-2)
      endif
    enddo

    x(i1,i2) = t1
  enddo
enddo

do i1 = 1, n
  t1 = 0.d0

  do k = 0, n - 1
    t2 = k
    if (i1 == 1) then
      t1 = t1 + c(k)
    else
      t1 = t1 + t2**(i1-1) * c(k)
    endif
  enddo

  v(i1) = t1
enddo

call dgefa (x, n, n, ipvt, info)
call dgesl (x, n, n, ipvt, v, job)
t1 = 0.d0

do i1 = 1, n
  v1(i1) = anint (v(i1))
  t1 = max (t1, abs (v(i1) - v1(i1)))
  dv(i1) = v(i1)
  dv1(i1) = v1(i1)
enddo

write (6, *) 'max deviation from integer value =', dble (t1)
write (6, *) 'raw regression coefficients ='
write (6, '(4f19.5)') dv
write (6, *) 'rounded regression coefficients ='
write (6, '(4f19.5)') dv1
d1 = 0.d0

do k = 0, n - 1
  t1 = 0.d0
  t2 = dble (k)

  do i1 = 1, n
    if (k == 0 .and. i1 == 1) then
      t1 = t1 + dv1(i1)
    else
      t1 = t1 + dv1(i1) * t2 ** (i1 - 1)
    endif
  enddo

  dc1(k) = t1
  d1 = d1 + abs (dc1(k) - dc(k))
enddo

write (6, *) 'regenerated data based on rounded regression coefficients ='
write (6, '(4f19.2)') (dc1(k), k = 0, n - 1)
write (6, *) 'total error from original data =', d1

stop
end

subroutine dgefa(a,lda,n,ipvt,info)
use ddmodule
double precision a, t
integer lda,n,ipvt(1),info
dimension a(lda,1)

!     dgefa factors a double precision matrix by gaussian elimination.

!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .

!    on entry

!        a       double precision(lda, n)
!                the matrix to be factored.

!        lda     integer
!                the leading dimension of the array  a .

!        n       integer
!                the order of the matrix  a .

!     on return

!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.

!        ipvt    integer(n)
!                an integer vector of pivot indices.

!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgesl or dgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.

!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.

!     subroutines and functions

!     blas daxpy,dscal,idamax

!     internal variables

integer idamax,j,k,kp1,l,nm1


!     gaussian elimination with partial pivoting

info = 0
nm1 = n - 1
if (nm1 .lt. 1) go to 70

do 60 k = 1, nm1
kp1 = k + 1

!        find l = pivot index

l = idamax(n-k+1,a(k,k),1) + k - 1
ipvt(k) = l

!        zero pivot implies this column already triangularized

if (a(l,k) .eq. 0.0d0) go to 40

!           interchange if necessary

   if (l .eq. k) go to 10
      t = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
10       continue

!           compute multipliers

   t = -1.0d0/a(k,k)
   call dscal(n-k,t,a(k+1,k),1)

!           row elimination with column indexing

   do 30 j = kp1, n
      t = a(l,j)
      if (l .eq. k) go to 20
         a(l,j) = a(k,j)
         a(k,j) = t
20          continue
      call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30       continue
go to 50
40    continue
   info = k
50    continue
60 continue

70 continue
ipvt(n) = n
if (a(n,n) .eq. 0.0d0) info = n
return
end

subroutine dgesl(a,lda,n,ipvt,b,job)
use ddmodule
double precision a, b, ddot, t
integer lda,n,ipvt(1),job
dimension a(lda,1),b(1)

!     dgesl solves the double precision system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgeco or dgefa.

!     on entry

!        a       double precision(lda, n)
!                the output from dgeco or dgefa.

!        lda     integer
!                the leading dimension of the array  a .

!        n       integer
!                the order of the matrix  a .

!        ipvt    integer(n)
!                the pivot vector from dgeco or dgefa.

!        b       double precision(n)
!                the right hand side vector.

!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.

!     on return

!        b       the solution vector  x .

!     error condition

!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgeco has set rcond .gt. 0.0
!        or dgefa has set info .eq. 0 .

!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call dgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue

!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.

!     subroutines and functions

!     blas daxpy,ddot

!     internal variables

integer k,kb,l,nm1

nm1 = n - 1
if (job .ne. 0) go to 50

!        job = 0 , solve  a * x = b
!        first solve  l*y = b

if (nm1 .lt. 1) go to 30
do 20 k = 1, nm1
   l = ipvt(k)
   t = b(l)
   if (l .eq. k) go to 10
      b(l) = b(k)
      b(k) = t
10       continue
   call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
20    continue
30    continue

!        now solve  u*x = y

do 40 kb = 1, n
   k = n + 1 - kb
   b(k) = b(k)/a(k,k)
   t = -b(k)
   call daxpy(k-1,t,a(1,k),1,b(1),1)
40    continue
go to 100
50 continue

!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b

do 60 k = 1, n
   t = ddot(k-1,a(1,k),1,b(1),1)
   b(k) = (b(k) - t)/a(k,k)
60    continue

!        now solve trans(l)*x = y

if (nm1 .lt. 1) go to 90
do 80 kb = 1, nm1
   k = n - kb
   b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
   l = ipvt(k)
   if (l .eq. k) go to 70
      t = b(l)
      b(l) = b(k)
      b(k) = t
70       continue
80    continue
90    continue
100 continue
return
end

subroutine daxpy(n,da,dx,incx,dy,incy)

!     constant times a vector plus a vector.
!     jack dongarra, linpack, 3/11/78.

use ddmodule
double precision dx, dy, da
dimension dx(1),dy(1)
integer i,incx,incy,ix,iy,m,mp1,n

if(n .le. 0)return
if (da .eq. 0.0d0) return
if(incx .eq. 1 .and. incy .eq. 1)go to 20

!        code for unequal increments or equal increments
!          not equal to 1

ix = 1
iy = 1
if(incx .lt. 0)ix = (-n+1)*incx + 1
if(incy .lt. 0)iy = (-n+1)*incy + 1
do 10 i = 1,n
  dy(iy) = dy(iy) + da*dx(ix)
  ix = ix + incx
  iy = iy + incy
10 continue
return

!        code for both increments equal to 1

20 continue
do 30 i = 1,n
  dy(i) = dy(i) + da*dx(i)
30 continue
return
end

function ddot(n,dx,incx,dy,incy)
use ddmodule
double precision ddot, dx, dy, dtemp

!     forms the dot product of two vectors.
!     jack dongarra, linpack, 3/11/78.

dimension dx(1),dy(1)
integer i,incx,incy,ix,iy,m,mp1,n

ddot = 0.0d0
dtemp = 0.0d0
if(n .le. 0)return
if(incx .eq. 1 .and. incy .eq. 1)go to 20

!        code for unequal increments or equal increments
!          not equal to 1

ix = 1
iy = 1
if(incx .lt. 0)ix = (-n+1)*incx + 1
if(incy .lt. 0)iy = (-n+1)*incy + 1
do 10 i = 1,n
  dtemp = dtemp + dx(ix)*dy(iy)
  ix = ix + incx
  iy = iy + incy
10 continue
ddot = dtemp
return

!        code for both increments equal to 1

20 continue
do 30 i = 1,n
  dtemp = dtemp + dx(i)*dy(i)
30 continue
ddot = dtemp
return
end

subroutine  dscal(n,da,dx,incx)

!     scales a vector by a constant.
!     jack dongarra, linpack, 3/11/78.

use ddmodule
double precision da, dx
dimension dx(1)
integer i,incx,m,mp1,n,nincx

if(n .le. 0)return
if(incx .eq. 1)go to 20

!        code for increment not equal to 1

nincx = n*incx
do 10 i = 1,nincx,incx
  dx(i) = da*dx(i)
10 continue
return

!        code for increment equal to 1

20 continue
do 30 i = 1,n
  dx(i) = da*dx(i)
30 continue
return
end

function idamax(n,dx,incx)

!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.

use ddmodule
double precision dx, dmax
dimension dx(1)
integer i,incx,ix,n

idamax = 0
if( n .lt. 1 ) return
idamax = 1
if(n .eq. 1)return
if(incx .eq. 1)go to 20

!        code for increment not equal to 1

ix = 1
dmax = abs(dx(1))
ix = ix + incx
do 10 i = 2,n
if(abs(dx(ix)) .le. dmax) go to 5
idamax = i
dmax = abs(dx(ix))
 5    ix = ix + incx
10 continue
return

!        code for increment equal to 1

20 dmax = abs(dx(1))
do 30 i = 2,n
if(abs(dx(i)) .le. dmax) go to 30
idamax = i
dmax = abs(dx(i))
30 continue
return
end

subroutine dmxpy (n1, y, n2, ldm, x, m)
use ddmodule
double precision m, x, y, s
dimension y(*), x(*), m(ldm,*)

!   purpose:
!     multiply matrix m times vector x and add the result to vector y.

!   parameters:

!     n1 integer, number of elements in vector y, and number of rows in
!         matrix m

!     y double precision(n1), vector of length n1 to which is added 
!         the product m*x

!     n2 integer, number of elements in vector x, and number of columns
!         in matrix m

!     ldm integer, leading dimension of array m

!     x double precision(n2), vector of length n2

!     m double precision(ldm,n2), matrix of n1 rows and n2 columns

!   simpilifed by DHB for MP demo.
! ----------------------------------------------------------------------

do 200 k = 1, n1
  s = 0.d0

  do 100 j = 1, n2
    s = s + m(k,j) * x(j)
100 continue

    y(k) = y(k) + s
200 continue

return
end

