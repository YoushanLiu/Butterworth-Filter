!! -----------------------------------------------------------------------------------------
!! Author: Youshan Liu, Institute of Geology and Geophysics, Chinese Academy of Sciences
!!
!! Date: 5/25/2016
!!
!! This is the main program.
!!
!! -----------------------------------------------------------------------------------------
program butterworth

implicit none

integer it, nt, order

real(8) ampl, PI
real(8) f0, t0, t2

real(8) dt, f, f1, f2
real(8) fp, fs, rp, rs
real(8) fp1, fp2, fs1, fs2
real(8) rp1, rp2, rs1, rs2

real(8), allocatable, dimension(:) :: a, b, r, s

!!=========================================================!
!! x-the data array
nt=1001
dt=1.d-3
allocate(r(0:nt))
allocate(s(0:nt))

!! Ricker wavelet
f0 = 30.d0
t0 = 1.d0/f0
PI = 4.d0*atan(1.d0)
ampl = (PI*f0)**2
do it = 0, nt, 1
   t2 = (dble(it*dt)-t0)**2
   r(it) = (1.d0 - 2.d0*ampl*t2) * exp(-ampl*t2)
end do
!!=========================================================!
!! lowpass
fp=2.d0; fs=10.d0
rp=3.d0; rs=40.d0
!! the unit of rp and rs is dB
!call buttordlp(fp, fs, rp, rs, order)
order=3
allocate(a(0:order))
allocate(b(0:order))
call buttlp(order, dt, fp, a, b)
!Butterworth filter
s = r
call buttfilt(order, nt, a, b, s)
!Zero-phase forward and reverse digitial IIR filtering
s = r
call filtfilt(order, nt, a, b, s)
deallocate(a, b)
!=========================================================!
! highpass
fp = 10.d0; fs = 5.d0
rp = 3.d0 ; rs = 40.d0
call buttordhp(fp, fs, rp, rs, order)
allocate(a(0:order))
allocate(b(0:order))
call butthp(order, dt, fp, a, b)
!Butterworth filter
s = r
call buttfilt(order, nt, a, b, s)
!Zero-phase forward and reverse digitial IIR filtering
s = r
call filtfilt(order, nt, a, b, s)
deallocate(a, b)
!=========================================================!
! bandpass
fp1 = 25.d0; fp2 = 45.d0
fs1 = 5.d0 ; fs2 = 60.d0
rp  = 3.d0 ; rs  = 40.d0
call buttordbp(fp1, fp2, fs1, fs2, rp, rs, order)
allocate(a(0:2*order))
allocate(b(0:2*order))
call buttbp(order, dt, fp1, fp2, a, b)
!Butterworth filter
s = r
call buttfilt(2*order, nt, a, b, s)
!Zero-phase forward and reverse digitial IIR filtering
s = r
call filtfilt(2*order, nt, a, b, s)
deallocate(a, b)
!=========================================================!
! bandstop
fp1 = 5.d0 ; fp2 = 60.d0
fs1 = 25.d0; fs2 = 45.d0
rp  = 3.d0 ; rs  = 40.d0
call buttordbs(fp1, fp2, fs1, fs2, rp, rs, order)
allocate(a(0:2*order))
allocate(b(0:2*order))
call buttbs(order, dt, fp1, fp2, a, b)
!Butterworth filter
s = r
call buttfilt(2*order, nt, a, b, s)
!Zero-phase forward and reverse digitial IIR filtering
s = r
call filtfilt(2*order, nt, a, b, s)
deallocate(a, b)
!=========================================================!
deallocate(r, s)

end program butterworth