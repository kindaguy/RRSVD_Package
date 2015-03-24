module dyn_array


implicit none
!
!!!!!!!!!!!!!!!!!!!!!!!!!
!NEW STUFF
!!!!!!!!!!!!!!!!!!!!!!!!!
!interface
!subroutine z_rrsvd(nrA,ncA,rawA,relevant,niter,tolerance,nrU,ncU,rawU, nrS,ncS,rawS,nrVT,ncVT,rawVT) bind(C)
!
!use,intrinsic::iso_c_binding
!integer, intent(in) :: nrA,ncA,relevant,niter
!real*8, intent(in):: tolerance
!integer, intent(out) :: nrU,ncU,nrVT,ncVT,nrS,ncS
!type(c_ptr),value :: rawA
!type(c_ptr) :: rawU,rawVT,rawS
!
!end subroutine z_rrsvd
!
!subroutine tamaImport(nrA,ncA,rawA,relevant,niter,tolerance,nrUim,ncUim,rawUim, nrSim,ncSim,rawSim,nrVTim,ncVTim,rawVTim)
!integer, intent(in) :: nrA,ncA,relevant,niter
!	real*8, intent(in):: tolerance
!integer, intent(inout) :: nrUim,ncUim,nrVTim,ncVTim,nrSim,ncSim
!
!double complex, allocatable, dimension(:), intent(in) ::rawA
!double complex, allocatable, dimension(:,:), intent(inout) :: rawUim,rawVTim
!real*8,allocatable, dimension(:), intent(inout):: rawSim
!end subroutine tamaImport
!
!end interface

contains

subroutine tamaImport(nrA,ncA,rawA,relevant,niter,tolerance,nrUim,ncUim,rawUim, nrSim,ncSim,rawSim,nrVTim,ncVTim,rawVTim)
!Interface parameters
!   use Dyn_array
use,intrinsic::iso_c_binding
integer, intent(in) :: nrA,ncA,relevant,niter
	real*8, intent(in):: tolerance
integer, intent(inout) :: nrUim,ncUim,nrVTim,ncVTim,nrSim,ncSim
double complex, allocatable, dimension(:), intent(in),target ::rawA
double complex, allocatable, dimension(:,:), intent(inout) :: rawUim,rawVTim
	real*8,allocatable, dimension(:), intent(inout):: rawSim
	real*8 :: start_t,stop_t


!Local variables
!the dimensions of the svd matrices will be set by the function rrsvd
integer :: debug,i,j
integer :: nrU,ncU,nrS,ncS,nrVT,ncVT

!Here we declare the pointers to the vectors
COMPLEX (C_DOUBLE_COMPLEX),pointer::rawU(:),rawVT(:)
real(C_DOUBLE),pointer::rawS(:)
!c binding
type(c_ptr):: prawU,prawS,prawVT

print *,'provaA',rawA(1)
print *,'provaA',rawA(nrA*ncA)
print *, 'rowsA',nrA
print *, 'colsA',ncA

!first call the rrsvd function
print *,'Chiamo rrsvd'

call z_rrsvd(nrA,ncA,c_loc(rawA),relevant,niter,tolerance,nrU,ncU,prawU,nrS,ncS,prawS,nrVT,ncVT,prawVT)
print *, 'fatto'

call cpu_time(start_t)
call c_f_pointer(prawU,rawU,[nrU*ncU])
call c_f_pointer(prawS,rawS,[nrS*ncS])
call c_f_pointer(prawVT,rawVT,[nrVT*ncVT])
call cpu_time(stop_t)
print *,'Pointer assignment time = ', stop_t - start_t
!
print *, 'nrU',nrU
print *, 'ncU',ncU
print *, 'nrS',nrS
print *, 'ncS',ncS
print *, 'nrVT',nrVT
print *, 'ncVT',ncVT
!
print *, 'rawU(1)' , rawU(1)
print *, 'rawS(1)' , rawS(1)
print *, 'rawVT(1)' , rawVT(1)
!   
!!Store U   
open(5,FILE = "RRSVD_decU.dat", form='unformatted',access='stream')
write(5) nrU
write(5) ncU

do i=1,nrU*ncU



write(5) rawU(i)        

end do

close(5)
	!
	!!Store S
	open(5,FILE = "RRSVD_decS.dat", form='unformatted',access='stream')
	write(5) nrS
	write(5) ncS
	!   
	do i=1,nrS*ncS

write(5) rawS(i)        

	end do
	!
close(5)
	!
	open(5,FILE = "RRSVD_decVT.dat", form='unformatted',access='stream')
	write(5) nrVT
	write(5) ncVT
	!   
	do i=1,nrVT*ncVT

write(5) rawVT(i)       

	end do
	!
close(5)
	!
	!!Here we copy the raW array into the matrices passed as parameters
	!
	nrUim = nrU
	ncUim = ncU
	nrSim = nrS
	ncSim = ncS
	nrVTim = nrVT
	ncVTim = ncVT
	call cpu_time(start_t)
	print *,'copying  the matrix Uraw into U'
	do i=1,nrU
	do j=1,ncU

rawUim(i,j)=rawU((j-1)*nrU+(i))
	end do
	end do
	print *,'...done!'
	!
	!!  deallocate(rawA)
	! deallocate(rawU)
call z_dealloc(prawU)
	print *,'copying  the matrix Sraw into S'
	do i=1,nrS
rawSim(i)=rawS(i)
	end do

	print *,'...done!'


!    deallocate(rawS)
	call d_dealloc(prawS);
	print *,'copying  the matrix VTraw into VT'
	do i=1,nrVT
	do j=1,ncVT

rawVTim(i,j)=rawVT((j-1)*nrVT+(i))
	end do
	end do
	print *,'...done!'
!deallocate(rawVT)
	call z_dealloc(prawVT);
	call cpu_time(stop_t)
	print *,'Copy time',stop_t - start_t
	!

	end subroutine tamaImport

end module dyn_array
