	program timeWrap

!	use lapack95
	use dyn_array

	implicit none

	integer:: i,info
	integer:: nrA 
	integer:: ncA
	integer:: nrU,ncU,nrVT,ncVT,nrS,ncS
	integer:: niter=4
	real*8::tolerance = 1e-3
	integer :: relevant = 100
	real*8 :: start_time, stop_time
	double complex, allocatable, dimension(:,:) :: rawUin, rawVTin
	real*8, allocatable, dimension(:) :: rawSin
	double complex, allocatable,dimension(:,:):: A
	double complex, allocatable, dimension(:) :: rawA
	real*8, allocatable, dimension(:)::ww
	!	double complex, allocatable, dimension(:) :: rawU,rawVT,rawS

    integer :: lwork
    double complex, allocatable, dimension(:) :: work
    real*8, allocatable, dimension(:) :: rwork

    character(Len=100) :: arg
    

    call getarg(1,arg)
    write (*,*) arg
  
    call cpu_time(start_time)
	
    open(5,FILE = arg, form='unformatted',access='stream')

	read(5) nrA
	read(5) ncA
	print *, 'nrA=' ,nrA
	print *, 'ncA=' ,ncA

allocate(rawA(nrA*ncA))
	do i=1,nrA*ncA

read(5) rawA(i)		

	end do

	close(5)
call cpu_time(stop_time)
	print *, 'Reading time: ',stop_time - start_time
	print *,'rawA(1)'  , rawA(1)
	nrU=nrA
	ncU=ncA
	nrS=ncA
	ncS=ncA
	nrVT=ncA
	ncVT=ncA

	allocate(raWUin(nrU,ncU))
	allocate(rawSin(nrS))
allocate(rawVTin(nrVT,ncVT))
allocate(A(nrA,ncA))
allocate(ww(ncA-1))
A=reshape(rawA,(/nrA,ncA/))
print *,'A reshaped(1,1): ',A(1,1)
call cpu_time(start_time)
!allocate the rwork array for zgesvd
allocate(rwork(5*nrA))
allocate(work(1))
!Set this for work dimension query
lwork = -1
call zgesvd('A','A',nrA,ncA,A,nrA,&
rawSin,rawUin,nrA,rawVTin,ncA,work,lwork,rwork,info)
lwork = work(1)
deallocate(work)
allocate(work(lwork))
call zgesvd('A','A',nrA,ncA,A,nrA,&
rawSin,rawUin,nrA,rawVTin,ncA,work,lwork,rwork,info)

deallocate(work)
deallocate(rwork)



!call gesvd(A,rawSin,rawUin,rawVTin,ww,'N',info)
call cpu_time(stop_time)
print *,'S(1): ',rawSin(1)
print *,'Time Standard SVD: ',stop_time - start_time

call cpu_time(start_time)
	print *,'Call to tamaImport'
call tamaImport(nrA,ncA,rawA,relevant,niter,tolerance,nrU,ncU,rawUin, nrS,ncS,rawSin,nrVT,ncVT,rawVTin)
print *,'tamaImport competed!'
call cpu_time(stop_time)
print *,'tamaImport time: ',stop_time - start_time

print *,'rawSin(1)',rawSin(1)
	end program
