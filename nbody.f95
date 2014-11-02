program nbody
!---------------Declarations
	IMPLICIT NONE
	integer :: i,j,u,h,coord,dim,r_num,n		!Dummy Index and Counters
	integer :: k,summ,count						!Random Num Counters
	real, parameter :: G = 1, e = 0.001
	real, allocatable :: posi(:,:), vel(:,:), acc(:,:), ray(:,:), mass(:)
	real, allocatable :: dumposi(:), dumvel(:)	!Dummy Initials
	real :: temp(5)
	character(len=12) :: fname
	
!---------------Random Number Generator	
	summ = 0
	do i=1,100000
	call system_clock(count)
	k = count - (count/2)*2
	summ = summ + k
	enddo
	write(*,*) summ
	
	call random_seed(k)
	do i=1,summ
	call random_number(temp(5))
	k = k + temp(5)*5
	call random_seed(k)
	enddo
	
!---------------Initial Conditions
	n = 15
	coord = 3
	dim = coord*n
	r_num = ((n-1)*n)/2
	
	allocate(mass(n))
	allocate(dumposi(dim))
	allocate(dumvel(dim))
	allocate(posi(n,coord))
	allocate(vel(n,coord))
	allocate(acc(n,coord))
	allocate(ray(n,n))
	
	do i=1,dim
	call random_number (temp(1))
	call random_number (temp(3))
		if (temp(1)<0.5) then
		temp(4) = -1
		else
		temp(4) = 1
		endif
	dumposi(i) = temp(3)*temp(4)*30
	enddo

	do i=1,dim
	call random_number (temp(1))
	call random_number (temp(3))
		if (temp(1)<0.5) then
		temp(4) = -1
		else
		temp(4) = 1
		endif
	dumvel(i) = temp(3)*temp(4)*.2
	enddo
	
	
	dumvel(1) = 0
	dumvel(2) = 0
	dumvel(3) = 0
	dumposi(1) = 0
	dumposi(2) = 0
	dumposi(3) = 0
	
	dumvel(4) = 0
	dumvel(5) = 1
	dumvel(6) = 0
	dumposi(4) = 1
	dumposi(5) = 0
	dumposi(6) = 0
	
	dumvel(7) = 0
	dumvel(8) = -.5
	dumvel(9) = 0
	dumposi(7) = -5
	dumposi(8) = 0
	dumposi(9) = 0
	
	j = 1
	do i=1,n
		do u=1,coord
		posi(i,u) = dumposi(j)
		vel(i,u) = dumvel(j)
		write(*,*) posi(i,u), vel(i,u)
		write(*,*) dumposi(j), dumvel(j)
		j = j + 1
		enddo
	enddo

	do i=1,5
		temp(i) = 0
	enddo
	
	do i=1,n
		call random_number(temp(5))
		mass(i) = temp(5)*.01
	enddo
	
	mass(1) = 1
	mass(2) = 0.1
	mass(3) = 0.1
	
!---------------File Declarations
	
	do i=1,n
	write(fname,20) i+100
	open(i,file=fname)
	enddo
	
20	format('file',i3.3,'.dat')

!---------------Write To File
	do i=1,1000000
		do j=1,n
!		write(j,*) i*e, posi(j,1), posi(j,2), posi(j,3)
		if (mod(i,5000) == 0) write(j,*) posi(j,1), posi(j,2), posi(j,3)
		enddo
		
!---------------Ray Calculations
	do j=1,n
		do u=1,n
		temp(1) = 0
			do h=1,coord
			temp(1) = temp(1) + (posi(j,h)-posi(u,h))**2
			enddo
		ray(j,u) = sqrt(temp(1))
		enddo
	enddo
	
!---------------Acceleration Calculations
	do j=1,n
	temp(1) = 0
		do h=1,coord
			do u=1,n
			if (u/=j) temp(1) = temp(1) + G*mass(u)*(posi(u,h)-posi(j,h))/(ray(j,u)**3)
			enddo
		acc(j,h) = temp(1)
		temp(1) = 0
		enddo
	enddo
	
	do j=1,n
		if (i==1) write(*,*) acc(j,1),acc(j,2),acc(j,3)
	enddo

!---------------Velocity Calculations
	do j=1,n
		do h=1,coord
		vel(j,h) = vel(j,h) + e*acc(j,h)
		enddo
	enddo

!---------------Position Calculations
	do j=1,n
		do h=1,coord
		posi(j,h) = posi(j,h) + e*vel(j,h)
		enddo
	enddo
	enddo
	
	do i=1,n
	close(i)
	enddo
	
end program nbody

