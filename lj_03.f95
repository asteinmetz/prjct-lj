program lj

implicit none

integer			:: i,j,k,t,u,v,seed
integer, parameter 	:: num = 50,steps = 1000, chunk = 100, size = 1000
real, parameter		:: tstep = 0.001, rm = 1, e = 1, m = 100
real, dimension(num,10) :: new
real, dimension(num,10) :: old
real, dimension(num) 	:: rad
character(10) 		:: mid
character(16) 		:: filename
character(8) 		:: format_string = '(I8.8)'

! PREAMBLE DATA FILES
do i = 1, (steps/chunk)
	write(mid,format_string) i*chunk
	filename = 'posi'//trim(mid)//'.dat'
	open(unit=i*chunk, file=filename)
end do
	

! GET RANDOM POSITION SEE
!write(*,*) "Seed"
!read(*,*) seed
seed = 15
call srand(seed)

! INITIAL POSITIONS AND VELOCITIES
old(:,:) = 0
new(:,:) = 0

do j = 1, num
	do k = 1, 3
		old(j,k) = rand()*50-rand()*40+500
!		old(j,k+3) = rand()
	end do
end do

! SANITY CHECK 1
write(*,*) "particle 1"
write(*,*) old(1,:)
write(*,*) "particle 2"
write(*,*) old(2,:)
write(*,*) "particle num"
write(*,*) old(num,:)

! PHYSICS
do i = 1, steps ! LOOP TIMESTEPS
	new(:,:) = 0
	do j = 1, num ! LOOP PARTICLES
		do v = 1, num ! LOOP PARTICLE
			do k = 1, 3 ! LOOP POSITION DIMENSIONS
				new(v,k) = old(v,k)+old(v,k+3)*tstep+(0.5)*old(j,k+6)*tstep**2
			end do
			rad(v) = (new(v,1)**2+new(v,2)**2+new(v,3)**2)**(0.5)
		end do

		do t = 1, num ! LOOP PARTICLES
			if (t .NE. j) then
				do k = 1, 3 ! LOOP ACCELERATION DIMENSIONS
					if (abs(rad(j)-rad(t)) .GT. (1)*rm/(2**(1/6))) then
						new(j,k+6) = new(j,k+6)+((12*e)/(m))*(rm**12/(rad(j)-rad(t))**14-rm**6/(rad(j)-rad(t))**8)*(new(j,k)-new(t,k))
					end if
				end do
			end if
		end do

		do k = 1, 3 ! LOOP VELOCITY DIMENSIONS
			if ((new(j,k) .LT. size) .AND. (new(j,k) .GT. 0)) then
				new(j,k+3) = old(j,k+3)+(0.5)*(old(j,k+6)+new(j,k+6))*tstep
			else
				new(j,k+3) = (-1)*(old(j,k+3)+(0.5)*(old(j,k+6)+new(j,k+6))*tstep)
			end if
				
		end do
		new(j,10) = (new(j,4)**2+new(j,5)**2+new(j,6)**2)**(0.5)
	end do
	old = new

if (i .EQ. 10) then
	open(unit=20, file="duck.dat")
	do j = 1, num
		write(20,*) old(j,:)
	end do
	close(20)
end if

! SANITY CHECK 2
	if (mod(i,chunk/100) .EQ. 0) then
		write(*,*) i, old(:,7)
	end if

! SAVE INTERMEDIATE DATA
	if (mod(i,chunk) .EQ. 0) then
		do j =1, num
			write(i,*) new(j,9)
		end do
		close(i)
	end if
end do

! SANITY CHECK 3
write(*,*) "particle 1"
write(*,*) new(1,:)
write(*,*) "particle 2"
write(*,*) new(2,:)
write(*,*) "particle num"
write(*,*) new(num,:)

end program lj