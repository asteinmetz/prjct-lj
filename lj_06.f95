program lj

implicit none

integer			:: i,j,k,t,u,v,seed
integer, parameter 	:: num=100,steps=1000000,chunk=100000
real			:: size=7,diff,repul,attra,energy,cool=-0.4
real, parameter		:: tstep=0.000000001,rm=1,e=299000,m=1,epermass=299000
real, dimension(num,11) :: new
real, dimension(num,11) :: old
character(10) 		:: mid
character(16) 		:: filename
character(8) 		:: format_string='(I8.8)'

! PREAMBLE DATA FILES
do i = 1, (steps/chunk)
	write(mid,format_string) i*chunk
	filename = 'posi'//trim(mid)//'.dat'
	open(unit=i*chunk, file=filename)
end do

! GET RANDOM POSITION SEED
!write(*,*) "Seed"
!read(*,*) seed
seed = 45
call srand(seed)

! INITIAL POSITIONS AND VELOCITIES
old(:,:) = 0
new(:,:) = 0

do j = 1, num
	do k = 1, 3
		old(j,k) = rand()*size
		old(j,k+3) = rand()*10-rand()*10
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
		do k = 1, 3 ! LOOP POSITION DIMENSIONS
			new(j,k) = old(j,k)+old(j,k+3)*tstep+(0.5)*old(j,k+6)*tstep**2
		end do
	end do

	do j = 1, num ! LOOP PARTICLES
		do t = 1, num ! LOOP PARTICLES
			if (t .NE. j) then
				diff = ((old(j,1)-old(t,1))**2+(old(j,2)-old(t,2))**2+(old(j,3)-old(t,3))**2)**(0.5)
				do k = 1, 3 ! LOOP ACCELERATION DIMENSIONS
					if ((diff .GT. (0.1)) .AND. (diff .LT. (2.25))) then ! COMPUTATIONAL CUTS
						repul = rm**12/diff**12
						attra = rm**6/diff**6
						new(j,11)  = new(j,11)+e*(repul-2*attra)
						new(j,k+6) = new(j,k+6)+(12*epermass/diff**2)*(repul-attra)*(new(j,k)-new(t,k))
					end if
				end do
			end if
		end do
	end do

	do j = 1, num ! LOOP PARTICLES
		do k = 1, 3 ! LOOP VELOCITY DIMENSIONS
			if ((old(j,k+3) .NE. old(j,k+3))) then ! REMOVES IMPOSSIBLE VALUES
				old(j,k) = rand()*size
				old(j,k+3) = rand()*10-rand()*10
				old(j,k+6) = 0
			end if
				
			if ((new(j,k) .LE. 0)) then ! FLOOR SCATTERING
				new(j,k+3) = (cool)*(old(j,k+3)+(0.5)*(old(j,k+6)+new(j,k+6))*tstep)
				new(j,k) = size/100000

			else if ((new(j,k) .GE. size)) then ! CEILING SCATTERGING
				new(j,k+3) = (-1)*(old(j,k+3)+(0.5)*(old(j,k+6)+new(j,k+6))*tstep)
				new(j,k) = size-size/100000

			else if ((new(j,k) .LT. size) .AND. (new(j,k) .GT. 0)) then ! INSIDE THE BOX NORMAL MOTION
				new(j,k+3) = old(j,k+3)+(0.5)*(old(j,k+6)+new(j,k+6))*tstep

			else ! CATCH-ALL FOR NONSENSE
				new(j,k) = rand()*size
				new(j,k+3) = rand()*10-rand()*10
				new(j,k+6) = 0
			end if
		end do
		new(j,10) = (new(j,4)**2+new(j,5)**2+new(j,6)**2)**(0.5)
	end do

	old(:,:) = new(:,:)
	if (i .EQ. 5000000) then
		size = 100
!		cool = -1
	end if

! INTERMEDIATE DUCK FILE
	if (i .EQ. 1) then
		open(unit=20, file="posiduck.dat")
		do j = 1, num
			write(20,*) old(j,:)
		end do
		close(20)
	end if

! SANITY CHECK 2
	if (mod(i,chunk/50) .EQ. 0) then
		write(*,*) i*tstep,new(1,:)
	end if

! SAVE INTERMEDIATE DATA
	if (mod(i,chunk) .EQ. 0) then
		do j =1, num
			write(i,*) new(j,:)
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
