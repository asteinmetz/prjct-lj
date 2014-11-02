program lj

implicit none

integer :: i,j,k,t,u,v,seed,num=5
real, parameter :: tstep = 0.1, rm = 1, e = 1, m = 1
real, dimension(5,9) 	:: new
real, dimension(5,9) 	:: old
real, dimension(5) 	:: rad
integer, dimension(3)	:: trip

write(*,*) "Seed"
read(*,*) seed
call srand(seed)

old(:,:) = 0
new(:,:) = 0

do j = 1, num
	do k = 1, 3
		old(j,k) = rand()*5
		old(j,k+3) = rand()*5
	end do
end do


write(*,*) "particle 1"
write(*,*) old(1,:)
write(*,*) "particle 2"
write(*,*) old(2,:)

do i = 1, 5 ! LOOP TIMESTEPS
	new(:,:) = 0
	do j = 1, num ! LOOP PARTICLES
		do k = 1, 3 ! LOOP POSITION DIMENSIONS
			new(j,k) = old(j,k)+old(j,k+3)*tstep+(0.5)*old(j,k+6)*tstep**2
		end do
		rad(j) = (new(j,1)**2+new(j,2)**2+new(j,3)**2)**(0.5)
		trip(:) = 0
		do t = 1, num ! LOOP PARTICLES
			if ( t .NE. j ) then
				do k = 1, 3 ! LOOP ACCELERATION DIMENSIONS
					if ( abs(new(t,k)-new(j,k)) .GT. rm/(2**(1/6)) ) then
						new(j,k+6) = new(j,k+6)+((12*e)/(m))*((-1)*rm**12/(rad(j)-rad(t))**14+rm**6/(rad(j)-rad(t))**8)*(new(t,k)-new(j,k))
					else
						trip(k) = 1
					end if
				end do
			end if
		end do
		do k = 1, 3 ! LOOP VELOCITY DIMENSIONS
			if ( trip(k) .EQ. 1 ) then
				new(j,k+3) = old(j,k+3)+(0.5)*(old(j,k+6)+new(j,k+6))*tstep
			else
				new(j,k+3) = (-1)*old(j,k+3)+(0.5)*(old(j,k+6)+new(j,k+6))*tstep
			end if
		end do
	end do
	write(*,*) "rad", rad
	old = new
end do
write(*,*) "particle 1"
write(*,*) new(1,:)
write(*,*) "particle 2"
write(*,*) new(2,:)

end program lj