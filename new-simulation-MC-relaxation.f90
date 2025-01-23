program MC_sim
implicit none 
integer :: N , Nc   !!!! Nc is number of configurations to be considered
double precision, allocatable :: x(:), y(:), z(:), R(:,:), V(:,:)
double precision, allocatable :: x1(:), y1(:), z1(:)
double precision :: sigma, L, e , Tc, box !kb=1.380649*10**(-23), Na=6.033*10**(23)
double precision  :: cond, Etotal, dE,dE_cond, rmin, rmax, Et_prev, tmax     
integer :: i, j, m, k, attempt_max, attempt
real (kind=8) :: kb, Na, B_factor 
real(8), allocatable :: rand_array(:)
double precision :: Rx, Ry, Rz

kb=1.380649d-23
Na=6.033d23
attempt_max=1000

! reading input from input.dat file
open(unit=10, file='input.dat',  status='old', action='read')

read(10, *) Nc, N,dE_cond , tmax, sigma, e, Tc

close(10)
allocate(x(N), y(N), z(N), R(N, N), V(N, N))
allocate(x1(N), y1(N), z1(N))


rmax=sigma*2.25d0
rmin= 3.0d0
L=10.229d0*sigma
box= L*L*L

!##### initial configuration and initial enegry 

allocate(rand_array(3*N))
attempt=0
10 call random_number(rand_array)

do i=1, N-1
x(i)= L*rand_array(3*(i-1) + 1)
y(i)=L*rand_array(3*(i-1) + 2)
z(i)=L*rand_array(3*(i-1) + 3)


do j=i+1, N


x(j)= L*rand_array(3*(j-1) + 4)

y(j)=L*rand_array(3*(j-1) + 5)

z(j)=L*rand_array(3*(j-1) + 6)



Rx=x(i)-x(j)
Ry=y(i)-y(j)
Rz=z(i)-z(j)

Rx= Rx - box*dfloat(idint(Rx/(box/2.d0)))
Ry= Ry - box*dfloat(idint(Ry/(box/2.d0)))
Rz= Rz - box*dfloat(idint(Rz/(box/2.d0)))

R(i,j)=sqrt(Rx**2 + Ry**2 +Rz**2)
end do
end do 

do i=1, N-1
do j=i+1, N

if (R(i,j) < rmin .and. attempt < attempt_max) then
attempt=attempt+1
goto 10
else if (attempt > attempt_max-1) then

write(*,*) " warning: Too many attempts to get valid positions."

!####### MC relaxation of the initial position #####
call pot(N,e,sigma,rmax, box, R, Etotal)
Et_prev=Etotal
dE=Etotal

do while (abs(dE)>dE_cond)
call MCposition(tmax,N,L,rmin ,x, y, z,R, x1, y1, z1)
call pot(N,e,sigma,rmax, box, R, Etotal)

dE=Etotal-Et_prev

B_factor=exp(-dE/(Tc))

call random_number(cond)
if (dE<= 0 .or. B_factor >= cond ) then
Et_prev=Etotal
x=x1
y=y1
z=z1
end if
end do  
exit 

end if

end do 
exit
end do






open(unit = 11, file="initial-conf.xyz", status="replace", action="write")
write(11 ,*) N
write(11,*)
write(*,*) N
write(*,*)
do i=1,N
write(11,*) 20.0d0, x(i), y(i), z(i)


!write(*,*) i , R(i,:)

end do
close(11)
print *, "initial configuration Data stored to 'initial-conf.xyz'"
 
! potential energy calcualtion 

call pot(N,e,sigma,rmax, box, R, Etotal)

Et_prev=Etotal


open(unit= 12, file='Enegry-output.res', status="replace", action="write")

m=0
write(12, *) m, Et_prev*0.008314





!!!!! Generation of new MC configurations !!!!!!!!!!!
dE=1
!do while (dE>-dE_cond) 
do i=1, Nc  !!!!! N of new configurations to be considered 

call MCposition(tmax, N, L,rmin, x, y, z, R, x1, y1, z1)

call pot(N,e,sigma,rmax, box, R, Etotal)    


!write(*,*) "Energy in J/mol", Etotal*0.008314


dE=Etotal-Et_prev

B_factor=exp(-dE/(Tc))

call random_number(cond)
if (dE<= 0 .or. B_factor >= cond ) then 
Et_prev=Etotal
m=m+1
write(12, *) m , Et_prev*0.008314

write(*,*) "Energy in J/mol", Etotal*0.008314, "dE is= ", dE

x=x1
y=y1
z=z1


end if 


end do 


end program MC_sim





Subroutine pot(N,e,sigma,rmax,box, R, Etotal)
integer, intent(in) :: N
double precision, intent(in) :: R(N,N) , e , box, sigma, rmax
double precision, intent(out) :: Etotal
double precision ::  V(N,N)

Etotal=0.d0

do i=1, N-1
do j=i+1, N
if ((i/=j) .AND. (R(i,j)<rmax)) then 
V(i,j)=4*e*( (sigma/R(i,j))**12 - (sigma/R(i,j))**6 )
Etotal=Etotal+V(i,j)

end if
end do
end do



end subroutine pot


subroutine MCposition(tmax,N,L,rmin ,x, y, z,R, x1, y1, z1)
integer, intent(in) :: N
double precision , intent(in) :: L, rmin
double precision , intent(in) :: x(N), y(N), z(N)
double precision , intent(out) :: R(N,N), x1(N), y1(N), z1(N)
double precision :: f,  Rx, Ry, Rz, tmax
integer :: i, ind, attempts
double precision :: rand
real(8), allocatable :: rand_array(:)

allocate(rand_array(3))


!Select a random particle     
call random_number(rand)
ind=int(rand*N)+1 
!write(*,*) "Moving particle: ", ind

!Current position
!12 x1=x
x1=x
y1=y
z1=z


call random_number(rand_array)

x1(ind)= x1(ind)+ (2.0*rand_array(1)-1.0)*tmax        !Generates a random displacement in the range [-tmax, +tmax]    
y1(ind)= y1(ind)+ (2.0*rand_array(2)-1.0)*tmax
z1(ind)= z1(ind)+ (2.0*rand_array(3)-1.0)*tmax 

!Applying Boundary Condition on selected random particle
if (x1(ind)<0) then 
x1(ind)=x1(ind)+L
else if (x1(ind)>L) then
x1(ind)=x1(ind)-L
end if


if (y1(ind)<0) then
y1(ind)=y1(ind)+L
else if (y1(ind)>L) then
y1(ind)=y1(ind)-L
end if

if (z1(ind)<0) then
z1(ind)=z1(ind)+L
else if (z1(ind)>L) then
z1(ind)=z1(ind)-L
end if


do i=1, N-1
do j=i+1, N

Rx=x1(i)-x1(j)
Ry=y1(i)-y1(j)
Rz=z1(i)-z1(j)

Rx= Rx - box*dfloat(idint(Rx/(L/2.d0)))
Ry= Ry - box*dfloat(idint(Ry/(L/2.d0)))
Rz= Rz - box*dfloat(idint(Rz/(L/2.d0)))

R(i,j)=sqrt(Rx**2 + Ry**2 +Rz**2)

end do
end do



end subroutine MCposition



