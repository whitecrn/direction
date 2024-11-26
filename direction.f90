module Math_Libraries
contains
recursive subroutine Fast_Sort(A,sort_type)
        real :: A(:)
        character(len=*),optional :: sort_type
        integer :: head,tail,i,s
        real :: temp
        real :: pivot
        if (.not. present(sort_type)) then
                sort_type='asc'
        end if
        head=0
        tail=size(A)+1
        pivot=A(1)
loop1:  do while (.true.)
                tail=tail-1
                head=head+1
tail_loop:      do while (.true.)
                        if (A(tail) <= pivot) then
                                exit tail_loop
                        end if
                        tail=tail-1
                end do tail_loop

head_loop:      do while (.true.) 
                        if (A(head) > pivot) then
                                exit head_loop
                        end if
                        head=head+1
                        if (head > size(A)) then
                                exit head_loop
                        end if
                end do head_loop

                if (head < tail) then
                        temp=A(head)
                        A(head)=A(tail)
                        A(tail)=temp

                else if (head == tail) then
                        temp=A(head)
                        A(head)=A(1)
                        A(1)=temp

                        exit loop1
                else if (head > tail) then
                        temp=A(tail)
                        A(tail)=A(1)
                        A(1)=temp

                        exit loop1
                end if
        end do loop1
                if (size(A(:tail-1)) > 1) then
                        call Fast_Sort(A(:tail-1))
                end if

                if (size(A(tail+1:)) > 1) then
                        call Fast_Sort(A(tail+1:))
                else 
                        return
                end if

        if (trim(adjustl(sort_type)) == 'desc') then
                call reverse_array(A)
        end if
        

contains
        
subroutine reverse_array(array)
        implicit none
        real :: array(:)
        integer :: s,i,temp

        s=size(array)

        do i=1,s/2
                temp=array(i)
                array(i)=array(s+1-i)
                array(s+1-i)=temp
        end do
end subroutine reverse_array

end subroutine Fast_Sort

subroutine minimum_image(displacement,lattice_constants)
        real :: displacement
        real :: lattice_constants
        displacement = displacement - nint(displacement / lattice_constant) * lattice_constant
end subroutine minimum_image

function cosin(x,y)
        real,dimension(3) :: x,y
        real :: cosin
        cosin = (dot_product(x,y))/(sqrt(x(1)**2+x(2)**2+x(3)**2)*sqrt(y(1)**2+y(2)**2+y(3)**2))
end function cosin

end module Math_Libraries

program diffusion
use omp_lib
use Math_Libraries
implicit none
type :: atom
        integer :: atom_type
        real :: x(3)
end type atom
integer :: i,j,k,l,m,n
integer :: alive
integer :: step,min_index,max_step
integer,dimension(5) :: atom_number
real,parameter :: cutoff_r=1.5
real,parameter :: cutoff_r2=2.9
real,dimension(6) :: x_b,x_Li
real :: ox,oy,oz,o_r,time,a,b,c,lix,liy,liz,min_angle
real,dimension(18) :: cos_angle
real,dimension(6) :: vector
type(atom),allocatable :: R(:,:),LiO(:)
integer,allocatable :: if_diff(:,:)
integer,allocatable :: co(:,:),direction(:,:)
real,allocatable :: dx(:,:),dy(:,:),dz(:,:),dr(:),d_new(:,:,:)
real,dimension(3,18) :: group_dirs
real,dimension(3) :: dir,vec,dirction
real,dimension(6) :: statistic_direction
step=0
atom_number(:)=0
!-----------------------------------------------------------------------|
!                                                                       |
!                       define all the possible direction               |
!                                                                       |
!-----------------------------------------------------------------------|
group_dirs(:,1:6)=reshape([&
1.0,0.0,0.0,-1.0,0.0,0.0,&
0.0,1.0,0.0,0.0,-1.0,0.0,&
0.0,0.0,1.0,0.0,0.0,-1.0 &
],[3,6])

group_dirs(:,7:12)=reshape([&
1.0,1.0,0.0,-1.0,-1.0,0.0,&
1.0,0.0,1.0,1.0,0.0,-1.0,&
0.0,1.0,1.0,0.0,1.0,-1.0 &
],[3,6])/sqrt(2.0)

group_dirs(:,13:18)=reshape([&
1.0,1.0,1.0,-1.0,1.0,1.0,&   
1.0,-1.0,1.0,1.0,1.0,-1.0,&  
-1.0,-1.0,-1.0,0.0,0.0,0.0 &                
],[3,6])/sqrt(3.0)
!-----------------------------------------------------------------------|
!                                                                       |
!                               read file                               |
!                                                                       |
!-----------------------------------------------------------------------|
open(unit=100,file='dump.atom',status='old',action='read')
open(unit=300,file='o.atom',status='old',action='read')
open(unit=500,file='data/direction.data',status='replace',action='write')
open(unit=550,file='data/direction_explict.data',status='replace',action='write')
open(unit=600,file='data/s_d.data',status='replace',action='write')
!-----------------------------------------------------------------------|
!                                                                       |
!                 confirm atom number and crystal constant              |
!                                                                       |
!-----------------------------------------------------------------------| 
do i=1,3
        read(100,*)
end do
read(100,*) atom_number(1)
read(100,*)
read(100,*) x_b(1),x_b(2)
read(100,*) x_b(3),x_b(4)
read(100,*) x_b(5),x_b(6)
a=x_b(2)-x_b(1)
b=x_b(4)-x_b(3)
c=x_b(6)-x_b(5)
rewind(100)
!-----------------------------------------------------------------------|
!                                                                       |
!                               confirm step                            |
!                                                                       |
!-----------------------------------------------------------------------| 
m=1
do while (.true.)
        do i=1,9
                read(100,*,iostat=alive)
                if (alive/=0) then
                        goto 55
                end if
        end do
        
        do i=1,atom_number(1)
                read(100,*) 
        end do
        step=step+1
 end do
!-----------------------------------------------------------------------|
!                                                                       |
!                               read Li atom                            |
!                                                                       |
!-----------------------------------------------------------------------| 
55 continue
rewind(100)
allocate(R(atom_number(1),step))
do i=1,step
        do j=1,9
                read(100,*)
        end do

        do j=1,atom_number(1)
                read(100,*) n,R(n,i)%atom_type,R(n,i)%x(1),R(n,i)%x(2),R(n,i)%x(3)
                R(n,i)%x(1)=R(n,i)%x(1)-x_b(1)
                R(n,i)%x(2)=R(n,i)%x(2)-x_b(3)
                R(n,i)%x(3)=R(n,i)%x(3)-x_b(5)
                if (i==1) then
                        atom_number(1+R(n,i)%atom_type)=atom_number(1+R(n,i)%atom_type)+1
                end if
        end do
end do

allocate(if_diff(atom_number(2),step-1))
allocate(dx(atom_number(2),step-1))
allocate(dy(atom_number(2),step-1))
allocate(dz(atom_number(2),step-1))
allocate(co(atom_number(2),step))
allocate(LiO(atom_number(1)))
allocate(direction(atom_number(2),step-1))
allocate(d_new(3,atom_number(2),step-1))
d_new(:,:,:)=0.0
direction(:,:)=0
co(:,:)=0
if_diff(:,:)=0
statistic_direction(:)=0
!-----------------------------------------------------------------------|
!                                                                       |
!                                O atom                                 |
!                                                                       |
!-----------------------------------------------------------------------|
do i=1,atom_number(5)
        read(300,*) n,LiO(n)%atom_type,LiO(n)%x(1),LiO(n)%x(2),LiO(n)%x(3)
        LiO(n)%atom_type=LiO(n)%atom_type
        LiO(n)%x(1)=LiO(n)%x(1)
        LiO(n)%x(2)=LiO(n)%x(2)
        LiO(n)%x(3)=LiO(n)%x(3)
end do
!-----------------------------------------------------------------------|
!                                                                       |
!                               Diffusion Ensurement                    |
!                                                                       |
!-----------------------------------------------------------------------|
!$omp parallel do schedule(static) private(j,k,m,dr,lix,liy,liz,ox,oy,oz,o_r,time)
do i=1,atom_number(2)
        do j=1,step-1
                dx(i,j)=R(i,j+1)%x(1)-R(i,j)%x(1)
                dy(i,j)=R(i,j+1)%x(2)-R(i,j)%x(2)
                dz(i,j)=R(i,j+1)%x(3)-R(i,j)%x(3)

                call minimum_image(dx(i,j),a)
                call minimum_image(dy(i,j),b)
                call minimum_image(dz(i,j),c)
                if (sqrt(dx(i,j)**2+dy(i,j)**2+dz(i,j)**2)>=cutoff_r) then
                        do k=atom_number(2)+1,atom_number(1)
                                ox=(R(i,j+1)%x(1)-LiO(k)%x(1))
                                oy=(R(i,j+1)%x(2)-LiO(k)%x(2))
                                oz=(R(i,j+1)%x(3)-LiO(k)%x(3))
                                call minimum_image(ox,a)
                                call minimum_image(oy,b)
                                call minimum_image(oz,c)

                                o_r=sqrt(ox**2+oy**2+oz**2)

                                if (o_r<=cutoff_r2) then
                                        co(i,j+1)=co(i,j+1)+1
                                end if
                        end do

                        if (co(i,j+1)==6) then
                                if_diff(i,j)=1
                        else 
                                if_diff(i,j)=2
                        end if
                end if
                write(400,*) i,j,if_diff(i,j)
                
        end do
        if (mod(i,100)==0) then
                call cpu_time(time)
                write(*,*) i,' is done ',' cpu time: ',time
        end if
end do
!$omp end parallel do
!-----------------------------------------------------------------------|
!                                                                       |
!                               Diffusion Direction                     |
!                                                                       |
!-----------------------------------------------------------------------|
!!$omp parallel do schedule(static) private(j,k,l,dir,vec,cos_angle,time)
do i=1,atom_number(2)
step_loop:do j=1,step-1
                if (if_diff(i,j)/=0) then
                        vec(1)=dx(i,j)
                        vec(2)=dy(i,j)
                        vec(3)=dz(i,j)
                        m=0
                        do k=1,3
                                do l=1,6
                                        m=m+1
                                        dir=group_dirs(:,(k-1)*6+l)
                                        cos_angle(m)=cosin(vec,dir)
                                end do
                        end do

                        min_angle=9999.0
                        do k=1,18
                                if (abs(abs(cos_angle(k))-1.0)<=min_angle .and. k/=18) then
                                        min_angle=abs(abs(cos_angle(k))-1.0)
                                        min_index=k
                                end if
                        end do
                        if (min_index>=1 .and. min_index<=6) then
                                k=1
                        else if (min_index>=7 .and. min_index<=12) then
                                k=2
                        else if (min_index>=13 .and. min_index<=18) then
                                k=3
                        end if
                                direction(i,j)=k
                                statistic_direction(k)=statistic_direction(k)+1
                                write(500,*) i,j,direction(i,j)
                else
                        cycle
                end if
        end do step_loop
end do
!!$omp end parallel do
!-----------------------------------------------------------------------|
!                                                                       |
!                               Diffusion Direction 2                   |
!                                                                       |
!-----------------------------------------------------------------------|
do i=1,atom_number(2)
        do j=1,step-1
                if (direction(i,j)==3 .and. if_diff(i,j)==1) then
                        max_step=min(3,step-1-j)
                        do k=1,max_step
                                if (direction(i,j+k)/=3) then
                                        exit
                                else
                                        d_new(1,i,j)=R(i,j+k+1)%x(1)-R(i,j)%x(1)
                                        d_new(2,i,j)=R(i,j+k+1)%x(2)-R(i,j)%x(2)
                                        d_new(3,i,j)=R(i,j+k+1)%x(3)-R(i,j)%x(3)
                                        call minimum_image(d_new(1,i,j),a)
                                        call minimum_image(d_new(2,i,j),b)
                                        call minimum_image(d_new(3,i,j),c)
                                        exit
                                end if
                        end do
                        vec(1)=d_new(1,i,j)
                        vec(2)=d_new(2,i,j)
                        vec(3)=d_new(3,i,j)
                        m=0
                        do k=1,3
                                do l=1,6
                                        m=m+1
                                        dir=group_dirs(:,(k-1)*6+l)
                                        cos_angle(m)=cosin(vec,dir)
                                end do
                        end do

                        min_angle=9999.0
                        do k=1,18
                                if (abs(abs(cos_angle(k))-1.0)<=min_angle .and. k/=18) then
                                        min_angle=abs(abs(cos_angle(k))-1.0)
                                        min_index=k
                                end if
                        end do

                        if (min_index>=1 .and. min_index<=6) then
                                k=4
                        else if (min_index>=7 .and. min_index<=12) then
                                k=5
                        else if (min_index>=13 .and. min_index<=18) then
                                k=6
                        end if

                        direction(i,j)=k
                        statistic_direction(k)=statistic_direction(k)+1
                        write(550,*) i,j,direction(i,j)
                end if
        end do
end do

write(600,*) (statistic_direction(i)/step,i=1,6)
stop
end program