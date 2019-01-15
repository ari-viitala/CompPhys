program rd
  !
  ! Modified random deposition simulation, homework
  !
  ! Compile as: gfortran -o RD_hw RD_hw.f90
  !
  implicit none ! This forces us to define all variables
  integer, allocatable ::  h(:) ! height
  real(kind(1.d0)), allocatable ::  wsq(:) ! roughness squared
  integer :: i, j, L, Nh, Nr, ii, jj
  character(len=40) :: w_file
  logical :: do_modified=.true.
  integer :: idate(8), iseed

  ! Initialization of rng by time:
  call date_and_time(VALUES=idate)
  write(*,*) 'date and time is ', idate
  iseed=idate(7)*1000+idate(8)
  write(*,*) 'seed is ', iseed
  call random_seed(iseed)
  
  write(*,*) 'Do modified RD? (false makes pure RD, true the modified one)'
  read(*,*) do_modified
  write(*,*) 'Number of sites (horizontal size of cluster) [try 100]'
  read(*,*) L
  !Note that we allocate two extra array elements. The actual heights are kept in h(1),h(2),...,h(L),
  !and h(0) and h(L+1) are always updated to contain the same value as h(L) and h(1) respectively.
  !This makes it easier to implement periodic boundary conditions, as the neighbours of any
  !site i (including i=1 and i=L) can be accessed as h(i-1) and h(i+1).
  allocate(h(0:L+1)) 
  write(*,*) 'How many drops (per column)? [try 10000]'
  read(*,*) Nh
  write(*,*) 'How many runs? [try 100]'
  read(*,*) Nr
  allocate(wsq(Nh))

  !Open a file for outputting the roughness values.
  if (do_modified) then 
     write(w_file,'(a,i0,a)') 'modified_w_average', L, '.txt'
  else
     write(w_file,'(a,i0,a)') 'rd_w_average', L, '.txt'     
  end if
  open(123,file=w_file)
  write(*,*) 'Data goes to file ', w_file

  wsq=0.d0

  !Loop over the number of runs.
  do jj=1, Nr

     write(*,*) 'Run ', jj, '/', Nr
     h=0  ! This contains heights.

     !Do Nh time steps.
     do ii=1, Nh

        !And for each time step do L drop trials.
        do j=1, L
           i=random_site(L)
           if (do_modified) then

              ! TODO: Increase h(i) if the neighbouring columns are not
              ! lower that site h(i)
              
              if ((h(i - 1) >= h(i)) .and. (h(i + 1) >= h(i))) then
                h(i) = h(i) + 1
              end if
              
              if (i==1) h(L+1)=h(1) ! These take care of 
              if (i==L) h(0)=h(L)   ! the periodic boundary conditions.
           else
              h(i)=h(i)+1 ! This is the rule for pure random deposition.
           end if
        end do

        wsq(ii)=wsq(ii)+rough2(h(1:L))
     end do

  end do
  wsq=wsq/real(Nr,kind(1.d0))
  do ii=1, Nh
     !Write the roughness values to a file. Take square root to output w and not w^2.
     write(123,*) ii, sqrt(wsq(ii))
  end do
  close(123)
  write(*,*) 'Data in ', w_file

contains

  !A function to get a random number in the interval [1,L]
  function random_site(L)
    integer :: random_site, L
    real(kind(1.d0)) :: r
    call random_number(r)
    random_site=min(int(r*L)+1,L)
  end function random_site

  !This function calculates the roughness value w for a given
  !height array h.
  function rough2(h)
    real(kind(1.d0)) :: rough2, eh, eh2
    integer :: h(:), L
    real(kind(1.d0)), dimension(size(h)) :: nh
    L=size(h)
    nh=h-minval(h)
    eh=sum(nh)/real(L,kind(1.d0))
    eh2=sum(nh**2)/real(L,kind(1.d0))
    rough2=eh2-eh**2
  end function rough2
end program rd
