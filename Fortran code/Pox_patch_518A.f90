program MainABM
  implicit none

  ! We assume dt = 1 (day) in this model, so dt is usually suppressed.

  integer, parameter :: M = 592, initial_patch = 518, L=1000, NR = 1, step = 1 !steps/yr
  integer, parameter :: parameters = 1, total_steps = L/step
  !real(4), parameter :: tau = .1;
  real(4) :: alpha, beta, delta, gamma, rho_bar, row_sum, lambda, epsilon, epsilonvec(4), lambdavec(2)
  real(4) :: mu_bar(M),A(M),Population(M),rhovec(2),p_bar_r(M),p_bar_v(M),mu_vec(M),Coordinates(M,2)
  real(4) :: PatchProbs(M)
  integer :: ii, jj, rr, kk, j, t, f, cols, rows, avec(3), bvec(2), ios, E0(M)
  integer :: aloop, bloop, rloop, lambdaloop, epsilonloop
  integer, dimension(M,M) :: S, E, I, R, D
  real(4), dimension(M,M) :: Q, Q_new
  real(4), dimension(M) :: H_trade, H_ceremony, T_trade, T_ceremony
  real(4) :: ss, f_trade, f_ceremony, nu
  integer, dimension(NR,L,parameters) :: Total_I_Tally
  integer, dimension(NR,M) :: IPP
  integer, dimension(NR,L) :: I_Tally
  integer, dimension(parameters,NR) :: Epi_Times, Total_E_All
  integer, dimension(NR,M,parameters) :: Total_IPP
  integer, dimension(NR) :: Total_unique_E,TT
  integer :: start_clock, current_clock, clock_rate
  real :: elapsed_time_seconds
  character(len=20) :: filename
  call random_seed()

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Load Q matrix ...probability of movement from patch to patch ~ M x M

  rows = M
  cols = M

  open(unit=8,file='Q.txt',status='old',action='read', iostat=ios)

  do ii = 1, rows
    read(8,*) Q(ii,1:cols)
  end do

  close(8)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Load Patch Probability vector
  open(unit=9, file = 'PatchProbs.txt',status='old',action='read',iostat=ios)

  read(9,*) PatchProbs

  close(9)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Load A vector....areas of each patch 1 x M

  open(unit=10, file='Areas.txt',status='old', action='read', iostat=ios)

  read(10,*) A
  close(10)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Load Population vector....1 x M

  open(unit=11, file='Populations.txt', status='old', action='read', iostat=ios)

  read(11,*) Population
  close(11)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Load Coordinate matrix ...location of each patch "center"...Mx2

  rows = M
  cols = 2

  open(unit=12,file='Coordinates.txt',status='old',action='read', iostat=ios)

  do ii = 1, rows
    read(12,*) Coordinates(ii,1:cols)
  end do

  close(12)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Initialize parameters
  E0 = 0
  E0(initial_patch) = 1

  avec = [3, 8, 12]
  bvec = [1, 2]
  rhovec = [1.0/5.0, 1.0/14.0]                              !~1/days 
  alpha = 1.0 - exp(-1.0/12.0)                              !~12 days
  gamma = 1.0 - exp(-1.0/18.0)                              !~18 days
  epsilonvec = [0.0, 0.01, 0.05, 0.1]
  lambdavec = [0.08, 0.1]
  beta = 0.334
  ss = 6.0
  f_trade = 0.25
  f_ceremony = 0.8
  nu = 0.75

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Simulation loop for different parameter sets (bvec and avec)

  kk = 1

  rloop = 1
  aloop = 2
  bloop = 1
  epsilonloop = 1
  lambdaloop = 1

  !do rloop = 1, size(rhovec)
  rho_bar = 1 - EXP(-rhovec(rloop))
  !do bloop = 1, size(bvec)  ! Change to loop over other parameter sets if needed
    !do aloop = 1, size(avec)  ! Adjust for other parameter combinations

H_trade    = f_trade    * Population / ss
H_ceremony = f_ceremony * Population / ss

T_trade    = avec(aloop) * H_trade**nu
T_ceremony = bvec(bloop) * H_ceremony**nu

mu_vec = (1.0 / 365.0) * (T_trade / Population) + (1.0 / 730.0) * (T_ceremony / Population)

mu_bar = 1 - EXP(-mu_vec)

!do epsilonloop = 1,size(epsilonvec)
  !do lambdaloop = 1,size(lambdavec)

epsilon = epsilonvec(epsilonloop)
lambda  = lambdavec(lambdaloop)

Q_new = DistanceDecay(Q, Coordinates, lambda, epsilon, M, PatchProbs)


open(unit=99, file="Q_row_209.txt")
write(99,*) Q_new(209, :)
close(99)


      call ABM_sub(S,E,I,R,D,NR,E0,A,M,L,mu_bar,rho_bar,alpha,gamma,beta,Q_new,Population,f,&
    p_bar_r,p_bar_v,step,total_steps,IPP,Total_unique_E,TT,I_Tally)

      
      Total_E_All(kk,:) = Total_unique_E
      Total_I_Tally(:,:,kk) = I_Tally
      Epi_Times(kk,:)   = TT
      Total_IPP(:,:,kk) = IPP

      kk = kk + 1

    !end do
  !end do

  !end do

    !end do !lambda loop
  !end do ! epsilon loop

  call save_4D_arrays(M,NR,total_steps,parameters,Total_IPP,Epi_Times,Total_I_Tally)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  contains
  subroutine ABM_sub(S,E,I,R,D,NR,E0,A,M,L,mu_bar,rho_bar,alpha,gamma,beta,Q_new,Population,f,&
    p_bar_r,p_bar_v,step,total_steps,IPP,Total_unique_E,TT,I_Tally)

  implicit none

    real(4), parameter :: dd = (1/38)*(1/365)    !~38 years expectancy / daily rate
    integer, parameter :: dt = 1                 ! time step
    real(4), intent(in) :: rho_bar, alpha, gamma, beta
    integer:: counter, TE
    integer, intent(in) :: step,total_steps
    integer, dimension(M), intent(in) :: E0
    real(4), dimension(M,M), intent(in) :: Q_new
    real(4), dimension(M), intent(in) :: A, mu_bar,p_bar_r,p_bar_v
    real(4), dimension(M) :: mu
    real(4), dimension(M,M) :: MuRho, alpha_mat, gamma_mat, delta_mat
    integer, intent(in) :: NR, M, L, f
    integer :: j, k, Total_Time,jj
    real(4) :: pi, PC,PD, rho ! Probability of Conception, Death
    real(4), dimension(M),     intent(in)     :: Population
    integer, dimension(NR)                    :: TimeEnd,TotalUniqueE
    integer, dimension(NR,M), intent(out)    :: IPP
    integer, dimension(M,M), intent(inout)    :: S, E, I, R, D
    integer, dimension(total_steps,M)         :: E_Mid
    integer, dimension(NR,total_steps), intent(inout) :: I_Tally
    integer, dimension(total_steps) :: TI_vec
    integer, dimension(NR), intent(inout) :: Total_unique_E,TT
    integer :: TR,TI,TD,T, A1, A2, A3, A4

    PD = 1 - EXP(-dd)              !Probability of a death -- dd is death rate
    PC = 0.0001098            !Probability of conception -- 0.5 female|1/30 f window|1/365
    alpha_mat = alpha
    gamma_mat = gamma
    delta_mat = delta
    pi = acos(-1.0)
    T = 90                              !Period for dispersal and transmission
    A1 = 0.0
    A2 = 0.0
    A3 = 0.0
    A4 = 0.0

    ! Loop for each realization (NR)
    do j = 1, NR

      ! Initialize Susceptible, Exposed, Infected, Recovered matrices
      S = 0
      E = 0
      I = 0
      R = 0
      D = 0
      E_Mid = 0
      TE = 0

      ! Initialize initial population (Seed) distributions
      S = diagmatrix(real(Population, kind=4)) !creates diagonal matrices....
      E = diagmatrix(real(E0,4))

      TI = sum(I) 
      TR = sum(R)
      counter = 1

      do k = 1, L

        mu = max(mu_bar*(1+A1*sin((2*pi*k)/T)),0.0)
        rho = max(rho_bar*(1+A2*sin((2*pi*k)/T)),0.0)
        MuRho = move_matrix(M,mu,rho)

        call Move(S, E, I, R, MuRho, Q_new)

        call Spreading(S, E, I, R, A, dt, M, TE)

        call RecoverDeath(S,E,I,R,alpha_mat,gamma_mat)

        print *, 'day =', k, 'I =',sum(I)

        if (mod(k,step)==0) then
              TI = sum(I)

        E_Mid(counter,:) = sum(E, dim=1)
        TI_vec(counter) = TI
        counter = counter + 1

        end if

        if (sum(I)==0 .AND. sum(E)==0) then
          Total_Time = k
          exit 
        end if

      end do

        Total_Time = k
        IPP(j,:) = sum(E_Mid, dim=1)
        Total_unique_E(j) = TE
        I_Tally(j,:) = TI_vec
        TT(j) = Total_Time
    end do
  end subroutine ABM_sub
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine Move(S, E, I, R, MuRho, Q)
  implicit none
    real(4), dimension(M,M), intent(in)    :: MuRho, Q
    integer, dimension(M,M), intent(inout) :: S, E, I, R
    integer, dimension(M,M)                :: S_movers, E_movers, I_movers, R_movers
    integer, dimension(M,M)                :: S_visit, E_visit, I_visit, R_visit
    integer, dimension(M,M)                :: S_return, E_return, I_return, R_return
    integer, dimension(M)                  :: vs, ve, vi, vr
    real(4), dimension(M)                  :: SDS,SDE,SDR
    !real(4):: factor

    ! Movers for each state (S, E, I, R) - binomial distribution for each patch
     S_movers = binornd2(S, MuRho)
     E_movers = binornd2(E, MuRho)
     R_movers = binornd2(R, MuRho)

    ! Reshuffling the population into new destinations (visitors and returners)
     vs = diag(S_movers)
     ve = diag(E_movers)
     vr = diag(R_movers)

    ! Generate multinomial distributions for the visitors' destination
     S_visit =  mnrnd(vs, Q)
     E_visit =  mnrnd(ve, Q)
     R_visit =  mnrnd(vr, Q)
     
    ! Returners - determine where they come back from
    SDS = sum(diagremove(real(S_movers,4)),dim=2)
    SDE = sum(diagremove(real(E_movers,4)),dim=2)
    SDR = sum(diagremove(real(R_movers,4)),dim=2)
    S_return = diagmatrix(SDS)
    E_return = diagmatrix(SDE)
    R_return = diagmatrix(SDR)

    ! Update population states with new destinations
     S = S - S_movers + S_visit + S_return
     E = E - E_movers + E_visit + E_return
     R = R - R_movers + R_visit + R_return

    S = max(S, 0)  ! Avoid rounding loss
    E = max(E, 0)
    I = max(I, 0)
    R = max(R, 0)

  end subroutine Move
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine RecoverDeath(S, E, I, R, alpha_mat, gamma_mat)

  implicit none

    integer, dimension(M,M), intent(inout) :: S, E, I, R
    integer, dimension(M,M) :: Inew_E, End_I, Dnew_I
    real(4), dimension(M,M), intent(in) :: alpha_mat, gamma_mat

    Inew_E = min(binornd2(E, alpha_mat), E)

    I = I + Inew_E
    E = E - Inew_E

    End_I  = min(binornd2(I, gamma_mat), I)

    R = R + End_I
    I = I - End_I

    E = max(E, 0)
    I = max(I, 0)
    R = max(R, 0)

  end subroutine RecoverDeath
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine Spreading(S, E, I, R, A, dt, M, TE)
  implicit none
    integer, dimension(M,M), intent(inout) :: S, E, I, R
    integer, intent(inout) :: TE
    integer, dimension(M,M) :: New_E_resident, New_E_visitor, New_E
    real(4) :: x, ii, tau
    integer, intent(in) :: M, dt
    real(4), dimension(M), intent(in) :: A
    integer, dimension(M) :: ones, zeros
    integer, dimension(M) :: xi_r, xi_v, eta, N_r, N_v, sumI
    real(4), dimension(M,M) :: P,pp
    integer, dimension(M,M) :: DMS
    ones = 1
    zeros = 0
     
    ! Calculate SumI, the sum of infectious individuals per patch
    SumI = sum(I, dim=2)

    ii = 0;

    do while (ii < dt)
    P = beta * S * spread(SumI, dim=2, ncopies=M) / spread(A, dim=2, ncopies=M)
    tau = adaptive_tau(S,P,dt-ii)
    New_E = poisson_rand(P*tau)
    New_E = min(New_E,S)

    S = S - New_E
    E = E + New_E

    ii = ii + tau
    TE = TE + sum(New_E)

    S = max(S,0)

    end do

  end subroutine Spreading
!%%%%%%%%%%%%%%%%%%%%%%_SAVE_OUTPUT_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine save_4D_arrays(M,NR,total_steps,parameters,Total_IPP,Epi_Times,Total_I_Tally)

  implicit none

  integer, intent(in) :: M, NR,parameters,total_steps
  character(len=30) :: filename1, filename2, filename3, filename4, filename5, filename6, filename7,filename8
  integer, dimension(parameters,NR) :: Epi_Times
  integer, dimension(NR,M) :: IPP
  integer, dimension(NR,L,parameters) :: Total_I_Tally
  integer, dimension(NR,M,parameters) :: Total_IPP
  integer :: unit
  
  filename1 = "518A_I_Tally.bin"
  filename2 = "518A_Unique_E.bin"
  filename3 = "518A_Epi_Lengths.bin"
  filename4 = "518A_Total_IPP.bin"

  open(unit=1, file=filename1, status="replace", form="unformatted", action="write", access="stream")
  write(1) Total_I_Tally
  close(1)

  open(unit=2, file=filename2, status="replace", form="unformatted", action="write", access="stream")
  write(2) Total_E_All
  close(2)

  open(unit=3, file=filename3, status="replace", form="unformatted", action="write", access="stream")
  write(3) Epi_Times
  close(3)

  open(unit=4, file=filename4, status="replace", form="unformatted", action="write", access="stream")
  write(4) Total_IPP
  close(4)

  end subroutine save_4D_arrays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mnrnd(n, p) result(result)
  implicit none
  integer, dimension(:), intent(in) :: n  ! Vector of number of trials
  real(4), dimension(:,:), intent(in) :: p  ! Matrix of probabilities (categories x trials)
  integer, dimension(size(p,1), size(n)) :: result  ! Resulting counts (categories x trials)

  integer :: i, j, k
  real(4) :: u, cumulative_prob
  real(4), dimension(size(p,1)) :: cumulative_probs

  ! Initialize result array
  result = 0.0

  ! Loop over trials
  do k = 1, size(n)
    ! Calculate cumulative probabilities for each trial
    cumulative_probs = 0.0
    cumulative_probs(1) = p(k,1)
    do j = 2, size(p,2)
      cumulative_probs(j) = cumulative_probs(j-1) + p(k,j)
    end do

    ! Generate n(k) samples for this trial
    do i = 1, n(k)
       ! Generate a random number between 0 and 1
       call random_number(u)
       !u = u * sum(p(k,:))
       ! Find which category the random value falls into
       do j = 1, size(p,1)
          if (u <= cumulative_probs(j)) then
             result(j, k) = result(j, k) + 1
             exit
          end if
       end do
    end do
  end do

end function mnrnd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diagmatrix(vector) result(matrix) !Constructs a diagonal matrix from vector input
  implicit none
  real(4), dimension(:), intent(in) :: vector   ! Input vector
  real(4), dimension(:,:), allocatable :: matrix ! Output matrix
  integer :: i, n

  n = size(vector)    ! Get the size of the input vector
  if (n>0) then
  allocate(matrix(n, n))  ! Allocate the output matrix as an n x n matrix
  !allocate(vector(n))

  else 
    print *, 'ERROR: Invalid array size n=', n
  end if

  ! Initialize the matrix to zero
  matrix = 0.0

  ! Set the diagonal elements
  do i = 1, n
     matrix(i, i) = vector(i)
  end do
end function diagmatrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function matrep(v) result(Matrix)
    implicit none
    real(4), dimension(:), intent(in) :: v  ! Input vector
    real(4), dimension(:,:), allocatable :: Matrix  ! Output matrix
    integer :: i, n
    ! Get the number of elements in the vector
    n = size(v)
    if (n>0) then
    allocate(Matrix(n,n))
    else 
      print *, 'ERROR: Invalid array size n=',n
    end if

    ! Fill each row of M with the vector v
    Matrix = spread(v, dim=1, ncopies=n) ! Repeat the vector n times and reshape it

  end function matrep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diagremove(matrix) result(output_matrix) ! Returns the matrix with diagonal elements set to zero
  implicit none
  real(4), dimension(:,:), intent(in) :: matrix  ! Input matrix
  real(4), dimension(size(matrix, 1), size(matrix, 2)) :: output_matrix
  integer :: i, n, m

  ! Get the number of rows and columns
  n = size(matrix, 1)  ! Number of rows
  m = size(matrix, 2)  ! Number of columns

  ! Initialize the output matrix with the input matrix
  output_matrix = matrix

  ! Loop through the diagonal elements and set them to zero
  do i = 1, min(n, m)  ! Only loop through the diagonal (min of rows and columns)
     output_matrix(i, i) = 0.0
  end do
end function diagremove
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine sdiag(matrix, sum_vector)
  implicit none
  integer, intent(in) :: matrix(:,:)
  integer, dimension(size(matrix,2)), intent(out) :: sum_vector
  integer :: i, j
  integer :: rows, cols

  ! Get the number of rows and columns
  rows = size(matrix, 1)
  cols = size(matrix, 2)

  ! Initialize the sum_vector to zero
  sum_vector = 0

  ! Loop over the columns of the matrix
  do j = 1, cols
     ! Loop over the rows
     do i = 1, rows
        ! Add to the sum for the current column, excluding the diagonal element
        if (i /= j) then
           sum_vector(j) = sum_vector(j) + matrix(i, j)
        end if
     end do
  end do
end subroutine sdiag
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diag(matrix) !This extracts the diagonal elements of a 2D array.
  implicit none
  integer, dimension(:,:), intent(in) :: matrix  ! Input matrix
  integer :: i
  integer :: rows, cols
  integer, dimension(:), allocatable :: diag

  ! Get the number of rows and columns of the matrix
  rows = size(matrix, 1)
  cols = size(matrix, 2)
  if ((rows>0 ) .and. (cols>0)) then
    allocate(diag(min(rows,cols)))
  else 
    print *, "ERROR: Invalid size"
  end if

  ! The diagonal will have the minimum of rows and columns
  if (size(diag) < min(rows, cols)) then
     print *, 'Error: diagonal array size is smaller than the matrix diagonal size, you dumbass!'
     stop
  end if

  ! Extract the diagonal elements
  do i = 1, min(rows, cols)
     diag(i) = matrix(i, i)
  end do

end function diag
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binornd2(trials, prob_success) result(out_matrix)
  implicit none
  integer, dimension(:,:), intent(in) :: trials  ! Input matrix for trials
  real(4), dimension(:,:), intent(in) :: prob_success  ! Input matrix for success probabilities
  integer, dimension(size(trials, 1), size(trials, 2)) :: out_matrix  ! Output matrix for binomial random numbers
  real(4) :: random_value
  integer :: i, j

  ! Initialize output matrix to zero
  out_matrix = 0

  ! Loop through matrix elements and generate binomial values
  do i = 1, size(trials, 1)
    do j = 1, size(trials, 2)

      ! Validate trials
      if (trials(i, j) < 0) then
          print *, "Error: Negative trials at (", i, ",", j, "):", trials(i, j)
          cycle
      end if

      ! Validate probability
      if (prob_success(i, j) < 0.0 .or. prob_success(i, j) > 1.0) then
          print *, "Error: Invalid probability at (", i, ",", j, "):", prob_success(i, j)
          cycle
      end if

      ! Call binomial generator
      call generate_binomial(trials(i, j), prob_success(i, j), out_matrix(i, j))

      ! Debug output
      !print *, "Generated binomial value at (", i, ",", j, "):", out_matrix(i, j)
      !print *, "Probability:", prob_success(i,j)

      ! Ensure valid output
      if (out_matrix(i, j) < 0) then
          print *, "Warning: Negative binomial output at (", i, ",", j, "):", out_matrix(i, j)
      end if

    end do
  end do

end function binornd2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function poisson_rand(lambda) result(poisson_matrix)
    use iso_fortran_env, only: real32, int32
    implicit none
    real(real32), dimension(:,:), intent(in) :: lambda  ! Matrix of lambda values
    integer(int32), dimension(size(lambda,1), size(lambda,2)) :: poisson_matrix
    integer(int32) :: i, j, k, rows, cols
    real(real32) :: u, L, p
    real(real32) :: sum_prob
    integer(int32) :: lambda_int

    ! Get matrix dimensions
    rows = size(lambda, 1)
    cols = size(lambda, 2)

    ! Loop through each element of the matrix
    do i = 1, rows
        do j = 1, cols
            lambda_int = nint(lambda(i,j))  ! Convert lambda to integer for computation
            L = exp(-lambda(i,j))          ! L = exp(-lambda)
            p = 1.0
            k = 0

            ! Inverse Transform Sampling to generate Poisson-distributed random variable
            call random_number(u)  ! Generate a uniform random number between 0 and 1
            sum_prob = u

            ! Sum until the cumulative probability exceeds a uniform random value
            do while (sum_prob >= L)
                k = k + 1
                call random_number(u)  ! Generate another uniform random number
                sum_prob = sum_prob * u
            end do

            poisson_matrix(i,j) = k  ! Assign the generated Poisson random variable
        end do
    end do
end function poisson_rand

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine generate_binomial(n, p, result)
  implicit none
  integer, intent(in) :: n
  real(4), intent(in) :: p
  integer, intent(out) :: result
  integer :: i
  real(4) :: rand_val

  if (n <= 0) then
    result = 0
    return
  end if

  result = 0  ! Initialize count of successes

  do i = 1, n
    call random_number(rand_val)
    if (rand_val < p) then
      result = result + 1  ! Count a success
    end if
  end do

end subroutine generate_binomial

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function randint(n, lower_bound, upper_bound) result(randintegers)
  implicit none
  integer, intent(in) :: n
  integer, dimension(:), intent(in) :: lower_bound, upper_bound
  integer, dimension(n) :: randintegers  ! This is the output variable
  real(4) :: random_value
  integer :: i

  ! Check for valid range
  !do i = 1, n
  !   if (lower_bound(i) > upper_bound(i)) then
  !     print *, 'Error: lower bound is greater than upper bound for index ', i
  !      stop
  !   end if
  !end do

  ! Generate n random integers in the specified range
  do i = 1, n
     call random_number(random_value)
     randintegers(i) = lower_bound(i) + int(random_value * max(0,(upper_bound(i) - lower_bound(i))))
  end do

end function randint

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function move_matrix(size, diagvalue, non_diagvalue) result(A) !nxn matrix with diagional and non diag 
        implicit none
        integer, intent(in) :: size
        real(4), dimension(:) :: diagvalue
        real(4), intent(in) :: non_diagvalue
        real(4), dimension(size,size):: A 
        integer :: i, j

        !print *, "INSIDE FUNCTION, diagvalue=", diagvalue

        ! Fill the matrix
        do i = 1, size
            do j = 1, size
                if (i == j) then
                    A(i, j) = diagvalue(i)  
                else
                    A(i, j) = non_diagvalue  
                end if
            end do
        end do
    end function move_matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function adaptive_tau(S_matrix, P_matrix, remaining_time) result(tau)
    use ieee_arithmetic
    implicit none
    real(4), parameter :: epsilon = 0.03
    integer, dimension(:,:), intent(in) :: S_matrix 
    real(4), dimension(:,:), intent(in) :: P_matrix ! probabilities
    real(4), intent(in) :: remaining_time
    real(4) :: tau1, tau2, tau
    real(4), dimension(size(S_matrix,1), size(S_matrix,2)) :: S_real

    ! Convert S_matrix to real to avoid integer division issues
    S_real = real(S_matrix, 4)

    ! Ensure element-wise division
    tau1 = minval(epsilon * (S_real / P_matrix))

    ! Ensure correct matrix multiplication
    tau2 = minval((epsilon**2) * (matmul(S_real, S_real) / P_matrix))

    ! Replace NaN or Inf values with positive infinity
    if (.not. ieee_is_finite(tau1)) tau1 = ieee_value(1.0, ieee_positive_inf)
    if (.not. ieee_is_finite(tau2)) tau2 = ieee_value(1.0, ieee_positive_inf)

    ! Compute the final tau value
    tau = min(tau1, tau2)
    tau = min(tau, remaining_time)
    tau = max(tau, 1e-6)  ! Ensure tau is not too small

end function adaptive_tau
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DistanceDecay(Q, Coordinates, lambda, epsilon, M, Pref) result(Qnew)
  implicit none
  integer, intent(in) :: M
  real, intent(in) :: Q(M,M), Coordinates(M,2), lambda, epsilon
  real, intent(in) :: Pref(M)
  real :: Qnew(M,M)

  real :: D(M,M), lat(M), lon(M)
  real :: long_dist_probs(M), sum_probs
  real :: dlat, dlon, a, c, R_earth, pi
  integer :: i, j

  pi = 3.141592653589793
  R_earth = 6371.0

  lat = Coordinates(:,2)  ! Latitude = column 2
  lon = Coordinates(:,1)  ! Longitude = column 1

  ! Compute pairwise Haversine distances
  do i = 1, M
    do j = 1, M
      dlat = (lat(j) - lat(i)) * pi / 180.0
      dlon = (lon(j) - lon(i)) * pi / 180.0
      a = sin(dlat/2.0)**2 + cos(lat(i)*pi/180.0) * cos(lat(j)*pi/180.0) * sin(dlon/2.0)**2
      c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
      D(i,j) = R_earth * c
    end do
  end do

  ! Construct distance-decay + preference movement probabilities
  do i = 1, M
    Qnew(i,:) = Q(i,:)  ! Start from original matrix

    do j = 1, M
      if (Q(i,j) > 0.0 .or. i == j) then
        long_dist_probs(j) = 0.0
      else
        long_dist_probs(j) = Pref(j) * exp(-lambda * D(i,j))
      end if
    end do

    sum_probs = sum(long_dist_probs)
    if (sum_probs > 0.0) then
      long_dist_probs = epsilon * (long_dist_probs / sum_probs)
    end if

    Qnew(i,:) = Qnew(i,:) + long_dist_probs

    ! Normalize full row
    if (sum(Qnew(i,:)) > 0.0) then
      Qnew(i,:) = Qnew(i,:) / sum(Qnew(i,:))
    end if
  end do

end function DistanceDecay
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end program MainABM
