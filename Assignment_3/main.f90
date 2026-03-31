! ==========================================
! Main program for 2D Co-rotational Truss FEA
! ==========================================
program main
    use kinds_mod
    use material_mod
    use solver_mod
    implicit none

    type(LinearElasticMaterial) :: my_mat
    type(CorotationalTrussFEA)  :: model
    
    ! Global Mesh Parameters
    integer, parameter :: n_bays = 6
    integer, parameter :: n_elem = (5 * n_bays) + 1  ! 5 elements per bay: bottom, top, left, right, diagonal
    integer, parameter :: n_node = 2 * n_bays + 2
    integer, parameter :: ndof   = n_node * 2
    
    ! Material and Geometric Properties
    real(wp), parameter :: L_bay   = 40.0_wp                        ! mm
    real(wp), parameter :: L_total = n_bays * L_bay                 ! mm
    real(wp), parameter :: H_total = 40.0_wp                        ! mm
    real(wp), parameter :: L_diag   = sqrt(L_bay**2 + H_total**2)   ! mm
    real(wp), parameter :: A_cross = 65.0_wp                        ! mm^2
    real(wp), parameter :: E_mod   = 200000.0_wp                    ! MPa (N/mm^2)

    ! Boundary Conditions and Loads
    real(wp), parameter :: P_load  = -90000.0_wp                    ! N (Downward)
    real(wp), parameter :: load_increment = 1000.0_wp               ! N (Downward)
    integer,  parameter :: n_steps = int(abs(P_load) / load_increment)
    real(wp), parameter :: load_location(2) = [L_total, H_total]    ! Apply load at top-right node
    real(wp), parameter :: fixed_X = 0.0_wp                         ! mm

    ! Mesh Arrays
    real(wp) :: X(n_node), Y(n_node)
    integer  :: conn(2, n_elem) 
    logical  :: is_fixed(ndof)
    real(wp) :: F_ref(ndof)
    real(wp) :: u_final(ndof)

    integer :: i, n1, n2, e
    character(len=100) :: out_folder = "Results/"
    
    ! System Clock Variables
    integer :: tick_start, tick_end, tick_rate
    real(wp) :: elapsed_time

    print *, "Starting 2D Co-rotational Truss Solver..."

    ! DEFINE NODAL COORDINATES (X, Y)
    ! Node numbering strategy: 
    ! Odd nodes (1, 3, 5...) on bottom chord
    ! Even nodes (2, 4, 6...) on top chord

    do i = 1, n_bays + 1
        X(2*i - 1) = real(i - 1, wp) * L_bay  ! Bottom chord
        Y(2*i - 1) = 0.0_wp
        X(2*i)     = real(i - 1, wp) * L_bay  ! Top chord
        Y(2*i)     = H_total
    end do
    
    ! DEFINE ELEMENT CONNECTIVITY
    e = 0
    do i = 1, n_bays
        e = e + 1
        conn(:, e) = [2*i - 1, 2*i + 1] ! Bottom chord element
    end do

    do i = 1, n_bays
        e = e + 1
        conn(:, e) = [2*i, 2*i + 2] ! Top chord element
    end do

    do i = 1, n_bays + 1
        e = e + 1
        conn(:, e) = [2*i - 1, 2*i] ! vertical element
    end do

    do i = 1, n_bays
        e = e + 1
        conn(:, e) = [2*i - 1, 2*i + 2] ! diagonal element
    end do


    ! BOUNDARY CONDITIONS & LOADS
    is_fixed = .false.
    F_ref = 0.0_wp

    ! We search the mesh. If a node is at X=0, we fix both its U and V DOFs.
    do i = 1, n_node
        if (abs(X(i) - fixed_X) < 1.0e-6_wp) then
            is_fixed(2*i - 1) = .true.  ! Fix X-DOF
            is_fixed(2*i)     = .true.  ! Fix Y-DOF
        end if
    end do


    ! INITIALIZE MATERIAL & SOLVER
    ! Units: MMGS, so E = 200,000 MPa = 200,000 N/mm^2
    my_mat = LinearElasticMaterial(E = 200000.0_wp)
    
    ! We pass the mesh data into our new solver class
    call model%init(n_node, n_elem, X, Y, conn, A=65.0_wp, mat_in=my_mat)


    ! EXECUTE INCREMENTAL NEWTON-RAPHSON

    call system_clock(tick_start, tick_rate)

    ! Load: 90 kN = 90,000 N. Applied in 1 kN increments -> 90 steps.
    ! Note: You will need to specify WHICH node/DOF gets the load inside solve_incremental.
    call model%solve_incremental(F_total = -90000.0_wp, & 
                                 n_steps = 90, & 
                                 max_iter = 50, &
                                 u_final = u_final, &
                                 out_folder = trim(out_folder))

    call system_clock(tick_end)
    elapsed_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)
    print '(A, F0.4, A)', ">>> Solver Time: ", elapsed_time, " seconds <<<"

end program main