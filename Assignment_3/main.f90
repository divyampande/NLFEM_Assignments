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
    
    ! Mesh Arrays
    real(wp) :: X(n_node), Y(n_node)
    integer  :: conn(2, n_elem) 
    real(wp) :: u_final(ndof)
    
    integer :: i
    character(len=100) :: out_folder = "Results/"
    
    ! System Clock Variables
    integer :: tick_start, tick_end, tick_rate
    real(wp) :: elapsed_time

    print *, "Starting 2D Co-rotational Truss Solver..."

    ! DEFINE NODAL COORDINATES (X, Y)
    ! Node numbering strategy: 
    ! Odd nodes (1, 3, 5...) on bottom chord
    ! Even nodes (2, 4, 6...) on top chord

    ! >> TODO: Populate X and Y arrays here based on L=240 (bottom) and L=240 (top)
    ! Example for first bay:
    ! X(1) = 0.0_wp;   Y(1) = 0.0_wp     ! Bottom Wall
    ! X(2) = 0.0_wp;   Y(2) = 40.0_wp    ! Top Wall
    ! X(3) = 40.0_wp;  Y(3) = 0.0_wp
    ! X(4) = 40.0_wp;  Y(4) = 40.0_wp
    ! ...
    

    ! DEFINE ELEMENT CONNECTIVITY
    ! >> TODO: Populate conn(1,:) and conn(2,:) mapping elements to nodes.
    ! Example:
    ! conn(:, 1) = [1, 3]  ! Bottom chord, Bay 1
    ! conn(:, 2) = [2, 4]  ! Top chord, Bay 1
    ! conn(:, 3) = [1, 2]  ! Wall vertical
    ! conn(:, 4) = [3, 4]  ! Vertical, Bay 1
    ! conn(:, 5) = [1, 4]  ! Diagonal, Bay 1
    ! ...

    ! 3. INITIALIZE MATERIAL & SOLVER
    ! Units: E = 200 GPa = 200,000 MPa (N/mm^2). Area = 65 mm^2
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