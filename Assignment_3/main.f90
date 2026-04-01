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
    integer, parameter :: n_elem = (4 * n_bays) + 1  ! 5 elements per bay: bottom, top, left, right, diagonal
    integer, parameter :: n_node = 2 * n_bays + 2
    integer, parameter :: ndof   = n_node * 2
    
    ! Material and Geometric Properties
    real(wp), parameter :: L_total = 240.0_wp                       ! mm
    real(wp), parameter :: L_bay   = L_total / real(n_bays, wp)     ! mm
    real(wp), parameter :: H_total = 40.0_wp                        ! mm
    real(wp), parameter :: L_diag   = sqrt(L_bay**2 + H_total**2)   ! mm
    real(wp), parameter :: A_cross = 65.0_wp                        ! mm^2
    real(wp), parameter :: E_mod   = 200000.0_wp                    ! MPa (N/mm^2)

    ! Boundary Conditions and Loads
    ! Dirichlet
    integer, parameter :: num_supports = 2  ! Increase to add more supports (e.g. roller at mid-span)
    
    ! Matrix format: [X_coord, Y_coord, Fix_X, Fix_Y] per column
    ! Fix flags: 1.0 = Fixed, 0.0 = Free
    real(wp), parameter :: support_data(4, num_supports) = reshape([ &
        0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, &
        0.0_wp, 40.0_wp, 1.0_wp, 1.0_wp  &
        ! If you wanted a roller at mid-span, you could add:
        ! 120.0_wp, 0.0_wp, 0.0_wp, 1.0_wp 
    ], [4, num_supports])

    ! Neumann
    integer, parameter :: num_loads = 1  ! Increase to add more loads (e.g. point load at mid-span)
    
    ! Matrix format: [X_coord, Y_coord, Force_X, Force_Y] per column
    real(wp), parameter :: load_data(4, num_loads) = reshape([ &
        240.0_wp, 40.0_wp, 0.0_wp, -90000.0_wp & ! Load 1: Tip, Downward
        ! 240.0_wp, 0.0_wp, 90000.0_wp, 0.0_wp & ! Load 2: Tip, Rightward
    ], [4, num_loads])
    integer, parameter :: n_steps = 90


    ! Variables for the search routines
    logical :: node_found
    integer :: j, i, e

    ! Mesh Arrays
    real(wp) :: X(n_node), Y(n_node)
    integer  :: conn(2, n_elem) 
    logical  :: is_fixed(ndof)
    real(wp) :: F_ref(ndof)
    real(wp) :: u_final(ndof)

    ! Output Variables
    character(len=100) :: out_folder = "outputs/"
    integer :: csv_nodes, csv_elems

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

    ! We search the mesh to find fixed nodes.
    is_fixed = .false.
    do j = 1, num_supports
        node_found = .false.
        do i = 1, n_node
            if (abs(X(i) - support_data(1,j)) < 1.0e-6_wp .and. &
                abs(Y(i) - support_data(2,j)) < 1.0e-6_wp) then
                
                ! Apply constraints based on flags
                if (support_data(3,j) == 1.0_wp) is_fixed(2*i - 1) = .true.
                if (support_data(4,j) == 1.0_wp) is_fixed(2*i)     = .true.
                
                node_found = .true.
                exit ! Node found, move to next support
            end if
        end do
        
        if (.not. node_found) then
            print *, "FATAL ERROR: Support specified at invalid coordinate:"
            print *, "X = ", support_data(1,j), " Y = ", support_data(2,j)
            print *, "Truss supports must be applied exactly at nodes."
            stop
        end if
    end do

    ! We search the mesh to find loaded nodes.
    F_ref = 0.0_wp
    
    do j = 1, num_loads
        node_found = .false.
        do i = 1, n_node
            if (abs(X(i) - load_data(1,j)) < 1.0e-6_wp .and. &
                abs(Y(i) - load_data(2,j)) < 1.0e-6_wp) then
                
                ! Add the load to the reference vector
                F_ref(2*i - 1) = F_ref(2*i - 1) + load_data(3,j)
                F_ref(2*i)     = F_ref(2*i)     + load_data(4,j)
                
                node_found = .true.
                exit
            end if
        end do
        
        if (.not. node_found) then
            print *, "FATAL ERROR: Load specified at invalid coordinate:"
            print *, "X = ", load_data(1,j), " Y = ", load_data(2,j)
            print *, "Truss loads must be applied exactly at nodes."
            stop
        end if
    end do

    ! EXPORT UNDEFORMED MESH FOR PYTHON VERIFICATION
    
    print *, "Exporting undeformed mesh data to CSV..."
    
    ! Export Nodal Data (Coordinates, Boundary Conditions, Loads)
    open(newunit=csv_nodes, file=trim(out_folder)//'undeformed_nodes.csv', status='replace')
    write(csv_nodes, '(A)') "Node_ID,X,Y,Fix_X,Fix_Y,Force_X,Force_Y"
    do i = 1, n_node
        write(csv_nodes, '(I0, 6(A, ES15.6))') &
            i, ",", X(i), ",", Y(i), &
            ",", merge(1.0_wp, 0.0_wp, is_fixed(2*i - 1)), &
            ",", merge(1.0_wp, 0.0_wp, is_fixed(2*i)), &
            ",", F_ref(2*i - 1), ",", F_ref(2*i)
    end do
    close(csv_nodes)

    ! Export Element Connectivity
    open(newunit=csv_elems, file=trim(out_folder)//'undeformed_elements.csv', status='replace')
    write(csv_elems, '(A)') "Elem_ID,Node_i,Node_j"
    do i = 1, n_elem
        write(csv_elems, '(I0, A, I0, A, I0)') i, ",", conn(1,i), ",", conn(2,i)
    end do
    close(csv_elems)
    
    print *, "Pre-processing complete."


    ! INITIALIZE MATERIAL & SOLVER
    
    print *, "Initializing Material and FEA Model..."
    
    ! Instantiate the material using the Input Deck parameters
    my_mat = LinearElasticMaterial(E = E_mod)
    
    ! Pass the mesh arrays and properties to the solver class
    call model%init(n_node, n_elem, X, Y, conn, A=A_cross, mat_in=my_mat)

    
    ! EXECUTE INCREMENTAL NEWTON-RAPHSON SOLVER
    
    print *, "Starting Incremental Newton-Raphson Solver..."
    call system_clock(tick_start, tick_rate)

    call model%solve_incremental(F_ref = F_ref, & 
                                 is_fixed = is_fixed, &
                                 n_steps = n_steps, & 
                                 max_iter = 50, &        
                                 u_final = u_final, &
                                 out_folder = trim(out_folder))

    
    call system_clock(tick_end)
    elapsed_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    ! Export Deformed Mesh for Visualization
    print *, "Exporting final deformed shape..."
    open(newunit=csv_nodes, file=trim(out_folder)//'deformed_nodes.csv', status='replace')
    
    ! We export the Original X, Y and the Displacements U, V
    write(csv_nodes, '(A)') "Node_ID,X_orig,Y_orig,U,V,X_def,Y_def"
    
    do i = 1, n_node
        ! Extract U and V for this node from the u_final vector
        write(csv_nodes, '(I0, 6(A, ES15.6))') &
            i, &
            ",", X(i), &
            ",", Y(i), &
            ",", u_final(2*i - 1), &  ! U displacement
            ",", u_final(2*i), &      ! V displacement
            ",", X(i) + u_final(2*i - 1), & ! Final Deformed X
            ",", Y(i) + u_final(2*i)        ! Final Deformed Y
    end do
    close(csv_nodes)
    
    print *, "=========================================="
    print '(A, F0.4, A)', " >>> Fortran Core Solver Time: ", elapsed_time, " seconds <<<"
    print *, "=========================================="
    print *, "Simulation Complete! Data exported to: ", trim(out_folder)

end program main