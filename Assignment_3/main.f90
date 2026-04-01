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

    ! DYNAMIC VARIABLES
    integer :: n_bays, n_elem, n_node, ndof
    real(wp) :: L_total, L_bay, H_total, L_diag, A_cross, E_mod
    integer :: num_supports, num_loads, n_steps
    
    ! Dynamically allocated arrays
    real(wp), allocatable :: support_data(:,:)
    real(wp), allocatable :: load_data(:,:)
    
    ! File IO variables
    integer :: in_unit, dot_idx, i, slash_idx
    character(len=256) :: input_file, base_name, prefix

    ! Variables for the search routines
    logical :: node_found
    integer :: j, e

    ! Mesh and Global Arrays
    real(wp), allocatable :: X(:), Y(:)
    integer, allocatable  :: conn(:,:)
    
    ! Solver Arrays
    real(wp), allocatable :: F_ref(:)
    logical, allocatable  :: is_fixed(:)
    real(wp), allocatable :: u_final(:)

    ! Output Variables
    integer :: csv_nodes, csv_elems

    ! System Clock Variables
    integer :: tick_start, tick_end, tick_rate
    real(wp) :: elapsed_time

    ! READ INPUT FILE
    ! Read filename from the terminal (e.g., ./fea_solver job1.inp)
    if (command_argument_count() >= 1) then
        call get_command_argument(1, input_file)
    else
        print *, "No input file specified. Defaulting to 'job1.inp'"
        input_file = "job1.inp"
    end if

    ! Find the last slash to strip any folder paths (e.g., "inputs/prob1.in" to "prob1.in")
    slash_idx = index(trim(input_file), '/', back=.true.)
    if (slash_idx > 0) then
        base_name = input_file(slash_idx+1 : len_trim(input_file))
    else
        base_name = trim(input_file)
    end if

    ! Extract base name for outputs (turns "job1.inp" into "job1")
    dot_idx = index(trim(base_name), '.', back=.true.)
    if (dot_idx > 0) then
        base_name = base_name(1 : dot_idx-1)
    end if
    
    prefix = "outputs/" // trim(base_name) // "_"

    ! Open and parse the text file
    print *, "Reading input file: ", trim(input_file)
    open(newunit=in_unit, file=trim(input_file), status='old')
    
    read(in_unit, *) n_bays
    read(in_unit, *) L_total
    read(in_unit, *) H_total
    read(in_unit, *) A_cross
    read(in_unit, *) E_mod
    read(in_unit, *) n_steps
    
    read(in_unit, *) num_supports
    allocate(support_data(4, num_supports))
    do i = 1, num_supports
        read(in_unit, *) support_data(1:4, i)
    end do
    
    read(in_unit, *) num_loads
    allocate(load_data(4, num_loads))
    do i = 1, num_loads
        read(in_unit, *) load_data(1:4, i)
    end do
    
    close(in_unit)

    ! CALCULATE DERIVED PARAMETERS
    n_elem = (4 * n_bays) + 1
    n_node = 2 * n_bays + 2
    ndof   = n_node * 2
    L_bay  = L_total / real(n_bays, wp)
    L_diag = sqrt(L_bay**2 + H_total**2)

    allocate(X(n_node))
    allocate(Y(n_node))
    allocate(conn(2, n_elem))
    allocate(F_ref(ndof))
    allocate(is_fixed(ndof))
    allocate(u_final(ndof))

    ! Initialize them to zero
    X = 0.0_wp; Y = 0.0_wp; conn = 0
    F_ref = 0.0_wp; is_fixed = .false.; u_final = 0.0_wp
    
    print *, "Initialization successful. N_DOFs: ", ndof

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
    ! We search the mesh to find fixed nodes.
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
    open(newunit=csv_nodes, file=trim(prefix)//'undeformed_nodes.csv', status='replace')
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
    open(newunit=csv_elems, file=trim(prefix)//'undeformed_elements.csv', status='replace')
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
                                 prefix = trim(prefix))

    
    call system_clock(tick_end)
    elapsed_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    ! Export Deformed Mesh for Visualization
    print *, "Exporting final deformed shape..."
    open(newunit=csv_nodes, file=trim(prefix)//'deformed_nodes.csv', status='replace')
    
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
    print *, "Simulation Complete! Data exported to: ", trim(prefix)

end program main