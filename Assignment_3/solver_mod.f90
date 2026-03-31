! ==========================================
! THE SOLVER
! ==========================================
module solver_mod
    use kinds_mod
    use material_mod
    implicit none

    ! The FEA Class
    type :: CorotationalTrussFEA
        integer :: n_node, n_elem, ndof
        real(wp) :: A
        class(Material), allocatable :: mat
        ! Arrays allocated during init
        real(wp), allocatable :: X0(:), Y0(:)    ! Original coordinates
        integer, allocatable  :: conn(:,:)       ! Connectivity
    contains
        procedure :: init => init_solver
        procedure :: assemble_system
        procedure :: solve_incremental
    end type CorotationalTrussFEA

contains

    subroutine init_solver(self, n_node, n_elem, X, Y, conn, A, mat_in)
        class(CorotationalTrussFEA), intent(inout) :: self
        integer, intent(in)  :: n_node, n_elem
        real(wp), intent(in) :: X(n_node), Y(n_node)
        integer, intent(in)  :: conn(2, n_elem)
        real(wp), intent(in) :: A
        class(Material), intent(in) :: mat_in

        self%n_node = n_node
        self%n_elem = n_elem
        self%ndof   = 2 * n_node
        self%A      = A

        if (allocated(self%X0))   deallocate(self%X0)
        if (allocated(self%Y0))   deallocate(self%Y0)
        if (allocated(self%conn)) deallocate(self%conn)
        if (allocated(self%mat))  deallocate(self%mat)

        allocate(self%X0, source=X)
        allocate(self%Y0, source=Y)
        allocate(self%conn, source=conn)

        allocate(self%mat, source=mat_in)

    end subroutine init_solver

    subroutine assemble_system(self, u, K_global, F_int)
        class(CorotationalTrussFEA), intent(in) :: self
        real(wp), intent(in)  :: u(self%ndof)          ! Current global displacements
        real(wp), intent(out) :: K_global(self%ndof, self%ndof) ! Dense global stiffness
        real(wp), intent(out) :: F_int(self%ndof)      ! Global internal force vector

        integer :: e, n1, n2, i, j
        integer :: dof(4)
        real(wp) :: x1, y1, x2, y2          ! Current deformed coordinates
        real(wp) :: x1_0, y1_0, x2_0, y2_0  ! Original coordinates
        real(wp) :: L0, L, c, s, delta, strain, N_force
        real(wp) :: stress(1), Et(1)
        
        real(wp) :: Ke(4,4), Kg(4,4), K_elem(4,4)
        real(wp) :: fe(4)

        ! Initialize global arrays to zero
        K_global = 0.0_wp
        F_int    = 0.0_wp

        ! ELEMENT LOOP
        do e = 1, self%n_elem
            ! Get Node IDs for this element
            n1 = self%conn(1, e)
            n2 = self%conn(2, e)

            ! Map local DOFs to global DOFs
            dof(1) = 2*n1 - 1
            dof(2) = 2*n1
            dof(3) = 2*n2 - 1
            dof(4) = 2*n2

            ! Original coordinates
            x1_0 = self%X0(n1); y1_0 = self%Y0(n1)
            x2_0 = self%X0(n2); y2_0 = self%Y0(n2)

            ! Original length (L0)
            L0 = sqrt((x2_0 - x1_0)**2 + (y2_0 - y1_0)**2)

            ! Current deformed coordinates (X0 + u)
            x1 = x1_0 + u(dof(1))
            y1 = y1_0 + u(dof(2))
            x2 = x2_0 + u(dof(3))
            y2 = y2_0 + u(dof(4))

            ! Current length (L) and angles (c, s)
            L = sqrt((x2 - x1)**2 + (y2 - y1)**2)
            c = (x2 - x1) / L
            s = (y2 - y1) / L

            ! Strain, stress, and internal force
            delta = L - L0
            strain = delta / L0

            call self%mat%get_stress_and_tangent([strain], stress, Et, 1)

            N_force = stress(1) * self%A

            ! Material Stiffness Matrix (Ke) and Geometric Stiffness Matrix (Kg)
            Ke(1,:) = (Et(1) * self%A) / L0 * [ c*c,  c*s, -c*c, -c*s]
            Ke(2,:) = (Et(1) * self%A) / L0 * [ c*s,  s*s, -c*s, -s*s]
            Ke(3,:) = (Et(1) * self%A) / L0 * [-c*c, -c*s,  c*c,  c*s]
            Ke(4,:) = (Et(1) * self%A) / L0 * [-c*s, -s*s,  c*s,  s*s]

            Kg(1,:) = (N_force / L) * [ s*s, -c*s, -s*s,  c*s]
            Kg(2,:) = (N_force / L) * [-c*s,  c*c,  c*s, -c*c]
            Kg(3,:) = (N_force / L) * [-s*s,  c*s,  s*s, -c*s]
            Kg(4,:) = (N_force / L) * [ c*s, -c*c, -c*s,  c*c]

            ! Element Internal Force Vector (fe)
            fe = N_force * [-c, -s, c, s]

            ! Global Matrices
            K_elem = Ke + Kg
            
            ! Assemble into global K and F_int
            do i = 1, 4
                F_int(dof(i)) = F_int(dof(i)) + fe(i)
            end do

            do j = 1, 4
                do i = 1, 4
                    K_global(dof(i), dof(j)) = K_global(dof(i), dof(j)) + K_elem(i, j)
                end do
            end do

        end do
    end subroutine assemble_system

    subroutine solve_incremental(self, F_ref, is_fixed, n_steps, max_iter, u_final, out_folder)
        class(CorotationalTrussFEA), intent(in) :: self
        real(wp), intent(in)  :: F_ref(self%ndof)       ! Reference load vector (100% load)
        logical, intent(in)   :: is_fixed(self%ndof)    ! Boundary condition flags
        integer, intent(in)   :: n_steps, max_iter
        real(wp), intent(out) :: u_final(self%ndof)
        character(len=*), intent(in) :: out_folder

        ! Linear Algebra Arrays
        real(wp) :: K_global(self%ndof, self%ndof)
        real(wp) :: F_int(self%ndof), F_ext(self%ndof), R(self%ndof)
        real(wp) :: du(self%ndof)
        
        ! LAPACK dgesv Variables
        integer :: ipiv(self%ndof)
        integer :: info

        ! Control & Output Variables
        integer  :: step, iter, i, csv_id
        real(wp) :: load_factor, norm_du
        
        ! INITIALIZATION
        u_final = 0.0_wp
        
        open(newunit=csv_id, file=trim(out_folder)//'history.csv', status='replace')
        ! Tracks the Y-displacement of the loaded node (last DOF) vs total applied load
        write(csv_id, '(A)') "Load_Factor,Tip_Disp_Y,Applied_Force_Y"
        write(csv_id, '(3(ES15.6, A))') 0.0_wp, ",", 0.0_wp, ",", 0.0_wp, ""

        ! Outer Loop
        do step = 1, n_steps
            load_factor = real(step, wp) / real(n_steps, wp)
            F_ext = F_ref * load_factor

            ! NEWTON-RAPHSON ITERATION (Inner Loop)
            do iter = 1, max_iter
                
                ! Current tangent stiffness and internal force
                call self%assemble_system(u_final, K_global, F_int)
                
                ! Residual
                R = F_ext - F_int
                
                ! Dirichlet Boundary Conditions
                do i = 1, self%ndof
                    if (is_fixed(i)) then
                        K_global(i, :) = 0.0_wp
                        K_global(:, i) = 0.0_wp
                        K_global(i, i) = 1.0_wp
                        R(i) = 0.0_wp
                    end if
                end do
                
                ! [K_global] {du} = {R}
                call dgesv(self%ndof, 1, K_global, self%ndof, ipiv, R, self%ndof, info)
                
                ! Error Handling for LAPACK
                if (info /= 0) then
                    print *, "FATAL ERROR: LAPACK dgesv failed with info = ", info
                    if (info > 0) then
                        print *, "The stiffness matrix is singular at step ", step
                        print *, "This means your structure is a mechanism. Check boundary conditions!"
                    end if
                    stop
                end if
                
                ! Update Displacements
                du = R
                u_final = u_final + du
                
                ! Convergence Check
                norm_du = sqrt(sum(du**2))
                if (norm_du < 1e-8_wp) then
                    exit  ! Converged
                end if
                
            end do
            
            ! Print step convergence cleanly
            print '(A, I0, A, I0, A, I0, A, I0, A)', &
                "Step ", step, "/", n_steps, &
                " converged after ", iter, "/", max_iter, " iterations."

            ! Log the Y-displacement of the top-right tip
            write(csv_id, '(ES15.6,A,ES15.6,A,ES15.6)') &
                load_factor, ",", u_final(self%ndof), ",", F_ext(self%ndof)

        end do
        
        close(csv_id)
    end subroutine solve_incremental

end module solver_mod