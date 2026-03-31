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

    subroutine solve_incremental(self, F_total, n_steps, max_iter, u_final, out_folder)
        class(CorotationalTrussFEA), intent(in) :: self
        real(wp), intent(in) :: F_total
        integer, intent(in) :: n_steps
        integer, intent(in) :: max_iter
        real(wp), intent(out) :: u_final(self%n_node)
        character(len=*), intent(in) :: out_folder

        ! LAPACK variables
        integer :: info
        
        ! Tridiagonal arrays (Lower, Main, Upper)
        real(wp) :: ld(self%n_node - 1)
        real(wp) :: d(self%n_node)
        real(wp) :: ud(self%n_node - 1)

        integer :: step, iter, csv_id
        real(wp) :: F_ext(self%n_node), f_int(self%n_node)
        real(wp) :: R(self%n_node)
        real(wp) :: du(self%n_node), force_inc, norm_du
        real(wp) :: strain_last(1), stress_last(1), Et_last(1), le

        u_final = 0.0_wp
        force_inc = F_total / real(n_steps, wp)

        ! Open CSV for history
        open(newunit=csv_id, file=trim(out_folder)//'history.csv', status='replace')
        write(csv_id, '(A)') "Force,Tip_Disp,Last_Strain,Last_Stress"
        write(csv_id, '(ES15.6,A,ES15.6,A,ES15.6,A,ES15.6)') 0.0_wp, ",", 0.0_wp, ",", 0.0_wp, ",", 0.0_wp

        le = self%L / real(self%n_elem, wp)

        do step = 1, n_steps
            F_ext = 0.0_wp
            F_ext(self%n_node) = force_inc * step

            do iter = 1, max_iter

                call self%assemble_system(u_final, ld, d, ud, f_int)
                R = F_ext - f_int
                
                d(1) = 1.0_wp
                ud(1) = 0.0_wp
                ld(1) = 0.0_wp   
                R(1) = 0.0_wp

                call dgtsv(self%n_node, 1, ld, d, ud, R, self%n_node, info)
                
                if (info /= 0) then
                    print *, "LAPACK Error: Matrix is singular! INFO = ", info
                    stop
                end if

                ! R is overwritten with the displacement increment
                du = R  
                u_final = u_final + du

                ! Convergence Check (L2 Norm of the full du vector)
                norm_du = sqrt(sum(du**2))
                if (norm_du < 1.0e-8_wp) exit
            end do
            print '(A, I0, A, I0, A, I0, A, I0, A)', &
                    "Step ", step, "/", n_steps, &
                    " converged after ", iter, "/", max_iter, " iterations."
            strain_last(1) = (u_final(self%n_node) - u_final(self%n_node-1)) / le
            call self%mat%get_stress_and_tangent(strain_last, stress_last, Et_last, 1)
            
            ! Write history data
            write(csv_id, '(ES15.6,A,ES15.6,A,ES15.6,A,ES15.6)') & 
                F_ext(self%n_node), ",", u_final(self%n_node), ",", strain_last(1), ",", stress_last(1)
        end do
        close(csv_id)
    end subroutine solve_incremental

end module solver_mod