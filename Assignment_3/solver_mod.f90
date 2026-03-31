! ==========================================
! THE SOLVER
! ==========================================
module solver_mod
    use kinds_mod
    use material_mod
    implicit none

    ! The FEA Class
    type :: CorotationalTrussFEA
        real(wp) :: L, A
        integer :: n_elem, n_node, n_gauss
        class(Material), allocatable :: mat  ! Polymorphic material object
        real(wp), allocatable :: gauss_pts(:), gauss_wts(:)
    contains
        procedure :: init => init_solver
        procedure :: assemble_system
        procedure :: solve_incremental
    end type CorotationalTrussFEA

contains

    subroutine init_solver(self, L, A, mat_in, n_elem, n_gauss)
        class(CorotationalTrussFEA), intent(inout) :: self
        real(wp), intent(in) :: L, A
        class(Material), intent(in) :: mat_in
        integer, intent(in) :: n_elem, n_gauss

        self%L = L
        self%A = A
        self%n_elem = n_elem
        self%n_node = n_elem + 1
        self%n_gauss = n_gauss

        ! Allocate the polymorphic class using the incoming object
        allocate(self%mat, source=mat_in) 

        select case (n_gauss)
            case (1)
                self%gauss_pts = [0.0_wp]
                self%gauss_wts = [2.0_wp]
            case (2)
                self%gauss_pts = [-1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp)]
                self%gauss_wts = [1.0_wp, 1.0_wp]
            case (3)
                self%gauss_pts = [-sqrt(0.6_wp), 0.0_wp, sqrt(0.6_wp)]
                self%gauss_wts = [5.0_wp/9.0_wp, 8.0_wp/9.0_wp, 5.0_wp/9.0_wp]
            case default
                print *, "FATAL ERROR: Invalid number of Gauss points requested."
                stop
        end select

    end subroutine init_solver

subroutine assemble_system(self, u, ld, d, ud, f_int)
        class(CorotationalTrussFEA), intent(in) :: self
        real(wp), intent(in) :: u(self%n_node)
        real(wp), intent(out) :: ld(self%n_node - 1) ! Lower diagonal
        real(wp), intent(out) :: d(self%n_node)      ! Main diagonal
        real(wp), intent(out) :: ud(self%n_node - 1) ! Upper diagonal
        real(wp), intent(out) :: f_int(self%n_node)

        integer :: i, g, n1, n2
        real(wp) :: le, detJ, invJ, B1, B2
        real(wp) :: strain(1), stress(1), Et(1), factor, f_factor
        real(wp) :: ke(2,2), fe(2)

        ! Initialize sparse arrays to zero
        ld = 0.0_wp
        d = 0.0_wp
        ud = 0.0_wp
        f_int = 0.0_wp

        le = self%L / real(self%n_elem, wp)
        detJ = le / 2.0_wp
        invJ = 1.0_wp / detJ
        B1 = -0.5_wp * invJ
        B2 =  0.5_wp * invJ

        ! Element Loop
        do i = 1, self%n_elem
            n1 = i
            n2 = i + 1
            ke = 0.0_wp
            fe = 0.0_wp

            ! Gauss Loop
            do g = 1, self%n_gauss
                strain(1) = B1 * u(n1) + B2 * u(n2)
                
                ! Call the polymorphic material object
                call self%mat%get_stress_and_tangent(strain, stress, Et, 1)

                factor = Et(1) * self%A * detJ * self%gauss_wts(g)
                ke(1,1) = ke(1,1) + B1 * B1 * factor
                ke(1,2) = ke(1,2) + B1 * B2 * factor
                ke(2,1) = ke(2,1) + B2 * B1 * factor
                ke(2,2) = ke(2,2) + B2 * B2 * factor

                f_factor = stress(1) * self%A * detJ * self%gauss_wts(g)
                fe(1) = fe(1) + B1 * f_factor
                fe(2) = fe(2) + B2 * f_factor
            end do

            ! Direct Sparse Assembly into Tridiagonal Bands
            d(n1)  = d(n1)  + ke(1,1)
            d(n2)  = d(n2)  + ke(2,2)
            
            ud(n1) = ud(n1) + ke(1,2)
            ld(n1) = ld(n1) + ke(2,1)
            
            f_int(n1) = f_int(n1) + fe(1)
            f_int(n2) = f_int(n2) + fe(2)
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