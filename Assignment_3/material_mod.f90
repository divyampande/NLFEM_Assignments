! ==========================================
! THE PHYSICS
! ==========================================
module material_mod
    use kinds_mod
    implicit none
    private
    public :: Material, LinearElasticMaterial, ArctanMaterial, BilinearYieldMaterial

    type, abstract :: Material
    contains
        procedure(get_stress_tangent_intf), deferred, pass(self) :: get_stress_and_tangent
    end type Material

    abstract interface
        subroutine get_stress_tangent_intf(self, strain, stress, Et, n)
            import :: Material, wp
            class(Material), intent(in) :: self
            integer, intent(in) :: n
            real(wp), intent(in) :: strain(n)
            real(wp), intent(out) :: stress(n), Et(n)
        end subroutine
    end interface
    
    ! Linear Elastic Material
    type, extends(Material) :: LinearElasticMaterial
        real(wp) :: E
    contains
        procedure, pass(self) :: get_stress_and_tangent => linear_eval
    end type LinearElasticMaterial

    ! Arctan Material
    type, extends(Material) :: ArctanMaterial
        real(wp) :: E, m
    contains
        procedure, pass(self) :: get_stress_and_tangent => arctan_eval
    end type ArctanMaterial

    ! Plastic Material
    type, extends(Material) :: BilinearYieldMaterial
        real(wp) :: E, yield_stress, yield_strain, Et_plastic
    contains
        procedure, pass(self) :: get_stress_and_tangent => bilinear_eval
    end type BilinearYieldMaterial

    ! Constructor interface
    interface BilinearYieldMaterial
        module procedure new_bilinear
    end interface

contains
    ! Linear Elastic Material Evaluation (Trivial, but we need it for the base case)
    subroutine linear_eval(self, strain, stress, Et, n)
        class(LinearElasticMaterial), intent(in) :: self
        integer, intent(in) :: n
        real(wp), intent(in) :: strain(n)
        real(wp), intent(out) :: stress(n), Et(n)

        integer :: i
        do i = 1, n
            stress(i) = self%E * strain(i)
            Et(i) = self%E
        end do
    end subroutine linear_eval

    ! Arctan Material Evaluation (Material nonlinearity without plasticity)
    subroutine arctan_eval(self, strain, stress, Et, n)
        class(ArctanMaterial), intent(in) :: self
        integer, intent(in) :: n
        real(wp), intent(in) :: strain(n)
        real(wp), intent(out) :: stress(n), Et(n)

        stress = self%E * atan(strain * self%m)
        Et = (self%E * self%m) / (1.0_wp + (strain * self%m)**2)

    end subroutine arctan_eval

    ! Constructor Function for Bilinear Material
    function new_bilinear(E, yield_stress, Et_plastic) result(mat)
        real(wp), intent(in) :: E, yield_stress, Et_plastic
        type(BilinearYieldMaterial) :: mat
        
        mat%E = E
        mat%yield_stress = yield_stress
        mat%yield_strain = yield_stress / E
        mat%Et_plastic = Et_plastic
    end function new_bilinear

    ! Bilinear Material Evaluation (Elastic-perfectly plastic with linear hardening)
    subroutine bilinear_eval(self, strain, stress, Et, n)
        class(BilinearYieldMaterial), intent(in) :: self
        integer, intent(in) :: n
        real(wp), intent(in) :: strain(n)
        real(wp), intent(out) :: stress(n), Et(n)
        
        integer :: i
        real(wp) :: excess, sign_val
        do i = 1, n
            if (abs(strain(i)) > self%yield_strain) then
                sign_val = sign(1.0_wp, strain(i))
                excess = abs(strain(i)) - self%yield_strain
                stress(i) = sign_val * (self%yield_stress + self%Et_plastic * excess)
                Et(i) = self%Et_plastic
            else
                stress(i) = self%E * strain(i)
                Et(i) = self%E
            end if
        end do
    end subroutine bilinear_eval
end module material_mod