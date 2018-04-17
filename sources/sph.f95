! ===========================================================================================================
! ===========================================================================================================

MODULE sph

    use math

    implicit none

contains

    ! -------------------------------------------------------------------------------------------------------
    ! créer un quadrillage à partir d'une subdivision de x et de y
    ! -------------------------------------------------------------------------------------------------------
    ! x1 : subdivision des abcisse
    ! x2 : subdivision des ordonnées
    subroutine meshgrid(x1, x2, x)
        ! paramètres
        real(rp), dimension(:), intent(in) :: x1, x2
        real(rp), dimension(:, :), allocatable, intent(out) :: x

        ! variables locales
        integer :: i, j, k

        allocate(x(size(x1) * size(x2), 2))
        k = 1

        do i = 1, size(x1)
            do j = 1, size(x2)
                x(k, :) = (/ x1(i), x2(j) /)
                k = k + 1
            end do
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! fonction pour noyau SPH
    ! -------------------------------------------------------------------------------------------------------
    function theta(q)
        ! paramètres
        real(rp), intent(in) :: q

        ! return
        real(rp) :: theta

        ! variables locales
        real(rp) :: C

        C = 7.0_rp

        if ((0.0_rp <= q) .and. (q <= 1.0_rp)) then
            theta = C * (1.0_rp - q)**4 * (1.0_rp + 4.0_rp*q)
        else
            theta = 0.0_rp
        end if
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Noyau SPH
    ! -------------------------------------------------------------------------------------------------------
    ! x : tableau des i particules (lignes) et de 2 coordonnées (colonnes)
    ! R : rayon SPH
    function W_SPH(x, R) result(W)
        ! paramètres
        real(rp), dimension(2), intent(in) :: x
        real(rp), intent(in) :: R

        ! return
        real(rp) :: W

        ! variables locales
        real(rp) :: norm

        call norme2(x, norm)
        if (norm <= R) then
            W = theta(norm / R) / (pi * R**2)
        else
            W = 0.0_rp
        end if
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Gradient en x de W(x, y, R) = W(x - y, R)
    ! -------------------------------------------------------------------------------------------------------
    subroutine dx_W_SPH(z, R, grad)
        ! paramètres
        real(rp), dimension(2), intent(in) :: z
        real(rp), intent(in) :: R
        real(rp), dimension(2), intent(out) :: grad

        ! variables locales
        !real(rp) :: u, v
        !real(rp), dimension(2) :: up, vp
        real(rp) :: norm
        real(rp) :: q

        call norme2(z, norm)

        if (norm <= R) then
            q = norm / R
            grad = (7.0_rp / (pi * R**2)) * (-20.0_rp / R) * (z / norm) * q * (1-q)**3
        else
            grad = (/ 0.0_rp, 0.0_rp /)
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Approximation régularisée    Permet vérification Somme(wj * W_SPH_ij) = 1 pour f = 1 partout
    ! -------------------------------------------------------------------------------------------------------
    subroutine AR(i, x, w, R, f, image)
        ! paramètres
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), intent(in) :: R
        real(rp), dimension(:), intent(in) :: f
        real(rp), intent(out) :: image

        ! variables locales
        integer :: j

        image = 0.0_rp

        do j = 1, size(w)
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                image = image + w(j) * f(j) * W_SPH(x(i, :) - x(j, :), R)
            end if
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Gradient régularisé      Permet vérification Somme(wj * \nabla_x W_SPH_ij) = 0 pour f = 1 partout
    ! -------------------------------------------------------------------------------------------------------
    subroutine GR(i, x, w, R, f, image)
        ! paramètres
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), intent(in) :: R
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(2), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(2) :: grad

        image = 0.0_rp

        do j = 1, i - 1
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * f(j) * grad
            end if
        end do

        do j = i + 1, size(w)
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * f(j) * grad
            end if
        end do
    end subroutine

    subroutine GR_m(i, x, w, R, f, image)
        ! paramètres
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), intent(in) :: R
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(2), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(2) :: grad

        image = 0.0_rp

        do j = 1, i - 1
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * (f(j) - f(i)) * grad
            end if
        end do

        do j = i + 1, size(w)
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * (f(j) - f(i)) * grad
            end if
        end do
    end subroutine

    subroutine GR_p(i, x, w, R, f, image)
        ! paramètres
        integer, intent(in) :: i
        real(rp), dimension(:, :), intent(in) :: x
        real(rp), dimension(:), intent(in) :: w
        real(rp), intent(in) :: R
        real(rp), dimension(:), intent(in) :: f
        real(rp), dimension(2), intent(out) :: image

        ! variables locales
        integer :: j
        real(rp), dimension(2) :: grad

        image = 0.0_rp

        do j = 1, i - 1
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * (f(j) + f(i)) * grad
            end if
        end do

        do j = i + 1, size(w)
            if (fnorme2(x(i, :) - x(j, :)) <= R) then
                call dx_W_SPH(x(i, :) - x(j, :), R, grad)
                image = image + w(j) * (f(j) + f(i)) * grad
            end if
        end do
    end subroutine

END MODULE sph
