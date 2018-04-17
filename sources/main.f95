! ===========================================================================================================
! Code SPH 2D
!
!
! ===========================================================================================================

PROGRAM main

    use math
    use donnees
    use sph

    implicit none

    ! variables pour "maillage" initial
    integer :: n
    real(rp) :: xmin, xmax, ymin, ymax
    real(rp), dimension(:), allocatable :: x1, x2
    real(rp) :: dx, dy

    ! nombre de particules
    integer :: np
    ! coordonnées des particules
    real(rp), dimension(:, :), allocatable :: x
    ! volume des particules
    real(rp), dimension(:), allocatable :: w
    ! rayon SPH
    real(rp) :: R

    ! variables temporaires
    real(rp)               :: real1
    real(rp), dimension(2) :: real2, real2_2, real2_3
    integer, dimension(2) :: sh

    integer :: i, k

    ! vecteurs n de kappa = (1 / Rc) div(n)  (rayon de courbure) normal à la surface libre
    real(rp), dimension(:, :), allocatable :: nvec

    integer :: cpt



    ! -------------------------------------------------------------------------------------------------------
    ! lecture et initialisation des variables
    ! -------------------------------------------------------------------------------------------------------
    open(unit = 10, file = "../entrees/constantes")
    read(10, *) n
    read(10, *) xmin
    read(10, *) xmax
    read(10, *) ymin
    read(10, *) ymax
    close(10)

    ! subdivision des deux dimensions
    x1 = linspace(xmin, xmax, n)
    x2 = linspace(ymin, ymax, n)

    ! création tableau de particules et pas d'espace
    call meshgrid(x1, x2, x)
    dx = x1(2) - x1(1)
    dy = x2(2) - x2(1)

    ! nombre de particules
    sh = shape(x)
    np = sh(1)

    ! vecteur des volumes
    allocate(w(np))
    w = dx * dy

    ! rayon noyau SPH
    R = 4.0_rp * dx

    call writeMat(x, "../sorties/x.dat")



    ! -------------------------------------------------------------------------------------------------------
    ! vérification du noyau SPH et des opérateurs régularisés sur la fonction f = 1
    ! -------------------------------------------------------------------------------------------------------
    ! vérif noyau sph
    allocate(nvec(np, 4))
    write (*, '("np. |  AR de 1  |  GR de 1 (x1)   GR de 1 (x2)   |  GR_m de 1                     |  &
        &GR_p de 1")')
    write (*, '("-------------------------------------------------------------------------------------------&
        &----------------------")')
    do i = 1, np
        ! approximation régularisée
        call AR(i, x, w, R, (/ (1.0_rp, k=1, np) /), real1)
        ! gradient régularisé
        call GR(i, x, w, R, (/ (1.0_rp, k=1, np) /), real2)
        ! formattage pour gnuplot
        nvec(i, :) = (/ x(i, :), (real2 / fnorme2(real2)) * 0.5 * dx /)

        call GR_m(i, x, w, R, (/ (1.0_rp, k=1, np) /), real2_2)
        call GR_p(i, x, w, R, (/ (1.0_rp, k=1, np) /), real2_3)

        if (((37 <= i) .and. (i <= 41)) .or. &
            ((48 <= i) .and. (i <= 52)) .or. &
            ((59 <= i) .and. (i <= 63)) .or. &
            ((70 <= i) .and. (i <= 74)) .or. &
            ((81 <= i) .and. (i <= 85))) then
            !write (*, *) i, "|", real1, "|", real2
            write (*, '(1I2,"  |  ",1F7.5,"  |",2E15.5E3,"  |",2E15.5E3,"  |",2E15.5E3)') &
                i, real1, real2, real2_2, real2_3
        end if
    end do

    ! -------------------------------------------------------------------------------------------------------
    ! vérification du noyau SPH et des opérateurs régularisés sur la fonction f(x, y) = x
    ! -------------------------------------------------------------------------------------------------------
    write (*, '(///,"np. |  f(x(i))  |  AR de x  |  GR de x (x1)   GR de x (x2)   |  GR_m de x              &
        &       |  GR_p de x")')
    write (*, '("-------------------------------------------------------------------------------------------&
        &----------------------------------")')
    do i = 1, np
        ! approximation régularisée
        call AR(i, x, w, R, x(:, 1), real1)
        ! gradient régularisé
        call GR(i, x, w, R, x(:, 1), real2)
        ! formattage pour gnuplot
        nvec(i, :) = (/ x(i, :), (real2 / fnorme2(real2)) * 0.5_rp * dx /)

        call GR_m(i, x, w, R, x(:, 1), real2_2)
        call GR_p(i, x, w, R, x(:, 1), real2_3)

        if (((37 <= i) .and. (i <= 41)) .or. &
            ((48 <= i) .and. (i <= 52)) .or. &
            ((59 <= i) .and. (i <= 63)) .or. &
            ((70 <= i) .and. (i <= 74)) .or. &
            ((81 <= i) .and. (i <= 85))) then
            !write (*, *) i, "|", real1, "|", real2
            write (*, '(1I2,"  |  ",1F7.5,"  |  ",1F7.5,"  |",2E15.5E3,"  |",2E15.5E3,"  |",2E15.5E3)') &
                i, x(i, 1), real1, real2, real2_2, real2_3
        end if
    end do

    call writeMat(nvec, "../sorties/nvec.dat")



    ! -------------------------------------------------------------------------------------------------------
    ! maillage d'une bulle
    ! -------------------------------------------------------------------------------------------------------
    cpt = 0
    do i = 1, np
        if (fnorme2(x(i, :) - (/ 0.5, 0.5 /)) < 0.5) then
            !++!
        end if
    end do



    ! *******************************************************************************************************
    deallocate(x1, x2)
    deallocate(x, w)

END PROGRAM main
