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

    ! PARTIE 1 : MAILLAGE D'UN PAVÉ 2D
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
    real(rp), dimension(:, :), allocatable :: nvec, grad_P, plot_vec


    ! PARTIE 2 : MAILLAGE D'UNE BULLE
    integer :: cpt
    real(rp) :: rayon
    real(rp), dimension(2) :: centre
    real(rp), dimension(100, 2) :: cc
    real(rp), dimension(100) ::xx

    ! VARIABLES ÉQUATION
    real(rp), dimension(:), allocatable :: P
    real(rp), dimension(:, :), allocatable :: sol, nor, fts, d_rwu_dt
    real(rp), dimension(2) :: temp



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
    call writeMat(x, "../sorties/x.dat")
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
        nvec(i, :) = (/ x(i, :), (real2 / fnorme2(real2)) * 0.5_rp * dx /)

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

    call writeMat(nvec, "../sorties/nvec.dat")



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



    ! -------------------------------------------------------------------------------------------------------
    ! maillage d'une bulle
    ! -------------------------------------------------------------------------------------------------------
    deallocate(x)

    ! données du disque à mailler
    centre = (/ 1.0_rp, 1.0_rp /)
    rayon = 1.0_rp


    ! maillage du disque et actualisation des paramètre dépendant du maillage
    call meshCircle(centre, rayon, n, x)
    dx = x(2, 1) - x(1, 1)
    R = 4.0_rp * dx
    call writeMat(x, "../sorties/x_circle.dat")
    sh = shape(x)
    np = sh(1)


    ! volume des particules
    deallocate(w, nvec)
    allocate(w(np), nvec(np, 4))
    !w = pi * rayon**2 / real(np, rp)
    w = dx**2


    ! vérification sur l'approx régularisée de la fonction 1 et de son gradient (0, 0)
    i = nint(np / 2.0_rp)
    call AR(i, x, w, R, (/ (1.0_rp, k=1, np) /), real1)
    call GR(i, x, w, R, (/ (1.0_rp, k=1, np) /), real2)
    print *, real1, real2


    ! contour de la bulle que nous avons maillé pour graphique
    xx = linspace(0.0_rp, 2.0_rp * pi, 100)
    cc(:, 1) = centre(1) + rayon * cos(xx)
    cc(:, 2) = centre(2) + rayon * sin(xx)
    call writeMat(cc, "../sorties/cc.dat")


    ! initialisation de la pression
    allocate(P(np))
    call init_pression(x, centre, rayon, P)
    allocate(sol(np, 3))
    sol(:, 1:2) = x
    sol(:, 3) = P
    call writeMat(sol, "../sorties/P.dat")


    ! approximation du gradient de pression avec l'opérateur GR_p
    allocate(grad_P(np, 4))
    do i = 1, np
        call GR_p(i, x, w, R, P, real2)
        ! formattage pour gnuplot
        grad_P(i, :) = (/ x(i, :), (real2 / fnorme2(real2)) * 0.5_rp * dx /)
    end do
    call writeMat(grad_P, "../sorties/grad_P.dat")


    ! tension de surface
    allocate(nor(np, 2))
    allocate(fts(np, 2))
    call normale(R, x, w, nor)
    do i = 1, np
        call F_TS(gamma_eau_air, i, x, w, nor, R, fts(i, :))
        nvec(i, :) = (/ x(i, :), (fts(i, :) / fnorme2(fts(i, :))) * 0.5_rp * dx /)
        !nvec(i, :) = (/ x(i, :), (nor(i, :) / fnorme2(nor(i, :))) * 0.5_rp * dx /)
    end do
    call writeMat(nvec, "../sorties/n_fts.dat")
    print *, "x"
    print *, fts(:, 1) - grad_P(:, 3)
    print *, "y"
    print *, fts(:, 2) - grad_P(:, 4)

    allocate(d_rwu_dt(np, 2))
    d_rwu_dt(:, 1) = w * grad_P(:, 3)
    d_rwu_dt(:, 2) = w * grad_P(:, 4)

    d_rwu_dt = d_rwu_dt + fts



    ! *******************************************************************************************************
    deallocate(x1, x2)
    deallocate(x, w)
    deallocate(grad_P)
    deallocate(fts)
    deallocate(nor)
    deallocate(d_rwu_dt)

END PROGRAM main
