module fft_poisson_test
    implicit none
    include 'fftw3.f'
    double precision, parameter :: pi = acos(-1.0d0)
    integer, parameter :: NXmin = 1, NXmax = 32     !x方向の計算領域の形状
    integer, parameter :: NYmin = 1, NYmax = 32     !y方向の計算領域の形状
    integer, parameter :: NZmin = 1, NZmax = 32     !z方向の計算領域の形状
    double precision, parameter :: Xmax = 2.0d0*pi, Ymax = 2.0d0*pi, Zmax = 2.0d0*pi  !各方向の計算領域の最大値
    double precision, parameter :: NU = 1.0d0       !動粘性係数
    double precision, parameter :: dt = 1.0d-2      !時間の刻み幅
    integer, save :: Ng                             !格子点数
    double precision, save :: ddt                   !時間の刻み幅の逆数
    double precision, save :: dX, dY, dZ            !各方向の刻み幅
    double precision, save :: ddX, ddY, ddZ         !各方向の刻み幅の逆数
    double precision, save :: ddX2, ddY2, ddZ2      !各方向の刻み幅の逆数の二乗
    double precision, save :: Xmin, Ymin, Zmin      !各方向の計算領域の最小値
    double precision, save :: plan_1d               !FFTplan
contains
!****************
!   格子を設定  *
!****************
    subroutine set_grid(X, Y, Z)
        double precision, intent(out) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        integer iX, iY, iZ
        !---各方向の刻み幅を計算---
        dX = Xmax / dble(NXmax)
        dY = Ymax / dble(NYmax)
        dZ = Zmax / dble(NZmax)
        !---計算用数値の設定---
        ddt = 1.0d0 / dt
        ddX = 1.0d0 / dX
        ddY = 1.0d0 / dY
        ddZ = 1.0d0 / dZ
        ddX2 = 1.0d0 / (dX**2)
        ddY2 = 1.0d0 / (dY**2)
        ddZ2 = 1.0d0 / (dZ**2)
        Ng = (NXmax - NXmin + 1)*(NYmax - NYmin + 1)*(NZmax - NZmin + 1)
        !---各方向の計算領域の最小値を計算---
        Xmin = dX*dble(NXmin)
        Ymin = dY*dble(NYmin)
        Zmin = dZ*dble(NZmin)
        !---格子点の設定---
        do iZ = NZmin-1, NZmax+1
            do iY = NYmin-1, NYmax+1
                do iX = NXmin-1, NXmax+1
                    X(iX, iY, iZ) = dX*dble(iX)
                    Y(iX, iY, iZ) = dY*dble(iY)
                    Z(iX, iY, iZ) = dZ*dble(iZ)
                enddo
            enddo
        enddo
    end subroutine set_grid
!************************************
!   各格子点における予測速度の計算  *
!************************************
    subroutine set_velocity(Vx_p, Vy_p, Vz_p, X, Y, Z)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Vy_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Vz_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        integer iX, iY, iZ
        do iZ = NZmin-1, NZmax
            do iY = NYmin-1, NYmax
                do iX = NXmin-1, NXmax
                    Vx_p(iX, iY, iZ) = sin((X(iX,iY,iZ) + 0.5d0*dX))
                    Vy_p(iX, iY, iZ) = sin((Y(iX,iY,iZ) + 0.5d0*dY))
                    Vz_p(iX, iY, iZ) = sin((Z(iX,iY,iZ) + 0.5d0*dZ))
                enddo
            enddo
        enddo
    end subroutine set_velocity
!************************************
!   ポアソン方程式の解析解の計算    *
!************************************
    subroutine cal_theory(Phi_th, X, Y, Z)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Phi_th(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        integer iX, iY, iZ
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                do iX = NXmin, NXmax
                    Phi_th(iX,iY,iZ) = -ddt*(cos(X(iX,iY,iZ)) + cos(Y(iX,iY,iZ)) + cos(Z(iX,iY,iZ)))
                enddo
            enddo
        enddo
    end subroutine cal_theory
!********************************
!   ポアソン方程式の右辺を計算  *
!********************************
    subroutine cal_poisson_rhs(Vx_p, Vy_p, Vz_p, divU)
        double precision, intent(in) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: divU(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        integer iX, iY, iZ
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                do iX = NXmin, NXmax
                    divU(iX,iY,iZ) = ddt*(ddX*(-Vx_p(iX-1,iY,iZ) + Vx_p(iX,iY,iZ)) &
                                    + ddY*(-Vy_p(iX,iY-1,iZ) + Vy_p(iX,iY,iZ)) &
                                    + ddZ*(-Vz_p(iX,iY,iZ-1) + Vz_p(iX,iY,iZ)))
                enddo
            enddo
        enddo
    end subroutine cal_poisson_rhs
!****************
!   plan生成    *
!****************
    subroutine FFT_init
        double precision AL(NXmin:NXmax)
        double complex BL(NXmin-1:NXmax/2)
        call dfftw_plan_dft_r2c_1d(plan_1d, NXmax, AL, BL, FFTW_ESTIMATE)
    end subroutine FFT_init
!************************
!   FFT(1次元)を実行    *
!************************
    subroutine FFT_1d_exe(RHS, RHSf)
        double precision, intent(in) :: RHS(NXmin:NXmax)
        double complex, intent(out) :: RHSf(NXmin-1:NXmax/2)
        call dfftw_execute(plan_1d)
    end subroutine FFT_1d_exe
!***************************************
!   ポアソン方程式の右辺をFFT(x方向)   *
!***************************************
    subroutine FFT_rhs(divU, divUf)
        double precision, intent(in) :: divU(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        double complex, intent(out) :: divUf(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        double complex RHSf(NXmin-1:NXmax/2)
        integer iY, iZ
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                call FFT_1d_exe(divU(NXmin:NXmax, iY, iZ), RHSf)
                divUf(NXmin-1:NXmax/2, iY, iZ) = RHSf(NXmin-1:NXmax/2)
            enddo
        enddo
    end subroutine FFT_rhs
!********************************
!   ポアソン方程式の反復計算    *
!********************************
    subroutine itr_poisson(divUf, Phif)
        double complex, intent(in) :: divUf(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        double complex, intent(out) :: Phif(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        double complex Er(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        double precision B0
        integer kX, iY, iZ, itr, itrmax
        itrmax = 10000
        do kX = NXmin-1, NXmax/2
            do itr = 1, itrmax
                do iZ = NZmin, NZmax
                    do iY = NYmin, NYmax
                    enddo
                enddo
            enddo
        enddo
    end subroutine itr_poisson
!****************************************
!   ポアソン方程式を解くサブルーチン    *
!****************************************
    subroutine cal_poisson(Vx_p, Vy_p, Vz_p, Phi, X, Y, Z)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision divU(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        double complex divUf(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        double complex Phif(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        !---右辺を計算---
        call cal_poisson_rhs(Vx_p, Vy_p, Vz_p, divU)
        !---右辺をFFT---
        call FFT_rhs(divU, divUf)
    end subroutine cal_poisson
!********************************
!   誤差を計算するサブルーチン  *
!********************************
    subroutine cal_error(Phi, Phi_th)
        double precision, intent(in) :: Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Phi_th(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        double precision Phi_num(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        double precision error(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax), error_sum, error_norm
        double precision C(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax), C_ave
        integer iX, iY, iZ
        !---積分定数の計算---
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                do iX = NXmin, NXmax
                    C(iX,iY,iZ) = Phi(iX,iY,iZ) - Phi_th(iX,iY,iZ)
                    C_ave = C_ave + C(iX,iY,iZ)
                enddo
            enddo
        enddo
        C_ave = C_ave / Ng
        !---誤差の計算---
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                do iX = NXmin, NXmax
                    Phi_num(iX,iY,iZ) = Phi(iX,iY,iZ) - C_ave
                    error(iX,iY,iZ) = Phi_num(iX,iY,iZ) - Phi_th(iX,iY,iZ)
                    error_sum = error_sum + error(iX,iY,iZ)**2
                enddo
            enddo
        enddo
        error_norm = sqrt(error_sum / Ng)
        open(11, file = 'chk_poisson_precision.dat', position = 'append')
        write(11, *) dX, error_norm
        close(11)
    end subroutine cal_error
end module fft_poisson_test

!************************
!   メインプログラム    *
!************************
program main
    use fft_poisson_test
    implicit none
    double precision X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Phi_th(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
    double precision Vx_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Vy_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Vz_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    call set_grid(X, Y, Z)
    call FFT_init
    call set_velocity(Vx_p, Vy_p, Vz_p, X, Y, Z)
    call cal_theory(Phi_th, X, Y, Z)
    call cal_poisson(Vx_p, Vy_p, Vz_p, Phi, X, Y, Z)
    call cal_error(Phi, Phi_th)
end program main


