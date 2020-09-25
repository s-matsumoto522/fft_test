module fft_1dpoisson
    implicit none
    include 'fftw3.f'
    double precision, parameter :: pi = acos(-1.0d0)
    integer, parameter :: Nxmin = 1, NXmax = 32     !x方向の計算領域
    double precision, parameter :: Xmax = 2.0d0*pi  !x方向の計算領域の最大値
    double precision, parameter :: dt = -1.0d0      !時間刻み幅
    integer, save :: Ng                             !全格子点数
    double precision, save :: dX                    !x方向の刻み幅
    double precision, save :: ddt                   !計算用数値
    double precision, save :: ddX
    double precision, save :: ddX2
    double precision, save :: Xmin                  !x方向の計算領域の最小値
    double precision, save :: X(NXmin-1:NXmax+1)    !x座標
contains
!***************
!   格子の設定  *
!***************
    subroutine set_grid
        integer iX
        !---各方向の刻み幅を計算---
        dX = Xmax / dble(NXmax)
        !---計算用数値の設定---
        ddt = 1.0d0 / dt
        ddX = 1.0d0 / dX
        ddX2 = 1.0d0 / (dX**2)
        Ng = (NXmax - NXmin + 1)
        !---各方向の計算領域の最小値の計算---
        Xmin = 0.0d0
        !---格子点の設定---
        do iX = NXmin-1, NXmax+1
            X(iX) = Xmin + dX*dble(iX - 1)
        enddo
    end subroutine set_grid
!*********************************
!   各格子点における予測速度の計算  *
!*********************************
    subroutine set_velocity(Vx_p)
        double precision, intent(in) :: Vx_p(NXmin-1:NXmax+1)
        integer iX
        do iX = NXmin-1, NXmax
            Vx_p(iX) = sin(X(iX) + 0.5d0*dX)
        enddo
    end subroutine set_velocity
!*********************************
!   ポアソン方程式の解析解を計算    *
!*********************************
    subroutine cal_theory(Phi_th)
        double precision, intent(out) :: Phi_th(NXmin:NXmax)
        integer iX
        do iX = NXmin, NXmax
            Phi_th(iX) = -ddt*cos(X(iX))
        enddo
    end subroutine cal_theory
!*****************************
!   ポアソン方程式の右辺を計算  *
!*****************************
    subroutine cal_poisson_rhs(Vx_p, divU)
        double precision, intent(in) :: Vx_p(NXmin-1:NXmax)
        double precision, intent(out) :: divU(NXmin:NXmax)
        integer iX
        do iX = NXmin, NXmax
            divU(iX) = ddt*(cos(X(iX)))
        enddo
    end subroutine cal_poisson_rhs
!***********************
!   FFT(1次元)を実行    *
!***********************
    subroutine FFT_1d_exe(RHS, RHSf)
        double precision, intent(in) :: RHS(NXmin:NXmax)
        complex(kind(0d0)), intent(out) :: RHSf(0:NXmax/2)
        double precision plan_1d_fft
        call dfftw_plan_dft_r2c_1d(plan_1d_fft, NXmax, RHS, RHSf, FFTW_ESTIMATE)
        call dfftw_execute(plan_1d_fft, RHS, RHSf)
        call dfftw_destroy_plan(plan_1d_fft)
    end subroutine FFT_1d_exe
!************************
!   IFFT(1次元)を実行   *
!************************
    subroutine IFFT_1d_exe(Pf, P)
        complex(kind(0d0)), intent(in) :: Pf(0:NXmax/2)
        double precision, intent(out) :: P(NXmin:NXmax)
        double precision plan_1d_ifft
        call dfftw_plan_dft_c2r_1d(plan_1d_ifft, NXmax, Pf, P, FFTW_ESTIMATE)
        call dfftw_execute(plan_1d_ifft, Pf, P)
        call dfftw_destroy_plan(plan_1d_ifft)
        !---規格化---
        P(:) = P(:) / dble(NXmax)
    end subroutine IFFT_1d_exe
!
!
!