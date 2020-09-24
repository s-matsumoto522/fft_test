module fft_2dpoisson
    implicit none
    include 'fftw3.f'
    double precision, parameter :: pi = acos(-1.d0)
    integer, parameter :: NXmin = 1, NXmax = 32     !x方向の計算形状
    integer, parameter :: NYmin = 1, NYmax = 32     !y方向の計算形状
    double precision, parameter :: Xmax = 2.0d0*pi, Ymax = 2.0d0*pi     !各方向の計算領域の最大値
    double precision, parameter :: dt = 1.0d0       !時間刻み幅
    integer, save :: Ng                             !総格子点数
    double precision, save :: dX, dY                !各方向の刻み幅
    double precision, save :: ddt                   !計算用数値
    double precision, save :: ddX, ddY
    double precision, save :: ddX2, ddY2
    double precision, save :: Xmin, Ymin            !各方向の計算領域の最小値
    double precision, save :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1)   !X座標
    double precision, save :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1)   !Y座標
contains
!****************
!   格子を設定  *
!****************
    subroutine set_grid
        integer iX, iY
        !---各方向の刻み幅を計算---
        dX = Xmax / dble(NXmax)
        dY = Ymax / dble(NYmax)
        !---計算用数値の設定---
        ddt = 1.0d0 / dt
        ddX = 1.0d0 / dX
        ddY = 1.0d0 / dY
        ddX2 = 1.0d0 / (dX**2)
        ddY2 = 1.0d0 / (dY**2)
        Ng = (NXmax - NXmin + 1)*(NYmax - NYmin + 1)
        !---各方向の計算領域の最小値を計算---
        Xmin = 0.0d0
        Ymin = 0.0d0
        !---格子点の設定---
        do iY = NYmin-1, NYmax+1
            do iX = NXmin-1, NXmax+1
                X(iX, iY) = Xmin + dX*dble(iX - 1)
                Y(iX, iY) = Ymin + dY*dble(iY - 1)
            enddo
        enddo
    end subroutine set_grid
!************************************
!   各格子点における予測速度の計算  *
!************************************
    subroutine set_velocity(Vx_p, Vy_p)
        double precision, intent(out) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax)
        double precision, intent(out) :: Vy_p(NXmin-1:NXmax, NYmin-1:NYmax)
        integer iX, iY
        do iY = NYmin-1, NYmax
            do iX = NXmin-1, NXmax
                Vx_p(iX, iY) = sin(X(iX, iY) + 0.5d0*dX)
                Vy_p(iX, iY) = sin(Y(iX, iY) + 0.5d0*dY)
            enddo
        enddo
    end subroutine set_velocity
!************************************
!   ポアソン方程式の解析解を計算    *
!************************************
    subroutine cal_theory(Phi_th)
        double precision, intent(out) :: Phi_th(NXmin:NXmax, NYmin:NYmax)
        integer iX, iY
        do iY = NYmin, NYmax
            do iX = NXmin, NXmax
                Phi_th(iX, iY) = -ddt*(cos(X(iX, iY)) + cos(Y(iX, iY)))
            enddo
        enddo
    end subroutine cal_theory
!********************************
!   ポアソン方程式の右辺を計算  *
!********************************
    subroutine cal_poisson_rhs(Vx_p, Vy_p, divU)
        double precision, intent(in) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax)
        double precision, intent(in) :: Vy_p(NXmin-1:NXmax, NYmin-1:NYmax)
        double precision, intent(out) :: divU(NXmin:NXmax, NYmin:NYmax)
        integer iX, iY
        do iY = NYmin, NYmax
            do iX = NXmin, NXmax
                divU(iX, iY) = ddt*(ddX*(-Vx_p(iX-1, iY) + Vx_p(iX, iY)) &
                                    + ddY*(-Vy_p(iX, iY-1) + Vy_p(iX, iY)))
            enddo
        enddo
    end subroutine cal_poisson_rhs
!************************
!   FFT(1次元)を実行    *
!************************
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
!****************************************
!   ポアソン方程式の右辺のFFT(x方向)    *
!****************************************
    subroutine FFT_rhs(divU, divUf)
        double precision, intent(in) :: divU(NXmin:NXmax, NYmin:NYmax)
        complex(kind(0d0)), intent(out) :: divUf(0:NXmax/2, NYmin:NYmax)
        double precision RHS(NXmin:NXmax)
        complex(kind(0d0)) RHSf(0:NXmax/2)
        integer iY, kX
        open(21, file = 'debug_divUf(2d).dat')
        open(22, file = 'debug_RHS(2d).dat')
        do iY = NYmin, NYmax
            RHS(NXmin:NXmax) = divU(NXmin:NXmax, iY)
            do kX = NXmin, NXmax
                write(22, *) RHS(kX)
            enddo
            write(22, *) ''
            call FFT_1d_exe(RHS, RHSf)
            divUf(0:NXmax/2, iY) = RHSf(0:NXmax/2)
            do kX = 0, NXmax/2
                write(21, *) divUf(kX, iY)
            enddo
            write(21, *) ''
        enddo
        close(21)
        close(22)
    end subroutine FFT_rhs
!****************************************
!   スカラーポテンシャルをIFFT(x方向)   *
!****************************************
    subroutine IFFT_Phif(Phif, Phi)
        complex(kind(0d0)), intent(in) :: Phif(0:Nxmax/2, NYmin-1:NYmax+1)
        double precision, intent(out) :: Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1)
        double precision P(NXmin:NXmax)
        complex(kind(0d0)) Pf(0:NXmax/2)
        integer iY, kX
        open(31, file = 'debug_Phif(2d).dat')
        open(32, file = 'debuf_Phi(2d).dat')
        do iY = NYmin, NYmax
            do kX = 0, NXmax/2
                write(31, *) Phif(kX, iY)
            enddo
            write(31, *) ''
            Pf(0:NXmax/2) = Phif(0:NXmax/2, iY)
            call IFFT_1d_exe(Pf, P)
            Phi(NXmin:NXmax, iY) = P(NXmin:NXmax)
            do kX = NXmin, NXmax
                write(32, *) Phi(kX, iY)
            enddo
            write(32, *)''
        enddo
        close(31)
        close(32)
    end subroutine IFFT_Phif
!****************************************
!   反復計算中の境界条件を設定(y方向)   *
!****************************************
    subroutine set_bc_itr(Phif)
        complex(kind(0d0)), intent(out) :: Phif(0:NXmax/2, NYmin-1:NYmax+1)
        Phif(0:NXmax/2, NYmin-1) = Phif(0:NXmax/2, NYmin)   !ノイマン条件
        Phif(0:NXmax/2, NYmax+1) = Phif(0:NXmax/2, NYmax)
    end subroutine set_bc_itr
!********************************************
!   スカラーポテンシャルの境界条件を設定    *
!********************************************
    subroutine set_bc(Phi)
        double precision, intent(out) :: Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1)
        Phi(NXmax+1, NYmin-1:NYmax+1) = Phi(NXmin, NYmin-1:NYmax+1)     !周期境界
        Phi(NXmin-1, NYmin-1:NYmax+1) = Phi(NXmax, NYmin-1:NYmax+1)     !周期境界
        Phi(NXmin:NXmax, NYmin-1) = Phi(NXmin:NXmax, NYmin)             !ノイマン条件
        Phi(NXmin:NXmax, NYmax+1) = Phi(NXmin:NXmax, NYmax)             !ノイマン条件
    end subroutine set_bc
!********************************
!   ポアソン方程式の反復計算    *
!********************************
    subroutine itr_poisson(divUf, Phif)
        complex(kind(0d0)), intent(in) :: divUf(0:NXmax/2, NYmin:NYmax)
        complex(kind(0d0)), intent(out) :: Phif(0:NXmax/2, NYmin-1:NYmax+1)
        complex(kind(0d0)) Er(0:NXmax/2, NYmin:NYmax)
        double precision B0, dB0, norm_er, norm_rhs, error
        integer kX, iY, itr, itrmax
        double precision beta, eps
        beta = 1.7d0        !過緩和係数
        itrmax = 1.0d5      !最大反復回数
        eps = 1.0d-6        !誤差の閾値
        Phif(:, :) = (0.0d0, 0.0d0)     !初期値の設定
        !---各波数毎にポアソン方程式をSOR法で解く----
        do kX = 0, NXmax/2
            !---計算用係数の計算---
            B0 = 2.0d0*(ddX2*(1.0d0 - cos(dble(kX)*dX)) + ddY2)
            dB0 = 1.0d0 / B0
            do itr = 1, itrmax
                !---誤差の初期値の設定---
                error = 0.0d0
                norm_er = 0.0d0
                norm_rhs = 0.0d0
                do iY = NYmin, NYmax
                    Er(kX,iY) = ddY2*(Phif(kX,iY-1) + Phif(kX,iY+1)) - divUf(kX,iY) - B0*Phif(kX,iY)
                    Phif(kX,iY) = Phif(kX,iY) + beta*dB0*Er(kX,iY)      !解の更新
                    norm_er = norm_er + sqrt(real(Er(kX,iY))**2 + aimag(Er(kX,iY))**2)
                    norm_rhs = norm_rhs + sqrt(real(divUf(kX,iY))**2 + aimag(divUf(kX,iY))**2)
                enddo
                norm_er = sqrt(norm_er / dble(NYmax))
                norm_rhs = sqrt(norm_rhs / dble(NYmax))
                error = norm_er / norm_rhs
                if(mod(itr, 100) == 0) then
                    write(*, *) 'itr, error = ', itr, error
                endif
                !---境界条件の設定---
                call set_bc_itr(Phif)
                !---収束判定---
                if(error < eps) then
                    write(*, *) 'poisson is converged.'
                    exit
                endif
            enddo
        enddo
    end subroutine itr_poisson
!****************************************
!   ポアソン方程式を解くサブルーチン    *
!****************************************
    subroutine cal_poisson(Vx_p, Vy_p, Phi)
        double precision, intent(in) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax)
        double precision, intent(in) :: Vy_p(NXmin-1:NXmax, NYmin-1:NYmax)
        double precision, intent(out) :: Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1)
        double precision divU(NXmin:NXmax, NYmin:NYmax)
        complex(kind(0d0)) divUf(0:NXmax/2, NYmin:NYmax)
        complex(kind(0d0)) Phif(0:NXmax/2, NYmin-1:NYmax+1)
        !---ポアソン方程式の右辺を計算---
        call cal_poisson_rhs(Vx_p, Vy_p, divU)
        !---右辺をFFT---
        call FFT_rhs(divU, divUf)
        !---反復計算---
        call itr_poisson(divUf, Phif)
        !---Phifを逆フーリエ変換---
        call IFFT_Phif(Phif, Phi)
        !---Phiの境界条件を設定---
        call set_bc(Phi)
    end subroutine cal_poisson
!********************************
!   誤差を計算するサブルーチン  *
!********************************
    subroutine cal_error(Phi, Phi_th)
        double precision, intent(in) :: Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1)
        double precision, intent(in) :: Phi_th(NXmin:NXmax, NYmin:NYmax)
        double precision Phi_num(NXmin:NXmax, NYmin:NYmax)
        double precision error(NXmin:NXmax, NYmin:NYmax), error_sum, norm_error
        double precision C(NXmin:NXmax, NYmin:NYmax), C_ave
        integer iX, iY
        !---積分定数の計算---
        C_ave = 0.0d0
        do iY = NYmin, NYmax
            do iX = NXmin, NXmax
                C(iX, iY) = Phi(iX, iY) - Phi_th(iX, iY)
                C_ave = C_ave + C(iX, iY)
            enddo
        enddo
        C_ave = C_ave / dble(Ng)
        write(*, *) C_ave
        !---誤差の計算---
        error_sum = 0.0d0
        do iY = NYmin, NYmax
            do iX = NXmin, NXmax
                Phi_num(iX, iY) = Phi(iX, iY) - C_ave
                error(iX, iY) = Phi_num(iX, iY) - Phi_th(iX, iY)
                error_sum = error_sum + error(iX, iY)**2
            enddo
        enddo
        norm_error = sqrt(error_sum / dble(Ng))
        open(11, file = 'chk_poisson2d_precision.dat', position = 'append')
        write(11, *) dX, norm_error
        close(11)
    end subroutine cal_error
end module fft_2dpoisson

!************************
!   メインプログラム    *
!************************
program main
    use fft_2dpoisson
    implicit none
    double precision Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1)
    double precision Phi_th(NXmin:NXmax, NYmin:NYmax)
    double precision Vx_p(NXmin-1:NXmax, NYmin-1:NYmax)
    double precision Vy_p(NXmin-1:NXmax, NYmin-1:NYmax)
    call set_grid
    call set_velocity(Vx_p, Vy_p)
    call cal_theory(Phi_th)
    call cal_poisson(Vx_p, Vy_p, Phi)
    call cal_error(Phi, Phi_th)
end program main