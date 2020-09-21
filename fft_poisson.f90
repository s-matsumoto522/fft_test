module fft_poisson_test
    implicit none
    include 'fftw3.f'
    double precision, parameter :: pi = acos(-1.0d0)
    integer, parameter :: NXmin = 1, NXmax = 64     !x方向の計算領域の形状
    integer, parameter :: NYmin = 1, NYmax = 64     !y方向の計算領域の形状
    integer, parameter :: NZmin = 1, NZmax = 64     !z方向の計算領域の形状
    double precision, parameter :: Xmax = 2.0d0*pi, Ymax = 2.0d0*pi, Zmax = 2.0d0*pi  !各方向の計算領域の最大値
    double precision, parameter :: NU = 1.0d0       !動粘性係数
    double precision, parameter :: dt = 1.0d0      !時間の刻み幅
    integer, save :: Ng                             !全格子点数
    double precision, save :: ddt                   !時間の刻み幅の逆数
    double precision, save :: dX, dY, dZ            !各方向の刻み幅
    double precision, save :: ddX, ddY, ddZ         !各方向の刻み幅の逆数
    double precision, save :: ddX2, ddY2, ddZ2      !各方向の刻み幅の逆数の二乗
    double precision, save :: Xmin, Ymin, Zmin      !各方向の計算領域の最小値
    double precision, save :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)  !X座標
    double precision, save :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)  !Y座標
    double precision, save :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)  !Z座標
contains
!****************
!   格子を設定  *
!****************
    subroutine set_grid
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
    subroutine set_velocity(Vx_p, Vy_p, Vz_p)
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
    subroutine cal_theory(Phi_th)
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
!************************
!   FFT(1次元)を実行    *
!************************
    subroutine FFT_1d_exe(RHS, RHSf)
        double precision, intent(in) :: RHS(NXmin:NXmax)
        complex(kind(0d0)), intent(out) :: RHSf(NXmin-1:NXmax/2)
        double precision plan_1d_fft
        call dfftw_plan_dft_r2c_1d(plan_1d_fft, NXmax, RHS, RHSf, FFTW_ESTIMATE)
        call dfftw_execute(plan_1d_fft, RHS, RHSf)
        call dfftw_destroy_plan(plan_1d_fft)
    end subroutine FFT_1d_exe
!************************
!   IFFT(1次元)を実行   *
!************************
    subroutine IFFT_1d_exe(Pf, P)
        complex(kind(0d0)), intent(in) :: Pf(NXmin-1:NXmax/2)
        double precision, intent(out) :: P(NXmin:NXmax)
        double precision plan_1d_ifft
        call dfftw_plan_dft_c2r_1d(plan_1d_ifft, NXmax, Pf, P, FFTW_ESTIMATE)
        call dfftw_execute(plan_1d_ifft, Pf, P)
        call dfftw_destroy_plan(plan_1d_ifft)
    end subroutine IFFT_1d_exe
!***************************************
!   ポアソン方程式の右辺をFFT(x方向)   *
!***************************************
    subroutine FFT_rhs(divU, divUf)
        double precision, intent(in) :: divU(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        complex(kind(0d0)), intent(out) :: divUf(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        double precision RHS(NXmin:NXmax)
        complex(kind(0d0)) RHSf(NXmin-1:NXmax/2)
        integer iY, iZ, kX
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                RHS(NXmin:NXmax) = divU(NXmin:NXmax, iY, iZ)
                call FFT_1d_exe(RHS, RHSf)
                divUf(NXmin-1:NXmax/2, iY, iZ) = RHSf(NXmin-1:NXmax/2)
            enddo
        enddo
    end subroutine FFT_rhs
!****************************************
!   スカラーポテンシャルをIFFT(x方向)   *
!****************************************
    subroutine IFFT_Phif(Phif, Phi)
        complex(kind(0d0)), intent(in) :: Phif(NXmin-1:NXmax/2, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision P(NXmin:NXmax)
        complex(kind(0d0)) Pf(NXmin-1:NXmax/2)
        integer iY, iZ
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                Pf(NXmin-1:NXmax/2) = Phif(NXmin-1:NXmax/2, iY, iZ)
                call IFFT_1d_exe(Pf, P)
                Phi(NXmin:NXmax, iY, iZ) = P(NXmin:NXmax)
            enddo
        enddo
        !---逆フーリエ変換で得た値を規格化---
        Phi(NXmin:NXmax,NYmin-1:NYmax+1,NZmin-1:NZmax+1) = Phi(NXmin:NXmax,NYmin-1:NYmax+1,NZmin-1:NZmax+1)/dble(NXmax+1)
    end subroutine IFFT_Phif
!********************************************
!   反復計算中の境界条件を設定(y, z方向)    *
!********************************************
    subroutine set_bc_itr(Phif)
        complex(kind(0d0)), intent(out) :: Phif(NXmin-1:NXmax/2, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        Phif(NXmin-1:NXmax/2,NYmin-1,NZmin:NZmax) = Phif(NXmin-1:NXmax/2,NYmin,NZmin:NZmax)
        Phif(NXmin-1:Nxmax/2,NYmax+1,NZmin:NZmax) = Phif(NXmin-1:NXmax/2,NYmax,NZmin:NZmax)
        Phif(NXmin-1:NXmax/2,NYmin:NYmax,NZmin-1) = Phif(NXmin-1:Nxmax/2,NYmin:NYmax,NZmin)
        Phif(NXmin-1:NXmax/2,NYmin:NYmax,NZmax+1) = Phif(NXmin-1:NXmax/2,NYmin:NYmax,NZmax)
    end subroutine set_bc_itr
!********************************************
!   スカラーポテンシャルの境界条件を設定    *
!********************************************
    subroutine set_bc(Phi)
        double precision, intent(out) :: Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        Phi(NXmax+1,NYmin:NYmax,NZmin:NZmax) = Phi(NXmin-1,NYmin:NYmax,NZmin:NZmax)         !周期境界
        Phi(NXmin:NXmax,NYmin-1,NZmin:NZmax) = Phi(NXmin:NXmax,NYmin,NZmin:NZmax)           !ノイマン条件
        Phi(NXmin:NXmax,NYmax+1,NZmin:NZmax) = Phi(NXmin:NXmax,NYmax,NZmin:NZmax)
        Phi(NXmin:NXmax,NYmin:NYmax,NZmin-1) = Phi(NXmin:NXmax,NYmin:NYmax,NZmin)
        Phi(NXmin:NXmax,NYmin:NYmax,NZmax+1) = Phi(NXmin:NXmax,NYmin:NYmax,NZmax)
    end subroutine set_bc
!********************************
!   ポアソン方程式の反復計算    *
!********************************
    subroutine itr_poisson(divUf, Phif)
        complex(kind(0d0)), intent(in) :: divUf(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        complex(kind(0d0)), intent(out) :: Phif(NXmin-1:NXmax/2, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        complex(kind(0d0)) Er(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        double precision B0, dB0, norm_er, norm_rhs, error
        integer kX, iY, iZ, itr, itrmax, N_2dg
        double precision beta, eps
        N_2dg = (NYmax - NYmin + 1)*(NZmax - NZmin + 1)     !y, z方向のみの格子点数
        beta = 1.7d0            !過緩和係数
        itrmax = 1.0d5          !最大反復回数
        eps = 1.0d-5            !誤差の閾値
        Phif(:, :, :) = (0.0d0, 0.0d0)   !初期値の計算
        !---各波数ごとのポアソン方程式をSOR法で解く---
        do kX = NXmin-1, NXmax/2
            B0 = 2.0d0*(ddX2*(1.0d0 - cos(dble(kX)*dX)) + ddY2 + ddZ2)
            dB0 = 1.0d0 / B0
            do itr = 1, itrmax
                !---誤差の初期値---
                error = 0.0d0
                norm_er = 0.0d0
                norm_rhs = 0.0d0
                do iZ = NZmin, NZmax
                    do iY = NYmin, NYmax
                        Er(kX,iY,iZ) = ddY2*(Phif(kX,iY-1,iZ) + Phif(kX,iY+1,iZ)) &
                                    + ddZ2*(Phif(kX,iY,iZ-1) + Phif(kX,iY,iZ+1)) &
                                    -divUf(kX,iY,iZ) - B0*Phif(kX,iY,iZ)
                        Phif(kX,iY,iZ) = Phif(kX,iY,iZ) + beta*dB0*Er(kX,iY,iZ)  !解の更新
                        norm_er = norm_er + sqrt(real(Er(kX,iY,iZ))**2 + aimag(Er(kX,iY,iZ))**2)
                        norm_rhs = norm_rhs + sqrt(real(divUf(kX,iY,iZ))**2 + aimag(divUf(kX,iY,iZ))**2)
                    enddo
                enddo
                norm_er = sqrt(norm_er / dble(N_2dg))
                norm_rhs = sqrt(norm_rhs / dble(N_2dg))
                error = norm_er / norm_rhs
                if(mod(itr, 100) == 0) then
                    write(*, *) 'itr, error = ', itr, error
                endif
                !---境界条件の設定---
                call set_bc_itr(Phif)
                !---収束判定---
                if(error < eps) then
                    write(*, *) 'converged'
                    exit
                endif
            enddo
        enddo
    end subroutine itr_poisson
!****************************************
!   ポアソン方程式を解くサブルーチン    *
!****************************************
    subroutine cal_poisson(Vx_p, Vy_p, Vz_p, Phi)
        double precision, intent(in) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision divU(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
        complex(kind(0d0)) divUf(NXmin-1:NXmax/2, NYmin:NYmax, NZmin:NZmax)
        complex(kind(0d0)) Phif(NXmin-1:NXmax/2, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        !---右辺を計算---
        call cal_poisson_rhs(Vx_p, Vy_p, Vz_p, divU)
        !---右辺をFFT---
        call FFT_rhs(divU, divUf)
        !---反復計算---
        call itr_poisson(divUf, Phif)
        !---逆フーリエ変換---
        call IFFT_Phif(Phif, Phi)
        !---境界条件を設定---
        call set_bc(Phi)
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
        C_ave = 0.0d0
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                do iX = NXmin, NXmax
                    C(iX,iY,iZ) = Phi(iX,iY,iZ) - Phi_th(iX,iY,iZ)
                    C_ave = C_ave + C(iX,iY,iZ)
                enddo
            enddo
        enddo
        C_ave = C_ave / dble(Ng)
        write(*, *) C_ave
        !---誤差の計算---
        error_sum = 0.0d0
        do iZ = NZmin, NZmax
            do iY = NYmin, NYmax
                do iX = NXmin, NXmax
                    Phi_num(iX,iY,iZ) = Phi(iX,iY,iZ) - C_ave
                    error(iX,iY,iZ) = Phi_num(iX,iY,iZ) - Phi_th(iX,iY,iZ)
                    error_sum = error_sum + error(iX,iY,iZ)**2
                enddo
            enddo
        enddo
        error_norm = sqrt(error_sum / dble(Ng))
        open(11, file = 'chk_fft_poisson_precision.dat', position = 'append')
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
    double precision Phi(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Phi_th(NXmin:NXmax, NYmin:NYmax, NZmin:NZmax)
    double precision Vx_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Vy_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Vz_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    call set_grid
    call set_velocity(Vx_p, Vy_p, Vz_p)
    call cal_theory(Phi_th)
    call cal_poisson(Vx_p, Vy_p, Vz_p, Phi)
    call cal_error(Phi, Phi_th)
end program main


