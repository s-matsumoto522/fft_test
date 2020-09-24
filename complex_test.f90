program main
    implicit none
    !double precision norm
    double complex im, im2
    im = (3.0d0, 2.0d0)
    im2 = (2.0d0, 1.0d0)
    write(*, *) im*2
    !im = 2.0d0*im
    !norm = real(im)**2 + aimag(im)**2
    !write(*, *) norm
end program main