module fmesh_wrapper

	use iso_c_binding, only: c_double, c_int
	use mesh, only: ptc2pdf, pdfgrad

	implicit none

contains

	subroutine c_ptc2pdf(Ng,dx,N,pg,x,y,spwt_input) bind(c)
		integer(c_int), intent(in) :: Ng, N
		real(c_double), intent(in) :: dx, x(N), y(N)
		real(c_double), intent(in), optional :: spwt_input(N)
		real(c_double), dimension(2*Ng+1,2*Ng+1), intent(out) :: pg
		call ptc2pdf(Ng,dx,N,pg,x,y,spwt_input)
	end subroutine

	subroutine c_pdfgrad(Ng,dx,pg,pg_dx,pg_dy) bind(c)
		integer(c_int), intent(in) :: Ng
		real(c_double), intent(in) :: dx, pg(2*Ng+1,2*Ng+1)
		real(c_double), dimension(2*Ng+1,2*Ng+1), intent(out) :: pg_dx,pg_dy
		pg_dx = 0.0d0
		pg_dy = 0.0d0
		pg_dx(2:2*Ng,:) = ( pg(3:2*Ng+1,:)-pg(1:2*Ng-1,:) )/2.0d0/dx
		pg_dy(:,2:2*Ng) = ( pg(:,3:2*Ng+1)-pg(:,1:2*Ng-1) )/2.0d0/dx
		call pdfgrad(Ng,dx,pg,pg_dx,pg_dy)
	end subroutine

end module