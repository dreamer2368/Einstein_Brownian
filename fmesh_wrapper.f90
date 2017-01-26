module fmesh_wrapper

	use iso_c_binding, only: c_double, c_int
	use mesh

	implicit none

contains

	subroutine c_ptc2pdf(Ng,dx,dy,N,pg,x,y,spwt) bind(c)
		integer(c_int), intent(in) :: Ng, N
		real(c_double), intent(in) :: dx, dy, x(N), y(N)
		real(c_double), intent(in) :: spwt(N)
		real(c_double), dimension(0:2*Ng,0:2*Ng), intent(out) :: pg
		call ptc2pdf(Ng,dx,dy,N,pg,x,y,spwt)
	end subroutine

	subroutine c_ptc2pdf_tsc(Ng,dx,dy,N,pg,x,y,spwt) bind(c)
		integer(c_int), intent(in) :: Ng, N
		real(c_double), intent(in) :: dx, dy, x(N), y(N)
		real(c_double), intent(in) :: spwt(N)
		real(c_double), dimension(0:2*Ng,0:2*Ng), intent(out) :: pg
		call ptc2pdf_tsc(Ng,dx,dy,N,pg,x,y,spwt)
	end subroutine

	subroutine c_pdf2ptc(Ng,dx,dy,N,pg,x,y,pp) bind(c)
		integer(c_int), intent(in) :: Ng, N
		real(c_double), intent(in) :: dx, dy, x(N), y(N)
		real(c_double), dimension(0:2*Ng,0:2*Ng), intent(in) :: pg
		real(c_double), dimension(N), intent(out) :: pp
		call pdf2ptc(Ng,dx,dy,N,pg,x,y,pp)
	end subroutine

	subroutine c_pdf2ptc_tsc(Ng,dx,dy,N,pg,x,y,pp) bind(c)
		integer(c_int), intent(in) :: Ng, N
		real(c_double), intent(in) :: dx, dy, x(N), y(N)
		real(c_double), dimension(0:2*Ng,0:2*Ng), intent(in) :: pg
		real(c_double), dimension(N), intent(out) :: pp
		call pdf2ptc_tsc(Ng,dx,dy,N,pg,x,y,pp)
	end subroutine

	subroutine c_pdfgrad(Ng,dx,dy,pg,pg_dx,pg_dy) bind(c)
		integer(c_int), intent(in) :: Ng
		real(c_double), intent(in) :: dx, dy, pg(2*Ng+1,2*Ng+1)
		real(c_double), dimension(0:2*Ng,0:2*Ng), intent(out) :: pg_dx,pg_dy
		call pdfgrad(Ng,dx,dy,pg,pg_dx,pg_dy)
	end subroutine

end module