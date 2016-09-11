module mesh

	implicit none

!	integer, parameter :: dp = selected_real_kind(15, 307)

contains

	subroutine ptc2pdf(Ng,dx,N,pg,x,y,spwt_input)
		integer, intent(in) :: Ng,N
		real*8, intent(in) :: dx, x(N), y(N)
		real*8, intent(in), optional :: spwt_input(N)
		integer :: i,g(2,2),xg,yg
		real*8 :: f(2,2)
		real*8 :: spwt(N)
		real*8, dimension(2*Ng+1,2*Ng+1), intent(out) :: pg
		if( .not. present(spwt_input) ) then
			spwt = 1.0d0/N
		else
			spwt = spwt_input
		end if

		pg = 0.0d0
		do i=1,N
			xg = FLOOR(x(i)/dx)
			yg = FLOOR(y(i)/dx)
			g(:,1) = (/xg,xg+1/) + Ng
			g(:,2) = (/yg,yg+1/) + Ng
			f = spwt(i)
			f(1,:) = f(1,:)*(1.0d0-x(i)+xg*dx)
			f(2,:) = f(2,:)*(x(i)-xg*dx)
			f(:,1) = f(:,1)*(1.0d0-y(i)+yg*dx)
			f(:,2) = f(:,2)*(y(i)-yg*dx)
			pg(g(:,1),g(:,2)) = pg(g(:,1),g(:,2)) + f
		end do
		pg = pg/dx
	end subroutine

	subroutine pdfgrad(Ng,dx,pg,pg_dx,pg_dy)
		integer, intent(in) :: Ng
		real*8, intent(in) :: dx, pg(2*Ng+1,2*Ng+1)
		real*8, dimension(2*Ng+1,2*Ng+1), intent(out) :: pg_dx,pg_dy
		pg_dx = 0.0d0
		pg_dy = 0.0d0
		pg_dx(2:2*Ng,:) = ( pg(3:2*Ng+1,:)-pg(1:2*Ng-1,:) )/2.0d0/dx
		pg_dy(:,2:2*Ng) = ( pg(:,3:2*Ng+1)-pg(:,1:2*Ng-1) )/2.0d0/dx
	end subroutine

end module