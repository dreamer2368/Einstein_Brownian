module mesh

	implicit none

!	integer, parameter :: dp = selected_real_kind(15, 307)

contains

	subroutine ptc2pdf(Ng,dx,dy,N,pg,x,y,spwt)
		integer, intent(in) :: Ng,N
		real*8, intent(in) :: dx, dy, x(N), y(N)
		real*8, intent(in) :: spwt(N)
		integer :: i,g(2,2),xg,yg
		real*8 :: f(2,2)
		real*8, dimension(0:2*Ng,0:2*Ng), intent(out) :: pg

		pg = 0.0d0
		do i=1,N
			xg = FLOOR(x(i)/dx)
			yg = FLOOR(y(i)/dy)
			g(:,1) = (/xg,xg+1/) + Ng
			g(:,2) = (/yg,yg+1/) + Ng
			if( (minval(g)<0).or.(maxval(g)>2*Ng) ) then
				cycle
			end if
			f = spwt(i)/dx/dy
			f(1,:) = f(1,:)*( 1.0d0 - (x(i)/dx-xg) )
			f(2,:) = f(2,:)*(x(i)/dx-xg)
			f(:,1) = f(:,1)*( 1.0d0 - (y(i)/dy-yg) )
			f(:,2) = f(:,2)*(y(i)/dy-yg)
			pg(g(:,1),g(:,2)) = pg(g(:,1),g(:,2)) + f
		end do
	end subroutine

	subroutine ptc2pdf_tsc(Ng,dx,dy,N,pg,x,y,spwt)
		integer, intent(in) :: Ng,N
		real*8, intent(in) :: dx, dy, x(N), y(N)
		real*8, intent(in) :: spwt(N)
		integer :: i,g(3,2),gx,gy
		real*8 :: f(3,3), fx(3),fy(3),hx,hy
		real*8, dimension(0:2*Ng,0:2*Ng), intent(out) :: pg

		pg = 0.0d0
		do i=1,N
			gx = FLOOR(x(i)/dx + 0.5d0)
			hx = x(i)/dx-gx
			fx(2) = 0.75d0 - hx*hx
			fx(1) = 0.5d0*(0.5d0-hx)*(0.5d0-hx)
			fx(3) = 0.5d0*(0.5d0+hx)*(0.5d0+hx)
			g(:,1) = (/gx-1,gx,gx+1/) + Ng
			if( gx.eq.-Ng ) then
				g(1,1) = 2*Ng
				fx(1) = 0.0d0
			elseif( gx.eq.Ng ) then
				g(3,1) = 0
				fx(3) = 0.0d0
			elseif( abs(gx)>Ng ) then
				cycle
			end if

			gy = FLOOR(y(i)/dy + 0.5d0)
			hy = y(i)/dy-gy
			fy(2) = 0.75d0 - hy*hy
			fy(1) = 0.5d0*(0.5d0-hy)*(0.5d0-hy)
			fy(3) = 0.5d0*(0.5d0+hy)*(0.5d0+hy)
			g(:,2) = (/ gy-1,gy,gy+1 /) + Ng
			if( gy.eq.-Ng ) then
				g(1,2) = 2*Ng
				fy(1) = 0.0d0
			elseif( gy.eq.Ng ) then
				g(3,2) = 0
				fy(3) = 0.0d0
			elseif( abs(gy)>Ng ) then
				cycle
			end if

			f = spwt(i)/dx/dy
			f(1,:) = f(1,:)*fx(1)
			f(2,:) = f(2,:)*fx(2)
			f(3,:) = f(3,:)*fx(3)
			f(:,1) = f(:,1)*fy(1)
			f(:,2) = f(:,2)*fy(2)
			f(:,3) = f(:,3)*fy(3)

			pg(g(:,1),g(:,2)) = pg(g(:,1),g(:,2)) + f
		end do
	end subroutine

	subroutine pdf2ptc(Ng,dx,dy,N,pg,x,y,pp)
		integer, intent(in) :: Ng, N
		real*8, intent(in) :: dx, dy, x(N), y(N)
		real*8, dimension(0:2*Ng,0:2*Ng), intent(in) :: pg
		real*8, dimension(N), intent(out) :: pp
		integer :: i,xg,yg,g(2,2)
		real*8 :: f(2,2)

		pp = 0.0d0
		do i=1,N
			xg = FLOOR(x(i)/dx)
			yg = FLOOR(y(i)/dy)
			g(:,1) = (/xg,xg+1/) + Ng
			g(:,2) = (/yg,yg+1/) + Ng
			if( (minval(g)<0).or.(maxval(g)>2*Ng) ) then
				cycle
			end if
			f = 1.0d0
			f(1,:) = f(1,:)*( 1.0d0 - (x(i)/dx-xg) )
			f(2,:) = f(2,:)*(x(i)/dx-xg)
			f(:,1) = f(:,1)*( 1.0d0 - (y(i)/dy-yg) )
			f(:,2) = f(:,2)*(y(i)/dy-yg)
			pp(i) = sum( pg(g(:,1),g(:,2))*f )
		end do
	end subroutine

	subroutine pdf2ptc_tsc(Ng,dx,dy,N,pg,x,y,pp)
		integer, intent(in) :: Ng, N
		real*8, intent(in) :: dx, dy, x(N), y(N)
		real*8, dimension(0:2*Ng,0:2*Ng), intent(in) :: pg
		real*8, dimension(N), intent(out) :: pp
		integer :: i,g(3,2),gx,gy
		real*8 :: f(3,3), fx(3),fy(3),hx,hy

		pp = 0.0d0
		do i=1,N
			gx = FLOOR(x(i)/dx + 0.5d0)
			hx = x(i)/dx-gx
			fx(2) = 0.75d0 - hx*hx
			fx(1) = 0.5d0*(0.5d0-hx)*(0.5d0-hx)
			fx(3) = 0.5d0*(0.5d0+hx)*(0.5d0+hx)
			g(:,1) = (/gx-1,gx,gx+1/) + Ng
			if( gx.eq.-Ng ) then
				g(1,1) = 2*Ng
				fx(1) = 0.0d0
			elseif( gx.eq.Ng ) then
				g(3,1) = 0
				fx(3) = 0.0d0
			elseif( abs(gx)>Ng ) then
				cycle
			end if

			gy = FLOOR(y(i)/dy + 0.5d0)
			hy = y(i)/dy-gy
			fy(2) = 0.75d0 - hy*hy
			fy(1) = 0.5d0*(0.5d0-hy)*(0.5d0-hy)
			fy(3) = 0.5d0*(0.5d0+hy)*(0.5d0+hy)
			g(:,2) = (/ gy-1,gy,gy+1 /) + Ng
			if( gy.eq.-Ng ) then
				g(1,2) = 2*Ng
				fy(1) = 0.0d0
			elseif( gy.eq.Ng ) then
				g(3,2) = 0
				fy(3) = 0.0d0
			elseif( abs(gy)>Ng ) then
				cycle
			end if

			f = 1.0d0
			f(1,:) = f(1,:)*fx(1)
			f(2,:) = f(2,:)*fx(2)
			f(3,:) = f(3,:)*fx(3)
			f(:,1) = f(:,1)*fy(1)
			f(:,2) = f(:,2)*fy(2)
			f(:,3) = f(:,3)*fy(3)
			pp(i) = sum( pg(g(:,1),g(:,2))*f )
		end do
	end subroutine


	subroutine pdfgrad(Ng,dx,dy,pg,pg_dx,pg_dy)
		integer, intent(in) :: Ng
		real*8, intent(in) :: dx,dy
		real*8, dimension(0:2*Ng,0:2*Ng), intent(in) :: pg
		real*8, dimension(0:2*Ng,0:2*Ng), intent(out) :: pg_dx,pg_dy
		pg_dx = 0.0d0
		pg_dy = 0.0d0
		pg_dx(1:2*Ng-1,:) = ( pg(2:2*Ng,:)-pg(0:2*Ng-2,:) )/2.0d0/dx
		pg_dy(:,1:2*Ng-1) = ( pg(:,2:2*Ng)-pg(:,0:2*Ng-2) )/2.0d0/dy
	end subroutine

end module