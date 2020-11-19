c
c	---------------------------------------------------
c	NK: changed to a integer function
c	subroutine elmnts(elem,emu,t,r,rd)
	function ilmnts(elem,emu,t,r,rd)
	implicit real*8 (a-h, o-z)
c
c	given r and rdot vectors determine the orbital elements
c
	dimension elem(6),r(3),rd(3),h(3)
	common /consts/ pi, twopi, xmu, dummy(4)
c
c	compute components of h
c
c	h(1) = r(2)*rd(3) - r(3)*rd(2)
c	h(2) = r(3)*rd(1) - r(1)*rd(3)
c	h(3) = r(1)*rd(2) - r(2)*rd(1)
c	hm = dsqrt(vdot(h,h))
	call vxp(r,rd,h)
	hm = vlen(h)
	elem(3) = dacos(h(3)/hm)
	hsi = hm*dsin(elem(3))
	c = -h(2)/hsi
	s = h(1)/hsi
	elem(4) = datan2(s,c)
	if(elem(4).lt.0.0) elem(4) = elem(4) + twopi
c	rm = dsqrt(vdot(r,r))
c	rdm = dsqrt(vdot(rd,rd))
	rm = vlen(r)
	rdm = vlen(rd)
	elem(1) = emu*rm/(2.*emu - rm*rdm*rdm)
	if(elem(1).lt.0) then
c		write(*,2)
c2		format(' orbit is hyperboic')
c		stop
		ilmnts = 1
		return
	endif
	elem(2) = dsqrt(1.0 - hm*hm/(emu*elem(1)))
	csu = (c*r(1) + s*r(2))/rm
	u = dacos(csu)
	if(r(3).lt.0.0) u = twopi - u
	cv = (hm*hm/(emu*rm) - 1.0)/elem(2)
	sv = hm*vdot(r,rd)/(elem(2)*emu*rm)
	v = datan2(sv,cv)
	if(v.lt.0) v = v + twopi
	elem(5) = u - v
	if(elem(5).lt.0.) elem(5) = elem(5) + twopi
	cea = rm*cv/elem(1) + elem(2)
	sea = rm*sv/(elem(1)*dsqrt(1.-elem(2)**2))
	ea = datan2(sea,cea)
c	if(ea.lt.0.) ea = ea + twopi
	per = elem(1)*dsqrt(elem(1)/emu)
	elem(6) = t - (ea-elem(2)*sea)*per
	per = twopi*per
3	continue
	if(elem(6).lt.0) then
		elem(6) = elem(6) + per
		go to 3
	endif
	ilmnts = 0
	return
	end

c
c	---------------------------------------------------
c	NK: changed to a function
c	subroutine rrdot(elem,emu,t,r,rd)
	function irrdot(elem,emu,t,r,rd)
	implicit real*8 (a-h,o-z)
c
c       given r and rdot compute orbital elements
c

c
c	calculation of eccentric anamoly from mean motion 
c	en = mean motion
c
	dimension elem(6),r(3),rd(3),rz(3),rzd(3)
	en = dsqrt(emu/elem(1))/elem(1)
	em = en*(t - elem(6))
c	NK: changed this as kepler is now a function returning int
	kep = kepler(em,elem(2),ea)
	if (kep.ne.0) then
		irrdot = 1
		return
	endif
c
c	calculate rectangular coordinates in orbital plane
c	
	srt = dsqrt(1.0 - elem(2)**2)
	edo = en/(1.0 - elem(2)*dcos(ea))
	rz(1) = elem(1)*(dcos(ea) - elem(2))
	rz(2) = elem(1)*srt*dsin(ea)
	rz(3) = 0.0
	rzd(1) = -elem(1)*edo*dsin(ea)
	rzd(2) = elem(1)*srt*edo*dcos(ea)
	rzd(3) = 0.0
c
c	now rotate to the ecliptic coordiante system
c
	call rotate(r,-elem(5),3,rz)
	call rotate(rd,-elem(5),3,rzd)
	call rotate(rz,-elem(3),1,r)
	call rotate(rzd,-elem(3),1,rd)
	call rotate(r,-elem(4),3,rz)
	call rotate(rd,-elem(4),3,rzd)
	irrdot = 0
	return
	end
c
c
c	---------------------------------------------------
c
c	NK: this was writing to stdout and messing up lc
c	so I converted it to a function so that non-zero
c	result flags an error
c	subroutine kepler(em,e,ea)
	function kepler(em,e,ea)
c
	implicit real*8 (a-h,o-z)
c
c	iterative solution of kepler's equation
c
	common /consts/ pi, twopi, xmu, dummy(4)
	emr = dmod(em,twopi)
	ea = emr
	do 1 i = 1, 99
		del = (ea - e*dsin(ea) - emr)/(1.0 - e*dcos(ea))
		ea = ea - del
c		if(dabs(del).le.(1.d-7*(1. + dabs(ea)))) return
		if(dabs(del).le.(1.d-7*(1. + dabs(ea)))) then
			kepler = 0
			return
		endif
1	continue
c	write(*,3) em, e
c	3	format(' Kepler failed for ',2f15.5)
c	2	return
	kepler = 1
	return
	end
c
c	---------------------------------------------------
c
	subroutine rotate(r,a,k,rz)
	implicit real*8 (a-h, o-z)
c
c	rotate vector rz by angle a about the k axis 
c	k = 1, 2, 3 for x, y, z respectively
c
	dimension r(3), rz(3)
	data aold, c, s / 0.0, 1.0, 0.0/
	if(a.ne.aold) then
		c = dcos(a)
		s= dsin(a)
		aold = a
	endif
	i = 1 + mod(k,3)
	j = 1 + mod(i,3)
	r(k) = rz(k)
	r(i) = c*rz(i) + s*rz(j)
	r(j) = - s*rz(i) + c*rz(j)
	return
	end
c
c	---------------------------------------------------
c
	subroutine vcopy(v1,v2)
c
c	set v2 = v1
c
	implicit real*8 (a-h,o-z)
	dimension v1(3), v2(3)
	v2(1) = v1(1)
	v2(2) = v1(2)
	v2(3) = v1(3)
	return
	end
c
c	---------------------------------------------------
c
	subroutine vsum(x,y,z)
c
c	z = x + y
c
	implicit real*8 (a-h, o-z)
	dimension x(3), y(3), z(3)
	do 1 i = 1, 3
		z(i) = x(i) + y(i)
1	continue
	return
	end
c
c	---------------------------------------------------
c
	subroutine vsub(x,y,z)
c
c	z = x - y
c
	implicit real*8 (a-h, o-z)
	dimension x(3), y(3), z(3)
	do 1 i = 1, 3
		z(i) = x(i) - y(i)
1	continue
	return
	end
c
c	---------------------------------------------------
c
	subroutine vzero(x)
c
c	x = 0
c
	implicit real*8 (a-h, o-z)
	dimension x(3)
	do 1 i = 1, 3
		x(i) = 0.0
1	continue
	return
	end
c
c	---------------------------------------------------
c
	function vdot(x,y)
	implicit real*8 (a-h,o-z)
c
c	computer the dot project of 2 vectors x and y
c
	dimension x(3),y(3)
	vdot = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
	return
	end
c
c	---------------------------------------------------
c
	subroutine vxp(x,y,z)
	implicit real*8 (a-h,o-z)
c
c	compute the vector cross product z = x cross y
c
	dimension x(3), y(3), z(3)
	z(1) = (x(2)*y(3) - x(3)*y(2)) 
        z(2) = (y(1)*x(3) - y(3)*x(1))  
        z(3) = (x(1)*y(2) - x(2)*y(1))
	return
	end
c
c	---------------------------------------------------
c
	function vlen(x)
c
c	return the length of vector x
c
	implicit real*8 (a-h,o-z)
	dimension x(3)
	vlen = dsqrt(vdot(x,x))
	return
	end
c
c	---------------------------------------------------
c
	function stp(x,y,z)
	implicit real*8 (a-h,o-z)
c
c	compute the scaler triple product of x dotted into
c	y cross z
c
	dimension x(3), y(3), z(3)
	stp = x(1)*(y(2)*z(3) - y(3)*z(2)) +
     &        x(2)*(z(1)*y(3) - z(3)*y(1)) + 
     &        x(3)*(y(1)*z(2) - y(2)*z(1))
	return
	end
c
c	---------------------------------------------------
c
        block data
	implicit real*8 (a-h,o-z)
	common /consts/ pi, twopi, xmu, xk, eps, c, rearth
	data pi / 3.1415926535897931/
        data twopi / 6.2831853071795862 /
	data xk / 0.01720209895 /
	data xmu / 0.0002959122370203 /
        data eps / 0.40927971 /
	data c / 173.1428 /
	data rearth / 4.254930d-5/
	end

