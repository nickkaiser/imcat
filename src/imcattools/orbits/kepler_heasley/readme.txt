Date: Tue, 11 Feb 2003 16:11:19 -1000 (HST)
From: Jim Heasley <heasley@hoku.ifa.hawaii.edu>
To: Nick Kaiser <kaiser@surf.ifa.hawaii.edu>
Subject: Celestial Mechiani

Nick,

The file appended below are the routines you asked for. They include:

        subroutine elemnts  Computes orbital elements given r and
                            rdot

        subroutine rrdot    Given orbital elements, compute r and
                            rdot

        subroutine kepler   Solves Kepler's equation

        subroutine rotate   Used by rrdot

        various vector utilities used by these codes (e.g, vdot,
            vlen, vxp)

        block data with constants

Regarding some units, r is in AU, rdot is AU/mean solar day,
emu = k*k where k = 0.01720209895 (and hence 
emu = 0.0002959122370203)


Jim Heasley

Date: Wed, 12 Feb 2003 09:01:39 -1000 (HST)
From: Jim Heasley <heasley@hoku.ifa.hawaii.edu>
To: Nick Kaiser <kaiser@ifa.hawaii.edu>
Subject: Re: Celestial Mechiani

You have it almost right:

The elements are a, e, i, Omega, omega, and T where here T is
the time of the (last?) perihelion passage. It is equivalent
to giving M at some specific time.

