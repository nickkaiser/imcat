# points.tcl - command to create dots, ellipses etc on canvas


proc dot {x y radius} {
	global c _xsmall _xbig _ysmall _ybig
	if {$x > $_xsmall && $x < $_xbig && $y > $_ysmall && $y < $_ybig} {
		set X [xpix $x]
		set Y [ypix $y]
		$c create oval [expr $X - $radius] [expr $Y - $radius] [expr $X + $radius] [expr $Y + $radius] \
			-outline black
	}
}



proc ellipse {x y a b phi} {
	global c _xsmall _xbig _ysmall _ybig
	if {$x > $_xsmall && $x < $_xbig && $y > $_ysmall && $y < $_ybig} {
		set X [xpix $x]
		set Y [ypix $y]
		set cos [expr cos($phi)]
		set sin [expr sin($phi)]
		set aX [expr $a * $cos]
		set aY [expr -$a * $sin]
		set bX [expr -$b * $sin]
		set bY [expr -$b * $cos]
		$c create line [expr $X + $aX] [expr $Y + $aY] \
				[expr $X + $bX] [expr $Y + $bY] \
				[expr $X - $aX] [expr $Y - $aY] \
				[expr $X - $bX] [expr $Y - $bY] \
			 	[expr $X + $aX] [expr $Y + $aY] \
				-smooth 1 -splinesteps 10
	}
}