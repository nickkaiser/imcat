proc line {x1 y1 x2 y2} {
	global c
	$c create line [xpix $x1] [ypix $y1] [xpix $x2] [ypix $y2]
}



proc box {} {
	global _x1 _y1 _x2 _y2
	line $_x1 $_y1 $_x2 $_y1
	line $_x2 $_y1 $_x2 $_y2
	line $_x2 $_y2 $_x1 $_y2
	line $_x1 $_y2 $_x1 $_y1
	ticks
}



proc ticks {} {
	global c _x1 _y1 _x2 _y2 _xbig _xsmall _ybig _ysmall _xrange _yrange
	set dx 1.e10
	while {$dx > $_xrange} {
		set dx [expr $dx / 10]
	}
	if {$_xrange < [expr 5 * $dx]} {set dx [expr $dx / 5]}
	set dy 1.e10
	while {$dy > $_yrange} {
		set dy [expr $dy / 10]
	}
	if {$_yrange < [expr 5 * $dy]} {set dy [expr $dy / 5]}
	set bigfracticklen 0.02
	set smallfracticklen 0.01
	set smalltickstep 0.2
	set ddx [expr $dx * $smalltickstep]
	set ddy [expr $dy * $smalltickstep]
	set bigyticklen [expr $bigfracticklen * [expr $_y2 - $_y1]]
	set bigxticklen [expr $bigfracticklen * [expr $_x2 - $_x1]]
	set smallyticklen [expr $smallfracticklen * [expr $_y2 - $_y1]]
	set smallxticklen [expr $smallfracticklen * [expr $_x2 - $_x1]]
	for {set ix [expr int(ceil($_xsmall / $dx))]} {$ix <= int(floor($_xbig / $dx))} {incr ix} {
		set x [expr $ix * $dx]
		line $x $_y1 $x [expr $_y1 + $bigyticklen]
		line $x $_y2 $x [expr $_y2 - $bigyticklen]
		$c create text [xpix $x] [ypix $_y1] -text [format "%g" $x] -anchor n
		for {set xx [expr $x + $ddx - $dx]} {$xx < ($x + $dx) && $xx < $_xbig} {set xx [expr $xx + $ddx]} {
			if {$xx > $_xsmall} {
				line $xx $_y1 $xx [expr $_y1 + $smallyticklen]
				line $xx $_y2 $xx [expr $_y2 - $smallyticklen]
			}
		}
	}
	for {set iy [expr int(ceil($_ysmall / $dy))]} {$iy <= int(floor($_ybig / $dy))} {incr iy} {
		set y [expr $iy * $dy]
		line $_x1 $y [expr $_x1 + $bigxticklen] $y
		line $_x2 $y [expr $_x2 - $bigxticklen] $y
		$c create text [xpix $_x1] [ypix $y] -text [format "%g " $y] -anchor e
		for {set yy [expr $y + $ddy - $dy]} {$yy < ($y + $dy) && $yy < $_ybig} {set yy [expr $yy + $ddy]} {
			if {$yy > $_ysmall} {
				line $_x1 $yy [expr $_x1 + $smallxticklen] $yy
				line $_x2 $yy [expr $_x2 - $smallxticklen] $yy
			}
		}
	}
}

proc xlabel {xlabeltext} {
	global c XLabelX XLabelY
	$c create text $XLabelX $XLabelY -text $xlabeltext -tags label
}



proc ylabel {ylabeltext} {
	global c YLabelX YLabelY
	$c create text $YLabelX $YLabelY -text $ylabeltext -tags label
}



proc plotlabel {plotlabeltext} {
	global c PlotLabelX PlotLabelY
	$c create text $PlotLabelX $PlotLabelY -text $plotlabeltext -tags label
}
