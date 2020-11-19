# plotinit.tcl
# commands to set up canvas dimension, user coords etc

proc setcanvassize {Width Height} {
	global CanvasWidth CanvasHeight
	set CanvasWidth $Width
	set CanvasHeight $Height
	catch {destroy $c}
}

proc erase {} {
	global c
	destroy $c
	createcanvas
}


proc createmenus {} {
	global mbar
	frame $mbar -relief raised -bd 2
	menubutton $mbar.file -text File -underline 0 -menu $mbar.file.menu
	menu $mbar.file.menu
	$mbar.file.menu add command -label "postscript->temp.ps" \
		-command "mkpostscript temp.ps"
	$mbar.file.menu add command -label "postscript->lpr" \
		-command "mkpostscript {|lpr -h -Pnetps1}"
	$mbar.file.menu add separator
	$mbar.file.menu add command -label "exit" -command "exit"
	pack $mbar.file -side left
}


proc createcanvas {} {
	global c CanvasWidth CanvasHeight
	canvas $c -width $CanvasWidth -height $CanvasHeight -background white
	$c bind label <1> "itemStartDrag $c %x %y"
	$c bind label <B1-Motion> "itemDrag $c %x %y"
}


proc mkpostscript {file} {
	global c
	set f [open $file w]
	puts $f [$c postscript -pagewidth 550]
}

proc setmargins {Left Bottom Right Top} {
	global BottomMargin TopMargin LeftMargin RightMargin
	set BottomMargin 	$Bottom
	set TopMargin 		$Top
	set LeftMargin 		$Left
	set RightMargin 	$Right
	setboxlimits
	setlabelpositions
}

proc setboxlimits {} {
	global BottomMargin TopMargin LeftMargin RightMargin
	global X1 Y1 X2 Y2 CanvasWidth CanvasHeight
	set X1 $LeftMargin 
	set X2 [expr $CanvasWidth - $RightMargin] 
	set Y1 [expr $CanvasHeight - $BottomMargin]
	set Y2 $TopMargin
}

proc setlabelpositions {} {
	global BottomMargin TopMargin LeftMargin RightMargin
	global X1 Y1 X2 Y2 CanvasWidth CanvasHeight
	global XLabelX XLabelY YLabelX YLabelY PlotLabelX PlotLabelY
	set XLabelX [expr ($X1 + $X2) / 2]
	set XLabelY [expr $Y1 + $BottomMargin / 2]
	set YLabelX [expr $LeftMargin / 3]
	set YLabelY [expr ($Y1 + $Y2) / 2]
	set PlotLabelX $XLabelX
	set PlotLabelY [expr $TopMargin / 2]
}


proc setuserlimits {left bottom right top} {
	global _x1 _y1 _x2 _y2 _xscale _yscale X1 Y1 X2 Y2
	global _xbig _xsmall _ybig _ysmall _xrange _yrange
	scan $left "%f" _x1
	scan $bottom "%f" _y1 
	scan $right "%f" _x2
	scan $top "%f" _y2
	set _xscale [expr [expr $X2 - $X1] / [expr $_x2 - $_x1]]
	set _yscale [expr [expr $Y2 - $Y1] / [expr $_y2 - $_y1]]
	if {$_x1 < $_x2} {
		set _xbig $_x2
		set _xsmall $_x1
	} else {
		set _xbig $_x1
		set _xsmall $_x2
	}
	if {$_y1 < $_y2} {
		set _ybig $_y2
		set _ysmall $_y1
	} else {
		set _ybig $_y1
		set _ysmall $_y2
	}
	set _xrange [expr $_xbig - $_xsmall]
	set _yrange [expr $_ybig - $_ysmall]
}

proc plotinit {} {
	global mbar c
	setcanvassize 800 800
	setmargins 100 70 30 60
	setuserlimits 0.0 0.0 1.0 1.0
	createmenus
	createcanvas
	pack $mbar $c -side top -fill x
}


proc xpix {x} {
	global X1 _x1 _xscale
	set X [expr $X1 + [expr $_xscale * [expr $x - $_x1]]]
}

proc ypix {y} {
	global Y1 _y1 _yscale
	set Y [expr $Y1 + [expr $_yscale * [expr $y - $_y1]]]
}

# Utility procedures to support dragging of items.

proc itemStartDrag {c x y} {
    global lastX lastY
    set lastX [$c canvasx $x]
    set lastY [$c canvasy $y]
}

proc itemDrag {c x y} {
    global lastX lastY
    set x [$c canvasx $x]
    set y [$c canvasy $y]
    $c move current [expr $x-$lastX] [expr $y-$lastY]
    set lastX $x
    set lastY $y
}

