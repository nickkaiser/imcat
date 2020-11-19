From squires Wed Mar 29 16:43:17 1995
Received: from hare.cita.utoronto.ca by chipmunk.cita.utoronto.ca id AA08163; Wed, 29 Mar 95 16:43:15 EST
From: Gordon Squires <squires>
Received: (squires@localhost) by hare.cita (8.6.10/8.6.10) id QAA07676 for kaiser@cita; Wed, 29 Mar 1995 16:43:14 -0500
Date: Wed, 29 Mar 1995 16:43:14 -0500
Message-Id: <199503292143.QAA07676@hare.cita>
To: kaiser@cita.utoronto.ca
Subject: mma macro
Status: R

Here are the commands to make 3-d surface plots. the
read-in routine is yours! the rest is straightforward....


#
# specific application to finite simulation plots
#

readf[file_] := (
	thefile = OpenRead[file];
	Read[thefile, Character];
	N1 = Read[thefile, Number];
	N2 = Read[thefile, Number];
	f = Read[thefile, Table[Table[Number, {j,N1}], {i,N2}]];
	Close[thefile];
)

readf["!imarith -m 1000 < me.fits  | fitstoascii"]

meplot = ListPlot3D[f, PlotRange -> {{0,63},{0,63},{-100,500}}, BoxRatios -> {1,1,1}, ViewPoint -> {1.3,-2,1}, PlotLabel -> "Input"]


******Note: the following  is handy for short fits files 
	    that have "MAGIC" pixels - 
            these are reset to "Indeterminate" by mma and not plotted:

readf["!fitstoascii < SS_n2000.fits"]

ss = ListPlot3D[Map[If[#==-32768,Indeterminate,#]&,f,{2}], PlotRange -> {{0,63},{0,63},{-200,500}}, BoxRatios -> {1,1,1}, ViewPoint -> {1.3,-2,1}, PlotLabel -> "SS"]



Display["input.ps", meplot]

Display["ss.ps", ss]

