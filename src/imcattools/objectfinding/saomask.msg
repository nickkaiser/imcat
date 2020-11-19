From:	AAOEPP::HDD          "Helen Davies" 30-JAN-1996 17:32:20.56
To:	NK
CC:	HDD
Subj:	

Return-path: <squires@magicbean.berkeley.edu>
Received: from nak.Berkeley.EDU by aaoepp.aao.gov.au with PMDF#10482; Tue, 23
 Jan 1996 11:26 +1000
Received: from magicbean.Berkeley.EDU (magicbean.Berkeley.EDU [128.32.92.103])
 by nak.berkeley.edu (8.6.10/8.6.10) with SMTP id QAA11961 for
 <NK@aaoepp2.aao.GOV.AU>; Mon, 22 Jan 1996 16:25:52 -0800
Received: by magicbean.Berkeley.EDU (5.0/SMI-SVR4) id AA01594; Mon, 22 Jan 1996
 16:26:11 -0800
Date: Mon, 22 Jan 1996 16:26:11 -0800
From: squires@magicbean.berkeley.edu (Gordon Squires)
Subject: mask script (might be a duplicate message)
To: NK@aaoepp.aao.gov.au
Message-id: <9601230026.AA01594@magicbean.Berkeley.EDU>
X-Envelope-to: NK
content-length: 2003

Hi. I have sometimes found it useful to use saoimage rather than
xim to look at images. saoimage has a feature (under "cursor") where
you can click on points (like the maskmask feature in xim) and I thought
it might be useful to have a script to convert the saoimage output
into a format that imcat can use to, eg., mask objects. 

You do the following: use "cursor" in saoimage to define rectangles that you
would want for making a mask, for example. Then write these to file (click
on "region" then "write". This will make a file that looks like:

# Mon Jan 22 12:15:15 1996
# shape x, y, [x dimension, y dimension], [angle]
 POINT(173,244)
 POINT(246,307)
 POINT(262,350)
 POINT(329,410)
 POINT(113,393)
 POINT(210,435)


If you run the following script on it, it will convert these to standard
imcat format (for use, as I said, with maskobjects, etc)

-Gordon-

---------------------------------


#!/bin/sh
 
# saomask :  converts sao mask coords to imcat style
 
case $# in
0)      echo "\nsaomask :  convert saoimage coords to imcat format\n\n\
        usage: saomask inputfile > output.msk \n\n\
        where the input file is a list of points made by saoimage.\n\
        The input file has the format:\n\
            # tmp.fits\n\
            # Sun Jan 21 15:38:45 1996\n\
            # shape x, y, [x dimension, y dimension], [angle]\n\
             POINT(516,422)\n\
             POINT(586,513)\n\
                .\n\
                .\n\
\n"     1>&2; exit 2
esac
 
echo "#  imin  jmin  imax  jmax"
 
awk '{ 
        FS ="("; 
        if ($1 ~ /POINT/) 
                {print $2}  
}' $1 |
awk '{
        FS =","; 
        print $1, $2
}' |
awk '{
        FS =")"; 
        print $1, $2
}' |
awk '{  
        if( NR == 1 || NR == 2)
                printf("   %s   %s", $1, $2)
        else {
                printf("\n   %s   %s", $1, $2)
                NR=1
        }
}'  |
awk '{
        printf("   %s   %s   %s   %s\n", $1, $2, $3, $4)
}'

--QAA09205.822356658/nak.berkeley.edu--


