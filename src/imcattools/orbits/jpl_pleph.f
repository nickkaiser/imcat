************************************************************************
******************************* JPL CODE *******************************
************************************************************************
C from Dave Tholen 06/24/03
C------------------------------ SUBROUTINE 'PLEPH' ---------------------
C
C++++++++++++++++++++++++++
C
       SUBROUTINE  PLEPH ( JD, TARG, CENT, RRD )
C
C++++++++++++++++++++++++++
C
C     THIS SUBROUTINE READS THE JPL PLANETARY EPHEMERIS
C     AND GIVES THE POSITION AND VELOCITY OF THE POINT 'TARG'
C     WITH RESPECT TO 'CENT'.
C
C     CALLING SEQUENCE PARAMETERS:
C
C       JD = D.P. JULIAN EPHEMERIS DATE AT WHICH INTERPOLATION
C            IS WANTED.
C
C       ** NOTE THE ENTRY DPLEPH FOR A DOUBLY-DIMENSIONED TIME **
C          THE REASON FOR THIS OPTION IS DISCUSSED IN THE
C          SUBROUTINE STATE
C
C     TARG = INTEGER NUMBER OF 'TARGET' POINT.
C
C     CENT = INTEGER NUMBER OF CENTER POINT.
C
C            THE NUMBERING CONVENTION FOR 'TARG' AND 'CENT' IS:
C
C                1 = MERCURY           8 = NEPTUNE
C                2 = VENUS             9 = PLUTO
C                3 = EARTH            10 = MOON
C                4 = MARS             11 = SUN
C                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
C                6 = SATURN           13 = EARTH-MOON BARYCENTER
C                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)
C                            15 = LIBRATIONS, IF ON EPH FILE
C
C             (IF NUTATIONS ARE WANTED, SET TARG = 14. FOR LIBRATIONS,
C              SET TARG = 15. 'CENT' WILL BE IGNORED ON EITHER CALL.)
C
C      RRD = OUTPUT 6-WORD D.P. ARRAY CONTAINING POSITION AND VELOCITY
C            OF POINT 'TARG' RELATIVE TO 'CENT'. THE UNITS ARE AU AND
C            AU/DAY. FOR LIBRATIONS THE UNITS ARE RADIANS AND RADIANS
C            PER DAY. IN THE CASE OF NUTATIONS THE FIRST FOUR WORDS OF
C            RRD WILL BE SET TO NUTATIONS AND RATES, HAVING UNITS OF
C            RADIANS AND RADIANS/DAY.
C
C            NOTE: IN MANY CASES THE USER WILL NEED ONLY POSITION
C                  VALUES FOR EPHEMERIDES OR NUTATIONS. FOR
C                  POSITION-ONLY OUTPUT, THE INTEGER VARIABLE 'IPV'
C                  IN THE COMMON AREA /PLECOM/ SHOULD BE SET = 1
C                  BEFORE THE NEXT CALL TO PLEPH. (ITS DEFAULT
C                  VALUE IS 2, WHICH RETURNS BOTH POSITIONS AND
C                  RATES.)
C
C     INSIDE is .TRUE. if the input Julian Ephemeris Date (JD) is within
C            the ephemeris time span.  If not, INSIDE is set to .FALSE.
C
C     $Header: testeph.pc,v 1.1 93/08/23 14:54:07 faith Exp $
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
       SAVE

C
C     CALLING SEQUENCE DECLARATIONS
C
       DOUBLE PRECISION JD
       INTEGER          TARG
       INTEGER          CENT
       DOUBLE PRECISION RRD(6)
       LOGICAL          INSIDE
C
C     MISCELLANEOUS DECLARATIONS
C
       DOUBLE PRECISION JD2(2),JDTOT,JED(2)
       DOUBLE PRECISION PV(6,13)
       DOUBLE PRECISION EMBF(2)
       DOUBLE PRECISION VE(2)
       DOUBLE PRECISION FAC
C
       LOGICAL FIRST
       LOGICAL BSAVE
C
       INTEGER LIST(12)
       INTEGER L(2)
       INTEGER TC(2)
       INTEGER LLST(13)
       INTEGER NEMB
C
C     COMMON AREA FOR CHARACTER DATA IN RECORD 1
C
       CHARACTER*6 TTL(14,3)
       CHARACTER*6 CNAM(400)
       COMMON/CHRHDR/CNAM,TTL
C
C     COMMON AREA FOR CONSTANTS AND POINTERS IN RECORD 1
C
       DOUBLE PRECISION SS(3)
       INTEGER NCON
       DOUBLE PRECISION CVAL(400)
       DOUBLE PRECISION AU
       DOUBLE PRECISION EMRAT
       INTEGER IPT(36)
       INTEGER DENUM
       INTEGER LPT(3)
       COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT,LPT
C
C     COMMON AREA FOR 'STATE'
C
       LOGICAL KM
       LOGICAL BARY
       DOUBLE PRECISION PVSUN(6)
       COMMON/STCOMM/KM,BARY,PVSUN
C
C     COMMON AREA FOR POSITION/VELOCITY SWITCH
C
       INTEGER IPV
       COMMON/PLECOM/IPV
C
C     DATA STATEMENTS
C
       DATA JED/2*0.D0/
       DATA PV/78*0.D0/
       DATA EMBF/-1.D0,1.D0/
       DATA LIST/12*0/
       DATA L/2*0/
       DATA TC/2*0/
       DATA LLST/1,2,10,4,5,6,7,8,9,10,11,11,3/
       DATA FIRST/.TRUE./
       DATA FAC/0.D0/
       DATA NEMB/1/
C
C
C     1ST TIME IN, BE SURE EPHEMERIS IS INITIALIZED
C
       IF(FIRST) THEN
         IPV=2
         FIRST=.FALSE.
         CALL EPHOPN
         VE(1)=1.D0/(1.D0+EMRAT)
         VE(2)=EMRAT*VE(1)
       ENDIF
C
C     INITIALIZE JED FOR 'STATE' AND SET UP COMPONENT COUNT
C
       JED(1)=JD
       JED(2)=0.D0
       GO TO 11

C     ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT
C          (SEE THE DISCUSSION IN THE SUBROUTINE STATE)

       ENTRY DPLEPH(JD2,TARG,CENT,RRD)

       JED(1)=JD2(1)
       JED(2)=JD2(2)

   11  JDTOT=JED(1)+JED(2)

       IF(JDTOT .GE. SS(1) .AND. JDTOT .LE. SS(2)) GO TO 96

         INSIDE = .FALSE.
         RETURN


   96  INSIDE = .TRUE.
       NCMP   = 3*IPV
C
C     CHECK FOR NUTATION CALL
C
       IF(TARG.EQ.14) THEN
         IF(IPT(35).GT.0) THEN
           LIST(11)=IPV
           CALL STATE(JED,LIST,PV,RRD)
           LIST(11)=0
           RETURN
         ELSE
           WRITE(6,297)
   297     FORMAT(' *****  NO NUTATIONS ON THE EPHEMERIS FILE  *****')
           STOP
         ENDIF
       ENDIF
C
C     CHECK FOR LIBRATIONS
C
       IF(TARG.EQ.15) THEN
         IF(LPT(2).GT.0) THEN
           LIST(12)=IPV
           CALL STATE(JED,LIST,PV,RRD)
           LIST(12)=0
           DO 7 I=1,NCMP
     7     RRD(I)=PV(I,11)
           RETURN
         ELSE
           WRITE(6,298)
   298     FORMAT(' *****  NO LIBRATIONS ON THE EPHEMERIS FILE  *****')
           STOP
         ENDIF
       ENDIF
C
C     CHECK FOR TARGET POINT = CENTER POINT
C
       IF(TARG.EQ.CENT) THEN
         DO 1 I=1,NCMP
     1   RRD(I)=0.D0
         RETURN
       ENDIF
C
C       FORCE BARYCENTRIC OUTPUT BY 'STATE'
C
       BSAVE=BARY
       BARY=.TRUE.
C
C       SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL
C
       TC(1)=TARG
       TC(2)=CENT
       LME=0
C
       DO 2 I=1,2
       L(I)=LLST(TC(I))
       IF(L(I).LT.11) LIST(L(I))=IPV
       IF(TC(I).EQ.3) THEN
         LME=3
         FAC=-VE(1)
       ELSEIF(TC(I).EQ.10) THEN
         LME=10
         FAC=VE(2)
       ELSEIF(TC(I).EQ.13) THEN
         NEMB=I
       ENDIF
     2 CONTINUE
C
       IF(LIST(10).EQ.IPV .AND. L(1).NE.L(2)) LIST(3)=IPV-LIST(3)
C
C       MAKE CALL TO STATE
C

       CALL STATE(JED,LIST,PV,RRD)
C
C       CASE: EARTH-TO-MOON
C
       IF(TARG.EQ.10 .AND. CENT.EQ.3) THEN
         DO 3 I=1,NCMP
     3   RRD(I)=PV(I,10)
C
C       CASE: MOON-TO-EARTH
C
       ELSEIF(TARG.EQ.3 .AND. CENT.EQ.10) THEN
         DO 4 I=1,NCMP
     4   RRD(I)=-PV(I,10)
C
C       CASE: EMBARY TO MOON OR EARTH
C
       ELSEIF((TARG.EQ.13 .OR. CENT.EQ.13) .AND. LIST(10).EQ.IPV) THEN
         DO 5 I=1,NCMP
     5   RRD(I)=PV(I,10)*FAC*EMBF(NEMB)
C
C       OTHERWISE, GET EARTH OR MOON VECTOR AND THEN GET OUTPUT VECTOR
C
       ELSE
         DO 6 I=1,NCMP
         PV(I,11)=PVSUN(I)
         PV(I,13)=PV(I,3)
         IF(LME.GT.0) PV(I,LME)=PV(I,3)+FAC*PV(I,10)
     6   RRD(I)=PV(I,TARG)-PV(I,CENT)
C
       ENDIF
C
C     CLEAR 'STATE' BODY ARRAY AND RESTORE BARYCENTER FLAG
C
       LIST(3)=0
       LIST(L(1))=0
       LIST(L(2))=0
       BARY=BSAVE
C
C     THAT'S ALL
C
       RETURN
       END

C---------------SUBROUTINE 'EPHOPN'  : VERSION #3--------------------
C
C     $Header: testeph.pc,v 1.1 93/08/23 14:54:07 faith Exp $
C
C++++++++++++++++++++++++
C
       SUBROUTINE EPHOPN
C
C++++++++++++++++++++++++
C
C     THIS SUBROUTINE OPENS THE JPL PLANETARY EPHEMERIS 12,
C     EXPECTED TO BE 'JPLEPH'.
C
       SAVE
C
C       COMMON AREA FOR CHARACTER DATA IN RECORD 1
C
       CHARACTER*6 TTL(14,3)
       CHARACTER*6 CNAM(400)
       COMMON/CHRHDR/CNAM,TTL
C
C       COMMON AREA FOR CONSTANTS AND POINTERS IN RECORD 1
C
       DOUBLE PRECISION SS(3)
       INTEGER NCON
       DOUBLE PRECISION CVAL(400)
       DOUBLE PRECISION AU
       DOUBLE PRECISION EMRAT
       INTEGER IPT(3,12)
       INTEGER DENUM
       INTEGER LPT(3)
       COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT,LPT
C
C
C       COMMON AREA FOR EPH DATA BUFFERS
C
       PARAMETER (IBDIM=2036)
       DOUBLE PRECISION BUF(IBDIM/2)
       INTEGER KSIZE
       COMMON/EPIB/BUF
       COMMON/EPKS/KSIZE
C
       INTEGER IRECSZ
C
       LOGICAL FIRST
C
C
       DATA FIRST/.TRUE./
C
C  *****************************************************************
C
C    **************  NOTE : POSSIBLE USER CHANGE  ************
C
C       THE NEXT STATEMENT DEFINES THE LENGTH OF THE RECORDS
C       IN THE DIRECT-ACCESS FILE. IF THE 'RECL' SPECIFICATION IN
C       THE FORTRAN COMPILER ON YOUR MACHINE EXPECTS THE LENGTH TO
C       BE IN WORDS, LEAVE THE STATEMENT AS IT IS. IF THE LENGTH
C       IS EXPECTED IN BYTES, CHANGE THE '1' TO A '4'.
C
       NRECL=4
C
C  *****************************************************************

       KSIZE = 2036

       IF(FIRST) THEN

         IF ( KSIZE .EQ. 0 ) THEN
            WRITE (*,*)  ' EPHOPN: KSIZE was not set. ' //
      .                  ' Process stopped.'
            STOP
         ENDIF

         IF(KSIZE .GT. IBDIM) THEN
           WRITE(6,266)IBDIM,KSIZE
  266      FORMAT(/'***** THE PARAMETER IBDIM IS SET AT ONLY',I6,
      .       '.  IT MUST BE INCREASED TO AT LEAST',I6,' *****')
           STOP
         ENDIF

         IRECSZ=NRECL*KSIZE


         OPEN(12,
      *       FILE='E:\DE405\DE405',
      *       ACCESS='DIRECT',
      *       FORM='UNFORMATTED',
      *       RECL=IRECSZ,
      *       STATUS='OLD')

         READ(12,REC=1)TTL,CNAM,SS,NCON,AU,EMRAT,IPT,DENUM,LPT
         READ(12,REC=2)CVAL

C      WRITE(*,'(14A6)')(TTL(I,1),I=1,14)

         FIRST=.FALSE.

       ENDIF
C
       RETURN
C
       END

C--------------------------- SUBROUTINE 'STATE' ------------------------
C
C     $Header: testeph.pc,v 1.1 93/08/23 14:54:07 faith Exp $
C
C++++++++++++++++++++++++++++++++
C
       SUBROUTINE STATE(JED,LIST,PV,NUT)
C
C++++++++++++++++++++++++++++++++
C
C     THIS SUBROUTINE READS AND INTERPOLATES THE JPL PLANETARY EPHEMERIS
C     FILE.
C
C     CALLING SEQUENCE PARAMETERS:
C
C     INPUT:
C
C         JED   DP 2-WORD JULIAN EPHEMERIS EPOCH AT WHICH INTERPOLATION
C               IS WANTED.  ANY COMBINATION OF JED(1)+JED(2) WHICH FALLS
C               WITHIN THE TIME SPAN ON THE FILE IS A PERMISSIBLE EPOCH.
C
C                A. FOR EASE IN PROGRAMMING, THE USER MAY PUT THE
C                   ENTIRE EPOCH IN JED(1) AND SET JED(2)=0.
C
C                B. FOR MAXIMUM INTERPOLATION ACCURACY, SET JED(1) =
C                   THE MOST RECENT MIDNIGHT AT OR BEFORE INTERPOLATION
C                   EPOCH AND SET JED(2) = FRACTIONAL PART OF A DAY
C                   ELAPSED BETWEEN JED(1) AND EPOCH.
C
C                C. AS AN ALTERNATIVE, IT MAY PROVE CONVENIENT TO SET
C                   JED(1) = SOME FIXED EPOCH, SUCH AS START OF INTEGRATION,
C                   AND JED(2) = ELAPSED INTERVAL BETWEEN THEN AND EPOCH.
C
C        LIST   12-WORD INTEGER ARRAY SPECIFYING WHAT INTERPOLATION
C               IS WANTED FOR EACH OF THE BODIES ON THE FILE.
C
C                         LIST(I)=0, NO INTERPOLATION FOR BODY I
C                                =1, POSITION ONLY
C                                =2, POSITION AND VELOCITY
C
C               THE DESIGNATION OF THE ASTRONOMICAL BODIES BY I IS:
C
C                         I = 1: MERCURY
C                           = 2: VENUS
C                           = 3: EARTH-MOON BARYCENTER
C                           = 4: MARS
C                           = 5: JUPITER
C                           = 6: SATURN
C                           = 7: URANUS
C                           = 8: NEPTUNE
C                           = 9: PLUTO
C                           =10: GEOCENTRIC MOON
C                           =11: NUTATIONS IN LONGITUDE AND OBLIQUITY
C                           =12: LUNAR LIBRATIONS (IF ON FILE)
C
C
C     OUTPUT:
C
C          PV   DP 6 X 11 ARRAY THAT WILL CONTAIN REQUESTED INTERPOLATED
C               QUANTITIES.  THE BODY SPECIFIED BY LIST(I) WILL HAVE ITS
C               STATE IN THE ARRAY STARTING AT PV(1,I).  (ON ANY GIVEN
C               CALL, ONLY THOSE WORDS IN 'PV' WHICH ARE AFFECTED BY THE
C               FIRST 10 'LIST' ENTRIES (AND BY LIST(12) IF LIBRATIONS ARE
C               ON THE FILE) ARE SET.  THE REST OF THE 'PV' ARRAY
C               IS UNTOUCHED.)  THE ORDER OF COMPONENTS STARTING IN
C               PV(1,I) IS: X,Y,Z,DX,DY,DZ.
C
C               ALL OUTPUT VECTORS ARE REFERENCED TO THE EARTH MEAN
C               EQUATOR AND EQUINOX OF EPOCH. THE MOON STATE IS ALWAYS
C               GEOCENTRIC; THE OTHER NINE STATES ARE EITHER HELIOCENTRIC
C               OR SOLAR-SYSTEM BARYCENTRIC, DEPENDING ON THE SETTING OF
C               COMMON FLAGS (SEE BELOW).
C
C               LUNAR LIBRATIONS, IF ON 12, ARE PUT INTO PV(K,11) IF
C               LIST(12) IS 1 OR 2.
C
C         NUT   DP 4-WORD ARRAY THAT WILL CONTAIN NUTATIONS AND RATES,
C               DEPENDING ON THE SETTING OF LIST(11).  THE ORDER OF
C               QUANTITIES IN NUT IS:
C
C                        D PSI  (NUTATION IN LONGITUDE)
C                        D EPSILON (NUTATION IN OBLIQUITY)
C                        D PSI DOT
C                        D EPSILON DOT
C
C     PARAMETERS IN COMMON:
C
C     COMMON AREA EPUNIT:
C
C        FILE   INTEGER OF THE UNIT CONTAINING THE EPHEMERIS. DEFAULT = 12
C
C     COMMON AREA STCOMM:
C
C          KM   LOGICAL FLAG DEFINING PHYSICAL UNITS OF THE OUTPUT
C               STATES. KM = .TRUE., KM AND KM/SEC
C                          = .FALSE., AU AND AU/DAY
C               DEFAULT VALUE = .FALSE.  (KM DETERMINES TIME UNIT
C               FOR NUTATIONS AND LIBRATIONS.  ANGLE UNIT IS ALWAYS RADIANS.)
C
C        BARY   LOGICAL FLAG DEFINING OUTPUT CENTER.
C               ONLY THE 9 PLANETS ARE AFFECTED.
C                        BARY = .TRUE. =\ CENTER IS SOLAR-SYSTEM BARYCENTER
C                             = .FALSE. =\ CENTER IS SUN
C               DEFAULT VALUE = .FALSE.
C
C       PVSUN   DP 6-WORD ARRAY CONTAINING THE BARYCENTRIC POSITION AND
C               VELOCITY OF THE SUN.
C
C
       SAVE
C
       DOUBLE PRECISION JED(2)
       INTEGER LIST(12)
       DOUBLE PRECISION PV(6,11)
       DOUBLE PRECISION NUT(4)
C
       DOUBLE PRECISION T(2)
       DOUBLE PRECISION AUFAC
       DOUBLE PRECISION JD(4)
       DOUBLE PRECISION S
C
       PARAMETER (IBDIM=2036)
       DOUBLE PRECISION BUF(IBDIM/2)
       INTEGER KSIZE
       COMMON/EPIB/BUF
       COMMON/EPKS/KSIZE
C
C
       LOGICAL FIRST
C
C     COMMON AREA FOR CHARACTER DATA IN RECORD 1
C
       CHARACTER*6 TTL(14,3)
       CHARACTER*6 CNAM(400)
       COMMON/CHRHDR/CNAM,TTL
C
C     COMMON AREA FOR CONSTANTS AND POINTERS IN RECORD 1
C
       DOUBLE PRECISION SS(3)
       INTEGER NCON
       DOUBLE PRECISION CVAL(400)
       DOUBLE PRECISION AU
       DOUBLE PRECISION EMRAT
       INTEGER L(3,12)
       INTEGER DENUM
       INTEGER LPT(3)
       COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,L,LPT
C
       LOGICAL KM
       LOGICAL BARY
       DOUBLE PRECISION PVSUN(3,2)
       COMMON/STCOMM/KM,BARY,PVSUN
C
       INTEGER NRL
C
C     DATA STATEMENTS
C
       DATA AUFAC/1.D0/
       DATA FIRST/.TRUE./
       DATA NRL/0/
C
C
C     1ST TIME IN, GET POINTER DATA, ETC., FROM EPH FILE
C
       IF(FIRST) THEN
         FIRST=.FALSE.
         CALL EPHOPN
         IF(KM) THEN
           T(2)=SS(3)*86400.D0
         ELSE
           T(2)=SS(3)
           AUFAC=1.D0/AU
         ENDIF
       ENDIF


C
C     MAIN ENTRY POINT -- CHECK EPOCH AND READ RIGHT RECORD
C
       S=JED(1)-.5D0
       CALL SPLIT(S,JD(1))
       CALL SPLIT(JED(2),JD(3))
       JD(1)=JD(1)+JD(3)+.5D0
       JD(2)=JD(2)+JD(4)
       CALL SPLIT(JD(2),JD(3))
       JD(1)=JD(1)+JD(3)
C
C     ERROR RETURN OF EPOCH OUT OF RANGE
C
       IF(JD(1).LT.SS(1) .OR. JD(1)+JD(4).GT.SS(2)) THEN
          WRITE (*,*) ' STATE: Epoch out of range.'
          STOP
       ENDIF
C
C     CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL
C
       NR=IDINT((JD(1)-SS(1))/SS(3))+3
       IF(JD(1).EQ.SS(2)) NR=NR-1
       T(1)=((JD(1)-(DBLE(NR-3)*SS(3)+SS(1)))+JD(4))/SS(3)
C
C     READ CORRECT RECORD IF NOT IN CORE
C
       IF(NR.NE.NRL) THEN

         NRL=NR
         READ(12,REC=NR,ERR=99)(BUF(K),K=1,KSIZE/2)

       ENDIF
C
C     INTERPOLATE SSBARY SUN
C
       CALL INTERP(BUF(L(1,11)),T,L(2,11),3,L(3,11),2,PVSUN)
       DO 6 I=1,6
6     PVSUN(I,1)=PVSUN(I,1)*AUFAC

C
C     CHECK AND INTERPOLATE WHICHEVER BODIES ARE REQUESTED
C

       DO 3 I=1,10

        IF(LIST(I).LE.0) GO TO 5

        IF(L(2,I).LE.0)CALL ERRPRT(I,'th body requested - not on file')

         CALL INTERP(BUF(L(1,I)),T,L(2,I),3,L(3,I),LIST(I),PV(1,I))

         DO J=1,LIST(I)*3
         IF(I.LE.9 .AND. .NOT.BARY) THEN
           PV(J,I)=PV(J,I)*AUFAC-PVSUN(J,1)
         ELSE
           PV(J,I)=PV(J,I)*AUFAC
         ENDIF
         ENDDO

     5 CONTINUE
     3 CONTINUE
C
C       DO NUTATIONS IF REQUESTED (AND IF ON FILE)
C
       IF(LIST(11).GT.0 .AND. L(2,12).GT.0)
      * CALL INTERP(BUF(L(1,12)),T,L(2,12),2,L(3,12),LIST(11),NUT)
C
C       GET LIBRATIONS IF REQUESTED (AND IF ON FILE)
C
       IF(LPT(2).GT.0 .AND. LIST(12).GT.0)
      * CALL INTERP(BUF(LPT(1)),T,LPT(2),3,LPT(3),LIST(12),PV(1,11))
C
C       THAT'S ALL
C
       RETURN
C
C       ERROR EXIT
C

    99 WRITE(6,299)
   299 FORMAT(' *****  ERROR RETURN IN CALL TO STATE FROM PLEPH  *****')
       RETURN
C
       END

C------------------------------ SUBROUTINE 'CONST' ---------------------
C
C+++++++++++++++++++++++++++++
C
       SUBROUTINE CONST(NAM,VAL,SSS,N)
C
C+++++++++++++++++++++++++++++
C
C     THIS SUBROUTINE OBTAINS THE CONSTANTS FROM THE EPHEMERIS FILE
C
C     CALLING SEQEUNCE PARAMETERS (ALL OUTPUT):
C
C       NAM = CHARACTER*6 ARRAY OF CONSTANT NAMES
C
C       VAL = D.P. ARRAY OF VALUES OF CONSTANTS
C
C       SSS = D.P. JD START, JD STOP, STEP OF EPHEMERIS
C
C         N = INTEGER NUMBER OF ENTRIES IN 'NAM' AND 'VAL' ARRAYS
C
C     THE ARRAYS 'NAM' AND 'VAL' MUST HAVE SUFFICIENTLY LARGE
C     DIMESNIONS TO ACCOMMODATE ALL ENTRIES ON THE FILE. THE
C     VALUE 400 IS A SAFE UPPER LIMIT FOR THIS DIMENSION.
C
       SAVE
C
       CHARACTER*6 NAM(*)
       DOUBLE PRECISION VAL(*)
       DOUBLE PRECISION SSS(3)
       INTEGER N
C
C       COMMON AREA FOR CHARACTER DATA IN RECORD 1
C
       CHARACTER*6 TTL(14,3)
       CHARACTER*6 CNAM(400)
       COMMON/CHRHDR/CNAM,TTL
C
C       COMMON AREA FOR CONSTANTS AND POINTERS IN RECORD 1
C
       DOUBLE PRECISION SS(3)
       INTEGER NCON
       DOUBLE PRECISION CVAL(400)
       DOUBLE PRECISION AU
       DOUBLE PRECISION EMRAT
       INTEGER IPT(36)
       INTEGER DENUM
       INTEGER LPT(3)
       COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT,LPT
C
C
C       OPEN THE FILE AND READ RECORD 1
C
       CALL EPHOPN
C
C       COPY CORRECT DATA TO CALLING SEQUENCE VARIABLES
C
       N=NCON
       DO 1 I=1,3
     1 SSS(I)=SS(I)
       DO 2 I=1,NCON
       NAM(I)=CNAM(I)
     2 VAL(I)=CVAL(I)
C
       RETURN
C
       END

C------------------------------- SUBROUTINE 'INTERP' -------------------
C
C     $Header: testeph.pc,v 1.1 93/08/23 14:54:07 faith Exp $
C
C+++++++++++++++++++++++++++++++++
C
       SUBROUTINE INTERP (BUF, T, NCF, NCM, NA, FL, PV)
C
C+++++++++++++++++++++++++++++++++
C
C     This subroutine differentiates and interpolates a set of
C     Chebyshev coefficients to give position and velocity.
C
C     Calling sequence parameters:
C
C       Input:
C
C         BUF   1st location of array of D.P. Chebyshev coefficients of position
C
C           T   T(1) IS DP FRACTIONAL TIME IN INTERVAL COVERED BY
C               COEFFICIENTS AT WHICH INTERPOLATION IS WANTED
C               (0 .LE. T(1) .LE. 1).  T(2) IS DP LENGTH OF WHOLE
C               INTERVAL IN INPUT TIME UNITS.
C
C         NCF   # OF COEFFICIENTS PER COMPONENT
C
C         NCM   # OF COMPONENTS PER SET OF COEFFICIENTS
C
C          NA   # OF SETS OF COEFFICIENTS IN FULL ARRAY
C               (I.E., # OF SUB-INTERVALS IN FULL INTERVAL)
C
C          FL   INTEGER FLAG: =1 FOR POSITIONS ONLY
C                             =2 FOR POS AND VEL
C
C
C       OUTPUT:
C
C         PV   INTERPOLATED QUANTITIES REQUESTED.  DIMENSION
C               EXPECTED IS PV(NCM,FL), DP.
C
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
       SAVE
C
C     DECLARATIONS OF SUBROUTINE ARGUMENTS
C
       DOUBLE PRECISION BUF(NCF,NCM,*)
       DOUBLE PRECISION T(2)
       INTEGER NCF
       INTEGER NCM
       INTEGER NA
       INTEGER FL
       DOUBLE PRECISION PV(NCM,*)
C
C
       INTEGER NP
       INTEGER NV
C
       DOUBLE PRECISION TWOT
       DOUBLE PRECISION PC(18)
       DOUBLE PRECISION VC(18)
C
       DATA NP/2/
       DATA NV/3/
       DATA TWOT/0.D0/
       DATA PC(1),PC(2)/1.D0,0.D0/
       DATA VC(2)/1.D0/
C
C       ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET
C       OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME
C       WITHIN THAT SUBINTERVAL.
C
       DNA=DBLE(NA)
       DT1=DINT(T(1))
       TEMP=DNA*T(1)
       L=IDINT(TEMP-DT1)+1
C
C         TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1)
C
       TC=2.D0*(DMOD(TEMP,1.D0)+DT1)-1.D0
C
C       CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,
C       AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.
C       (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE
C       CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)
C
       IF(TC.NE.PC(2)) THEN
         NP=2
         NV=3
         PC(2)=TC
         TWOT=TC+TC
       ENDIF
C
C       BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED
C       AND ARE STORED IN THE ARRAY 'PC'.
C
       IF(NP.LT.NCF) THEN
         DO 1 I=NP+1,NCF
         PC(I)=TWOT*PC(I-1)-PC(I-2)
     1   CONTINUE
         NP=NCF
       ENDIF
C
C       INTERPOLATE TO GET POSITION FOR EACH COMPONENT
C
       DO 2 I=1,NCM
       PV(I,1)=0.D0
       DO 3 J=NCF,1,-1
       PV(I,1)=PV(I,1)+PC(J)*BUF(J,I,L)
     3 CONTINUE
     2 CONTINUE
       IF(FL.LE.1) RETURN
C
C       IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH
C       DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.
C
       VFAC=(DNA+DNA)/T(2)
       VC(3)=TWOT+TWOT
       IF(NV.LT.NCF) THEN
         DO 4 I=NV+1,NCF
         VC(I)=TWOT*VC(I-1)+PC(I-1)+PC(I-1)-VC(I-2)
     4   CONTINUE
         NV=NCF
       ENDIF
C
C       INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT
C
       DO 5 I=1,NCM
       PV(I,2)=0.D0
       DO 6 J=NCF,2,-1
       PV(I,2)=PV(I,2)+VC(J)*BUF(J,I,L)
     6 CONTINUE
       PV(I,2)=PV(I,2)*VFAC
     5 CONTINUE
C
       RETURN
C
       END

C------------------------------- SUBROUTINE 'SPLIT' --------------------
C
C     $Header: testeph.pc,v 1.1 93/08/23 14:54:07 faith Exp $
C
C+++++++++++++++++++++++++
C
       SUBROUTINE SPLIT ( TT, FR )
C
C+++++++++++++++++++++++++
C
C     This subroutine breaks a D.P. number into a D.P. integer
c     and ha D.P. fractional part.
c
C     Calling sequence parameters:
C
C       TT = D.P. input number
C
C       FR = D.P. 2-word output array.
C            FR(1) contains integer part
C            FR(2) contains fractional part
C
C            For negative input numbers, FR(1) contains the next
C            more negative integer; FR(2) contains a positive fraction.
C
       DOUBLE PRECISION  TT
       DOUBLE PRECISION  FR(2)


       FR(1) = DINT(TT)
       FR(2) = TT - FR(1)

       IF ( (TT .LT. 0.D0) .AND. (FR(2) .NE. 0.D0) ) THEN

C
C        Make adjustments for negative input number.
C
          FR(1) = FR(1) - 1.D0
          FR(2) = FR(2) + 1.D0

       ENDIF

       RETURN
       END


C      $Header: testeph.pc,v 1.1 93/08/23 14:54:07 faith Exp $
C
        SUBROUTINE  ERRPRT (I, MSG)

        CHARACTER*(*)  MSG
        INTEGER        I

        WRITE (*,200)  I, MSG
200    FORMAT('ERROR #',I8,2X,A50)

        STOP ' ERROR '
        END

       subroutine selcon(nams,nns,vals)

c  Input

c           nams :  a list of names (character*6)

c           nns  :  number of names in the list

c  Output

c           vals :  values corresponding to the input names



       double precision vlc(400),sss(3),vals(1)

       character*6 nmc(400),nams(1)

       call const(nmc,vlc,sss,nnc)

       do 2 i=1,nns

       do j=1,nnc
       if(nams(i) .eq. nmc(j)) go to 1
       enddo

       vals(i)=0.0
       go to 2

    1  vals(i)=vlc(j)

    2  continue

       return
       end
