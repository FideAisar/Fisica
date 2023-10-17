      PROGRAM LHAFileConvert
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      
********LESHOUCHES USER PROCESS EVENT COMMON BLOCK *******************
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500) 
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      REAL*8 XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP  
      DIMENSION IDUP(MAXNUP),
     &   ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &   VTIMUP(MAXNUP),SPINUP(MAXNUP) 

C...Lines to read in assumed never longer than 200 characters. 
      PARAMETER (MAXLEN=200)
      CHARACTER*(MAXLEN) STRING,INPUTFILE,OUTPUTFILE

C...Format for reading lines.
      CHARACTER*6 STRFMT
      STRFMT='(A000)'
      WRITE(STRFMT(3:5),'(I3)') MAXLEN

C...Get input and output files
      WRITE(*,*)'File to be converted?'
      READ(*,*) INPUTFILE
      OPEN(unit=23,file=INPUTFILE,status='old')

      WRITE(*,*)'Output File?'
      READ(*,*) OUTPUTFILE
      OPEN(unit=24,file=OUTPUTFILE,status='new')
      
      
C...Loop until finds line beginning with "<event". 
  100 READ(23,STRFMT,END=130,ERR=130) STRING
  
  
      IBEG=0
  110 IBEG=IBEG+1
C...Allow indentation.
      IF(STRING(IBEG:IBEG).EQ.' '.AND.IBEG.LT.MAXLEN-5) GOTO 110 
      IF(STRING(IBEG:IBEG+5).NE.'<event') GOTO 100

      IFLAG=1
      DO WHILE(IFLAG.EQ.1)
        READ(23,'(2(I4,1x),4(E13.6,1x)',END=130,ERR=130)
     &        NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
        DO i=1,NUP
          READ(23,'(6(I4,1x),5(E13.6,1x),2(E10.1))',END=130,ERR=130)
     &  	    IDUP(i),ISTUP(i),MOTHUP(1,i),MOTHUP(2,i),
     &  	    ICOLUP(1,i),ICOLUP(2,i),(PUP(j,i),j=1,5),
     &  	    VTIMUP(i),SPINUP(i)
        ENDDO
C...Skip #pdf , <\event> and <event lines
          READ(23,'(A)',END=130,ERR=130) STRING
          READ(23,'(A)',END=130,ERR=130) STRING
          READ(23,'(A)',END=130,ERR=130) STRING
      
	WRITE(24,99999)(IDUP(i),i=1,8)
	WRITE(24,99998)(ISTUP(i),i=1,8)
	DO i=1,8
	  WRITE(24,99997)(PUP(j,i),j=1,4)
	ENDDO
	WRITE(24,99996)((ICOLUP(j,i),j=1,2),i=1,8)
        WRITE(24,99995)SCALUP,AQEDUP,AQCDUP
      ENDDO
      
99999 FORMAT('IDUP  ',8I4)
99998 FORMAT('ISTUP ',8I4)
99997 FORMAT(4E17.9)
99996 FORMAT('ICOLUP ',16I4)
99995 FORMAT('SCALUP ',E17.9,' AQEDUP ',E17.9,' AQCDUP ',E17.9)  

C...Error exit, typically when no more events.
  130 WRITE(*,*) ' Failed to read LHEF event information.'
      WRITE(*,*) ' Will assume end of file has been reached.'

      RETURN
      END

