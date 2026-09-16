C                                                                               
C                                                                               
      PROGRAM PSM_EE                                                            
C                                                                               
C*************************************************************                  
C                                                            *                  
C               PROGRAM: PROJECTED SHELL MODEL               *                  
C               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~               *                  
C        WITH 0, 2n, 2p AND 2n X 2p QP-CONFIGURATIONS        *                  
C                     (EVEN-EVEN NUCLEI)                     *                  
C                   MUNICH, AUG.'90, Y.SUN                   *                  
C                                                            *                  
C     ADOPTED TO PC JUL.'91 AND REVISED FEB.'94, Z.Y.ZHU     *                  
C          DEBUGGED AND REWRITTEN, SEP.'96, K.HARA           *                  
C                                                            *                  
C  NOTE:                                                     *                  
C  ~~~~~                                                     *                  
C  1. ORIGINAL FILES BY K.HARA AND S.IWASAKI, MUNICH '79:    *                  
C     HFOR2, HSUB2, HDAT2, HED2                              *                  
C                                                            *                  
C  2. THE "CDC NOS" EDIT PARAMETERS                          *                  
C       ($MD,$NDEG,$NBCS,$JMAX,$JOTMAX,$NRANK,               *                  
C        $NRANT,$KMAX,$MINT,$NATA)                           *                  
C     HAVE BEEN REPLACED BY INCLUDE PARAMETERS               *                  
C       (NMD,NNDEG,NNBCS,NJMAX,NJOTMAX,NNRANK,               *                  
C        NNRANT,NKMAX,NMINT,NNATA)                           *                  
C     WHICH DEFINE THE DIMENSIONS OF VARIOUS ARRAYS.         *                  
C                                                            *                  
C  3. THE DIRECT DISK ACCESS METHOD HAS BEEN REMOVED FROM    *                  
C     THE ORIGINAL PROGRAM. TESTED ON WORKSTATION AND PC.    *                  
C     REQUIRES AT LEAST 16MB RAM OR USE OF VIRTUAL MEMORY    *                  
C     FOR LESS RAM FOR PC, SEE THE DOCUMENT PSM.DOC.         *                  
C                                                            *                  
C  4. SIMILAR MODIFICATIONS HAVE BEEN DONE IN OTHER THREE    *                  
C     MUNICH PROGRAMS FOR OE, EO AND OO NUCLEI.              *                  
C                                                            *                  
C*************************************************************                  
C                                                                               
C     CALLs to DOSDAT and DOSTIM are specific to the NDP Fortran                
c      CALL DOSDAT(MONTH,MDAY,MYEAR,MWEEK)                                      
c      CALL DOSTIM(IHR,IMIN,ISEC,I100TH)                                        
C                                                                               
      CALL STORAG(1)                                                            
C                                                                               
      WRITE(*,15)                                                               
      WRITE(1,15)                                                               
      CALL DATA                                                                 
C                                                                               
c      WRITE(*,20)                                                              
c      WRITE(1,20)                                                              
      CALL GFV(127)                                                             
C                                                                               
c      WRITE(*,25)                                                              
c      WRITE(1,25)                                                              
      CALL LEGEND                                                               
C                                                                               
c      WRITE(*,30)                                                              
c      WRITE(1,30)                                                              
      CALL SPHSET                                                               
C                                                                               
      WRITE(*,35)                                                               
      WRITE(1,35)                                                               
      CALL NILBCS                                                               
C                                                                               
      WRITE(*,40)                                                               
      WRITE(1,40)                                                               
      CALL INTERA                                                               
C                                                                               
      WRITE(*,45)                                                               
      WRITE(1,45)                                                               
      CALL BASE00                                                               
C                                                                               
c      WRITE(*,*)                                                               
c      WRITE(*,*) '----------------------------------------------'              
c      WRITE(*,*) 'MODEL SPACE PREPARED: 0.[STOP] OR 1.[CONTINUE]'              
c    1 READ*, KEY                                                               
c      IF (KEY.GT.1) GOTO 1                                                     
c      WRITE(*,*) '----------------------------------------------'              
c      IF (KEY.EQ.0) GOTO 2                                                     
C                                                                               
      CALL STORAG(2)                                                            
C                                                                               
c      WRITE(*,50)                                                              
c      WRITE(1,50)                                                              
      CALL QPBASE                                                               
C                                                                               
      WRITE(*,55)                                                               
      WRITE(1,55)                                                               
      CALL KERNEL                                                               
C                                                                               
      WRITE(*,60)                                                               
      WRITE(1,60)                                                               
      CALL INTEG                                                                
C                                                                               
    2 CONTINUE                                                                  
C                                                                               
   10 FORMAT(/,20('*'), ' STORAG')                                              
   15 FORMAT(/,20('*'), ' DATA')                                                
   20 FORMAT(/,20('*'), ' GFV')                                                 
   25 FORMAT(/,20('*'), ' LEGEND')                                              
   30 FORMAT(/,20('*'), ' SPHSET')                                              
   35 FORMAT(/,20('*'), ' NILBCS')                                              
   40 FORMAT(/,20('*'), ' INTERA')                                              
   45 FORMAT(/,20('*'), ' BASE00')                                              
   50 FORMAT(/,20('*'), ' QPBASE')                                              
   55 FORMAT(/,20('*'), ' KERNEL')                                              
   60 FORMAT(/,20('*'), ' INTEG')                                               
C                                                                               
c      CALL DOSTIM(JHR,JMIN,JSEC,J100TH)                                        
c      WRITE(1,100) MONTH,MDAY,MYEAR,MWEEK                                      
c      WRITE(1,200) IHR,IMIN,ISEC,I100TH,JHR,JMIN,JSEC,J100TH                   
C                                                                               
  100 FORMAT(//' EXECUTION DAY: ',I2,'-',I2,'-',I4,10X,                         
     +       ' DAY OF THE WEEK:',I2)                                            
  200 FORMAT(/,'    STARTED AT ',I2,':',I2,':',I2,'.',I2,                       
     +     '  AND  FINISHED AT ',I2,':',I2,':',I2,'.',I2)                       
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE STORAG(IO)                                                     
C                                                                               
C**************************                                                     
C     OPEN I/O FILES      *                                                     
C**************************                                                     
C                                                                               
      CHARACTER*31 FILE0,FILE1,FILE2,FILE3,FILE4                                
C                                                                               
      SAVE FILE2,FILE3,FILE4                                                    
C                                                                               
      COMMON/NHELEM/NHOUT                                                       
C                                                                               
      IF (IO.GT.1) GOTO 6                                                       
C                                                                               
C     Default I/O file names:                                                   
      FILE0='EE_DATA'                                                           
      FILE1='EE_OUT1'                                                           
      FILE2='EE_OUT2'                                                           
      FILE3='EE_OUT3'                                                           
      FILE4='EE_OUT4'                                                           
C     Other I/O file names can be used by setting INP=1; 0 for default          
      INP=0                                                                     
      IF (INP.EQ.0) GOTO 5                                                      
      WRITE(*,*)                                                                
      WRITE(*,*) '----------------------------------------------'               
      WRITE(*,*) 'ENTER THE NAME OF INPUT FILE: [PATH\]FILE_NAME'               
      READ(*,10) FILE0                                                          
   10 FORMAT(A31)                                                               
      WRITE(*,*) '----------------------------------------------'               
      L=31                                                                      
      DO 1 M=1,31                                                               
      IF (FILE0(L:L).NE.' ') GOTO 2                                             
      L=L-1                                                                     
    1 CONTINUE                                                                  
      STOP                                                                      
    2 CONTINUE                                                                  
      K=0                                                                       
      DO 3 M=1,4                                                                
      IF (FILE0(L-K:L-K).EQ.'.') GOTO 4                                         
      K=K+1                                                                     
    3 CONTINUE                                                                  
      K=-2                                                                      
    4 CONTINUE                                                                  
      L=L-K-1                                                                   
      FILE1=FILE0(1:L)//'.O1 '                                                  
      FILE2=FILE0(1:L)//'.O2 '                                                  
      FILE3=FILE0(1:L)//'.O3 '                                                  
      FILE4=FILE0(1:L)//'.O4 '                                                  
    5 CONTINUE                                                                  
C                                                                               
      OPEN(4,FILE=FILE0,STATUS='OLD')                                           
C                                                                               
      OPEN(1,FILE=FILE1,STATUS='UNKNOWN')                                       
      CLOSE(1,STATUS='DELETE')                                                  
      OPEN(1,FILE=FILE1,STATUS='NEW')                                           
      WRITE(1,20) FILE0                                                         
   20 FORMAT(/,' INPUT FILE: ',A31)                                             
C                                                                               
      RETURN                                                                    
C                                                                               
    6 CLOSE(4)                                                                  
      OPEN(2,FILE=FILE2,STATUS='UNKNOWN')                                       
      OPEN(3,FILE=FILE3,STATUS='UNKNOWN')                                       
      CLOSE(2,STATUS='DELETE')                                                  
      CLOSE(3,STATUS='DELETE')                                                  
      OPEN(2,FILE=FILE2,STATUS='NEW')                                           
      OPEN(3,FILE=FILE3,STATUS='NEW')                                           
C     Output hamiltonian and norm matrix (NHOUT=0 for no, 1 for yes)            
      NHOUT=0                                                                   
      IF (NHOUT.NE.0) THEN                                                      
       OPEN(4,FILE=FILE4,STATUS='UNKNOWN')                                      
       CLOSE(4,STATUS='DELETE')                                                 
       OPEN(4,FILE=FILE4,STATUS='NEW')                                          
       WRITE(4,30)                                                              
   30  FORMAT(' NORM AND HAMILTONIAN MATRIX ELEMENTS')                          
      ENDIF                                                                     
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE DATA                                                           
C                                                                               
C====================================================================           
C                                                                               
C NPNUMB(2) NUMBER OF NEUTRONS AND PROTONS                                      
C EPS(2)    DEFORMATION PARAMETERS EPS2 AND EPS4                                
C                                                                               
C AG(2)     MONOPOLE PAIRING FORCE STRENGTH FOR NEUTRONS AND PROTONS            
C AGQ(2)    QUADRUPOLE PAIRING FORCE STRENGTH FOR NEUTRONS AND PROTONS          
C                                                                               
C EP(NNBCS,2)   NILSSON ENERGIES FOR NEUTRON AND PROTON                         
C MXSH(2)       NUMBER OF MAJOR SHELLS FOR NEUTRON AND PROTRON                  
C NILS(3,2)     DIMENSION OF THE SHELL FOR NEUTRON AND PROTRON                  
C NPSH(3,2)     MAJOR QUANTUM NUMBER N FOR THE SHELLS                           
C                                                                               
C====================================================================           
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /BCSBCS/ ANP(2),ALAMB(2),DELTA(2),EP(NNBCS,2),                     
     1                HSP(2),QQ0(2),HH0(2),QP0(2),LNP(2)                        
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /GAUSSL/ X(NNDEG),Y(NNDEG),WEIT(NNDEG),NDEG                        
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /JJJJJJ/ JJ(NJOTMX),JOTMAX                                         
      COMMON /MAXMAX/ EMAX(3),KMAX(3),NMAX(3),JMAX                              
      COMMON /NUCDAT/ AMAS,EPS(2),NPNUMB(2),NPCORE(2),LSF                       
      COMMON /PARAM / AG(2),AGQ(2),DELE,HOLNG2(2),QCL,HFB,                      
     1                AX(3),AXH(3),XQQ(2),GQP(2)                                
C                                                                               
      DIMENSION GNP(2),NPC(6)                                                   
C                                                                               
      DATA NPC/2,8,20,40,70,112/                                                
      DELE=0.5                                                                  
C                                                                               
C------------------------------------------                                     
C     READ DATA SKIPPING THE FIRST LINE                                         
C------------------------------------------                                     
C                                                                               
      READ(4,*)                                                                 
      READ(4,101) NPNUMB                                                        
      READ(4,102) EPS                                                           
      READ(4,101) NMAX                                                          
      READ(4,101) KMAX                                                          
      READ(4,102) EMAX                                                          
      DO 1 IT=1,2                                                               
      READ(4,101) MS,(NPSH(N,IT),N=1,MS),NSINT                                  
      MXSH(IT)=MS                                                               
      NNINT(IT)=NPSH(NSINT,IT)                                                  
      MIXINT(IT)=NSINT+(IT-1)*MXSH(1)                                           
    1 CONTINUE                                                                  
      READ(4,101) JINI,JFIN,JDIF                                                
      READ(4,102) GNP                                                           
      READ(4,102) GQ                                                            
      READ(4,101) NDEG                                                          
      READ(4,101) LSF                                                           
C                                                                               
      MASS=NPNUMB(1)+NPNUMB(2)                                                  
      AMAS=MASS                                                                 
C                                                                               
      M=JINI                                                                    
      N=JFIN                                                                    
      JINI=MIN(M,N)                                                             
      JFIN=MAX(M,N)                                                             
      JOTMAX=(JFIN-JINI)/JDIF+1                                                 
C                                                                               
      B=JFIN                                                                    
      M=NINT(1.3*B)                                                             
      N=2*(MAX(M,NDEG)/2)                                                       
      NDEG=N                                                                    
C                                                                               
      NBCS=0                                                                    
      DO 5 IT=1,2                                                               
      ND=NNINT(IT)+1                                                            
      ND=(ND*(ND+1))/2                                                          
      NLINT(IT)=ND                                                              
      ITINT(IT)=IT                                                              
      NPBCS=0                                                                   
      NNMAX=0                                                                   
      NNMIN=100                                                                 
      MS=MXSH(IT)                                                               
      DO 3 N=1,MS                                                               
      NPBCS=NPBCS+(NPSH(N,IT)+1)*(NPSH(N,IT)+2)/2                               
      NNMIN=MIN(NPSH(N,IT),NNMIN)                                               
      NNMAX=MAX(NPSH(N,IT),NNMAX)                                               
    3 CONTINUE                                                                  
      NBCS=MAX(NPBCS,NBCS)                                                      
      NPCORE(IT)=0                                                              
      IF (NNMIN.GT.0) NPCORE(IT)=NPC(NNMIN)                                     
C                                                                               
      B=NPNUMB(IT)-NPNUMB(3-IT)                                                 
      AG(IT)=(GNP(1)-GNP(2)*B/AMAS)/AMAS                                        
C                                                                               
C     ANOTHER POSSIBLE FORM OF PAIRING STRENGTH                                 
C                                                                               
C      AG(IT)=GNP(IT)/AMAS                                                      
C                                                                               
      AGQ(IT)=GQ*AG(IT)                                                         
    5 CONTINUE                                                                  
C                                                                               
C------------------------------------------                                     
C     WRITE DOWN THE INPUT DATA                                                 
C------------------------------------------                                     
C                                                                               
      WRITE(1,*)                                                                
      WRITE(1,208) NPNUMB,MASS,NPCORE                                           
      WRITE(1,209) EPS                                                          
      WRITE(1,205) NMAX                                                         
      WRITE(1,206) KMAX,EMAX                                                    
      DO 6 IT=1,2                                                               
      MS=MXSH(IT)                                                               
      IF (IT.EQ.1) THEN                                                         
       WRITE(1,201) MS,(NPSH(N,IT),N=1,MS),NNINT(IT),NLINT(IT)                  
      ELSE                                                                      
       WRITE(1,202) MS,(NPSH(N,IT),N=1,MS),NNINT(IT),NLINT(IT)                  
      ENDIF                                                                     
    6 CONTINUE                                                                  
      WRITE(1,207) JINI,JFIN,JDIF,JOTMAX                                        
      WRITE(1,210) GNP,GQ                                                       
      WRITE(1,211) NDEG                                                         
      WRITE(1,212) LSF                                                          
      WRITE(1,213) DELE                                                         
C                                                                               
C---------------------------------------                                        
C     CONSISTENCY CHECK                                                         
C---------------------------------------                                        
C                                                                               
      WRITE(1,*)                                                                
      IF (NPNUMB(1).LE.2.OR.NPNUMB(2).LE.2) THEN                                
        WRITE(1,100)                                                            
        WRITE(*,100)                                                            
  100   FORMAT(' STOP: PARTICLE NUMBER TOO SMALL')                              
        STOP                                                                    
      ENDIF                                                                     
C                                                                               
      NERR=0                                                                    
      NP1=MOD(NPNUMB(1),2)                                                      
      NP2=MOD(NPNUMB(2),2)                                                      
      IF (NP1.EQ.0.AND.NP2.EQ.0) GOTO 7                                         
      WRITE(1,*) ' INCONSISTENCY IN INPUT DATA: NPNUMB(N,Z)',NPNUMB             
      NERR=NERR+1                                                               
    7 CONTINUE                                                                  
C                                                                               
      IF (NDEG.GT.NNDEG) THEN                                                   
        WRITE(1,*) ' SET NNDEG IN EE_PARA FILE TO',NDEG                         
        NERR=NERR+1                                                             
      ENDIF                                                                     
C                                                                               
      IF (NBCS.GT.NNBCS) THEN                                                   
        WRITE(1,*) ' SET NNBCS IN EE_PARA FILE TO',NBCS                         
        NERR=NERR+1                                                             
      ENDIF                                                                     
C                                                                               
      MS1=MXSH(1)                                                               
      MS2=MXSH(2)                                                               
      JMAX=MAX(NPSH(MS1,1),NPSH(MS2,2))+1                                       
      IF (JMAX.GT.NJMAX) THEN                                                   
        WRITE(1,*) ' SET NJMAX IN EE_PARA FILE TO',JMAX                         
        NERR=NERR+1                                                             
      ENDIF                                                                     
C                                                                               
      MD=JMAX*(JMAX+1)/2                                                        
      IF (MD.GT.NMD) THEN                                                       
        WRITE(1,*) ' SET NMD IN EE_PARA FILE TO',MD                             
        NERR=NERR+1                                                             
      ENDIF                                                                     
C                                                                               
      NN=NMAX(1)                                                                
      KK=KMAX(1)                                                                
      DO 8 N=2,3                                                                
      IF (NMAX(N).GT.NN) NN=NMAX(N)                                             
      IF (KMAX(N).GT.KK) KK=KMAX(N)                                             
    8 CONTINUE                                                                  
C                                                                               
      IF (NN.GT.NNRANK) THEN                                                    
        WRITE(1,*) ' SET NNRANK IN EE_PARA FILE TO',NN                          
        NERR=NERR+1                                                             
      ENDIF                                                                     
C                                                                               
      IF (KK.GT.NKMAX) THEN                                                     
        WRITE(1,*) ' SET NKMAX IN EE_PARA FILE TO',KK                           
        NERR=NERR+1                                                             
      ENDIF                                                                     
C                                                                               
      IF (JOTMAX.GT.NJOTMX) THEN                                                
        WRITE(1,*) ' SET NJOTMX IN EE_PARA FILE TO',JOTMAX                      
        NERR=NERR+1                                                             
      ENDIF                                                                     
C                                                                               
      IF (NERR.GT.0) THEN                                                       
        WRITE(*,*) ' STOP: ERROR SUMMARY WRITTEN TO OUTPUT FILE 1'              
        STOP                                                                    
      ENDIF                                                                     
C                                                                               
      N=0                                                                       
      DO 9 M=JINI,JFIN,JDIF                                                     
      N=N+1                                                                     
      JJ(N)=M                                                                   
    9 CONTINUE                                                                  
C                                                                               
  101 FORMAT(20X,6I5)                                                           
  102 FORMAT(20X,3F10.5)                                                        
  201 FORMAT(' SINGLE-PARTICLE SPACE NEUTRONS:',I2,' SHELLS',3I3,               
     1       '  INTRUDER',I2,' LEVELS',I3)                                      
  202 FORMAT('                        PROTONS:',I2,' SHELLS',3I3,               
     1       '  INTRUDER',I2,' LEVELS',I3)                                      
  205 FORMAT(' LIMIT OF QP STATES: NMAX =',3I4,6X,'(2N, 2P, 2NX2P)')            
  206 FORMAT(' WINDOWS TO K AND E: KMAX =',3I4,6X,'EMAX =',3F6.2)               
  207 FORMAT(' RANGE OF SPIN: JINI, JFIN, JDIF =',3I4,6X,'JOTMAX =',I3)         
  208 FORMAT(' PARTICLE NUMBERS (N,Z,A):',3I4,7X,'CORE (N,Z):',2I4)             
  209 FORMAT(' DEFORMATION PARAMETERS: EPS2 = ',F6.3,7X,'EPS4 = ',F6.3)         
  210 FORMAT(' PAIRING: G1, G2 =',2F7.2,5X,'G(QUADRUPOLE) =',F5.2,              
     1       '*G(MONOPOLE)')                                                    
  211 FORMAT(' BETA-MESH: NDEG =',I3,3X,'A RULE OF THUMB: ',                    
     1       ' NDEG .GT. 1.3*JFIN')                                             
  212 FORMAT(' LS FORCE: LSF =',I2,3X,'1 FOR NILSSON ET AL',                    
     1       ' AND 2 FOR ZHANG ET AL')                                          
  213 FORMAT(' EFFECTIVE CHARGE: DELE =',2F5.2)                                 
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE NILSON(NSHELL,IT)                                              
C                                                                               
C***********************************************************                    
C     CALCULATION OF NILSSON SINGLE PARTICLE LEVELS        *                    
C***********************************************************                    
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      PARAMETER (NJMAX1=NJMAX+1)                                                
      PARAMETER (NIE=127)                                                       
C                                                                               
      COMMON /GAUSSL/ X(NNDEG),Y(NNDEG),WEIT(NNDEG),NDEG                        
      COMMON /GFVIV / IV(-NIE:NIE)                                              
      COMMON /GFVSQ / SQ(0:NIE)                                                 
      COMMON /LEG   / PLEG(NNDEG)                                               
      COMMON /NIL   / VAL(NMD),VEC(NMD,NMD),RSQ(NJMAX1,NJMAX1),                 
     1                W0(2),EPS1,EPS2,NS,ND                                     
      COMMON /NUCDAT/ AMAS,EPS(2),NPNUMB(2),NPCORE(2),LSF                       
      COMMON /SPHJLM/ J1(NMD),L1(NMD),M1(NMD)                                   
      COMMON /SPH   / ES(NMD),Q(NMD,NMD,5),R(NMD,NMD,3),P(NMD,NMD),             
     1                SP(NMD,NMD,3)                                             
C                                                                               
      DIMENSION AKA(2),BKA(2),AMU(2),BMU(2),CKA(0:8,2),CMU(0:8,2),              
     1          FKA(3,2),FMU(3,2),KAMU(2),Z(NMD)                                
C                                                                               
      SAVE AKA,BKA,AMU,BMU,CKA,CMU,FKA,FMU,KAMU,W                               
C                                                                               
      DATA W/41.4678/,                                                          
C                                                                               
C CHOICE OF THE LS FORCE PARAMETERS KAPPA AND MU (CONTROLLED BY LSF)            
C                                                                               
C LSF=1: NILSSON ET AL FOR RARE-EARTH NUCLEI, NUCL. PHYS. A131 (1969) 1         
     1   AKA/0.06433, 0.07441/, BKA/0.004D-3, 0.068D-3/,                        
     2   AMU/0.59010, 0.51180/, BMU/1.080D-3,-0.560D-3/,                        
C LSF=2: ZHANG ET AL FOR GENERAL NUCLEI, J. PHYS. G13 (1987) L75                
     1   CKA/.12, .12, .105, .09, .070, .062, .062, .062, .062,                 
     2       .12, .12, .105, .09, .065, .060, .054, .054, .054/,                
     3   CMU/.00, .00, .00, .25, .39, .43, .34, .26, .26,                       
     4       .00, .00, .00, .30, .57, .54, .52, .69, .69/,                      
C FOR FINE TUNING OF KAPPA AND MU, FACTORS MAY BE MULTIPLIED:                   
     5   FKA/3*1.0,3*1.0/,FMU/3*1.0,3*1.0/,KAMU/2*0/                            
C                                                                               
      NS=NSHELL+1                                                               
      ND=(NS*(NS+1))/2                                                          
      II=1                                                                      
      IF (IV(NSHELL).NE.+1) II=2                                                
      DO 10 I=1,ND                                                              
   10 L1(I)=NS-((NS-J1(I))/2)*2                                                 
      DO 15 I1=II,NS,2                                                          
      DO 15 I2=II,I1,2                                                          
      RSQ(I2,I1)=0.                                                             
      LL=I1-I2                                                                  
      IF (LL.GT.4) GOTO 15                                                      
      IF (LL-2) 11,12,13                                                        
   11 RSQ(I2,I1)=NS+0.5                                                         
      GOTO 15                                                                   
   12 RSQ(I2,I1)=SQ(NS-I2)*SQ(NS+I2+1)                                          
      GOTO 15                                                                   
   13 RSQ(I2,I1)=(NS+3.5+2.*I2)                                                 
     1   *SQ(NS-I2)*SQ(NS-I2-2)                                                 
     2   /(SQ(NS+I2+1)*SQ(NS+I2+3))                                             
   15 RSQ(I1,I2)=RSQ(I2,I1)                                                     
      B=NPNUMB(IT)-NPNUMB(3-IT)                                                 
      C=(1.+B/AMAS)**(1./3.)                                                    
      WC=W*C/AMAS**(1./3.)                                                      
      EPS1=EPS(1)*2./3.                                                         
      EPS2=EPS(2)                                                               
      NDEH=NDEG/2                                                               
      SS=0.                                                                     
      DO 1 I=1,NDEH                                                             
      COSBET=X(I)                                                               
      CALL LEGPOL(5,COSBET)                                                     
    1 SS=SS+Y(I)/(1.-EPS1*PLEG(3)+2.*EPS2*PLEG(5))**1.5                         
      SS=(SS/(1.+0.5*EPS1)/SQRT(1.-EPS1))**(1./3.)                              
    2 W0(IT)=WC*SS                                                              
C                                                                               
      IF (LSF.EQ.1) THEN                                                        
        A=AKA(IT)-BKA(IT)*AMAS                                                  
        B=AMU(IT)-BMU(IT)*AMAS                                                  
      ELSE                                                                      
        A=CKA(NSHELL,IT)                                                        
        B=CMU(NSHELL,IT)                                                        
      ENDIF                                                                     
      IF (IT.EQ.1) THEN                                                         
       WRITE(1,201) NSHELL,A,B                                                  
      ELSE                                                                      
       WRITE(1,202) NSHELL,A,B                                                  
      ENDIF                                                                     
  201 FORMAT(' NEUTRON SHELL, KAPPA AND MU: ',I3,2F10.5)                        
  202 FORMAT(' PROTON  SHELL, KAPPA AND MU: ',I3,2F10.5)                        
C                                                                               
      KAMU(IT)=KAMU(IT)+1                                                       
      IF (FKA(KAMU(IT),IT).NE.1.0.OR.FMU(KAMU(IT),IT).NE.1.0) THEN              
        A=A*FKA(KAMU(IT),IT)                                                    
        B=B*FMU(KAMU(IT),IT)                                                    
        WRITE(1,203) FKA(KAMU(IT),IT),FMU(KAMU(IT),IT)                          
        WRITE(1,204) A,B                                                        
      ENDIF                                                                     
  203 FORMAT('   FINE TUNING FACTORS      :',2F10.5)                            
  204 FORMAT('   KAPPA AND MU (EFFECTIVE) :',2F10.5)                            
C                                                                               
      A=A*WC                                                                    
      C=NSHELL+1.5                                                              
      D=W0(IT)*C+0.5*A*(1.5+B*(C-1.5)*(C+1.5))                                  
      B=B-1.                                                                    
      F=EPS1*W0(IT)                                                             
      E=EPS2*W0(IT)                                                             
C                                                                               
      DO 8 K=1,ND                                                               
      L = L1(K) - 1                                                             
      J = 2 * J1(K) - 1                                                         
      M = 2 * M1(K) - 1                                                         
      SJ=J1(K)-0.5                                                              
      LS=L1(K)-1                                                                
      SL=LS                                                                     
      SE=D-A*(SJ*(SJ+1.)+B*SL*(SL+1.))                                          
      ES(K)=SE                                                                  
      DO 7 L=K,ND                                                               
      SH=0.                                                                     
      IF (M1(K).NE.M1(L)) GOTO 7                                                
      JJ=J1(K)-J1(L)                                                            
      IF (JJ.GT.4) GOTO 7                                                       
      SH=RSQ(L1(K),L1(L))*E*P(K,L)                                              
      IF (JJ.GT.2) GOTO 7                                                       
      SH=SH-RSQ(L1(K),L1(L))*F*Q(K,L,3)                                         
    7 VEC(L,K)=SH                                                               
    8 VEC(K,K)=VEC(K,K)+SE                                                      
C                                                                               
      CALL RSDIAG(NMD,ND,VEC,VAL,Z,1)                                           
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE BCS(IT)                                                        
C                                                                               
C CALLED BY SUB. NILBCS                                                         
C                                                                               
C********************************************************                       
C                                                       *                       
C  IT    1  NEUTRON                                     *                       
C        2  PROTON                                      *                       
C                                                       *                       
C  MODE  0  SPHERICAL NUCLEUS                           *                       
C        1  DEFORMED NUCLEUS                            *                       
C                                                       *                       
C  INOUT 0  INITIAL VALUES ARE CALCULATED               *                       
C        1  INITIAL VALUES ARE GIVEN EXTERNALLY         *                       
C           INOUT RETURNS THE NUMBER OF ITERATIONS      *                       
C                 WHEN BCS CONVERGED ELSE 0             *                       
C                                                       *                       
C  L  NUMBER OF LEVELS                                  *                       
C  N  NUMBER OF PARTICLES                               *                       
C                                                       *                       
C  C  LAMBDA                                            *                       
C  D  DELTA                                             *                       
C  G  FORCE CONSTANT FOR THE MONOPOLE PAIRING           *                       
C                                                       *                       
C  EP(L,IT)     SINGLE PARTICLE ENERGIES                *                       
C  OMEGA(L,IT)  PAIR DEGENERACY [FOR SPHERICAL CASE]    *                       
C                                                       *                       
C********************************************************                       
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /BCSBCS/ ANP(2),ALAMB(2),DELTA(2),EP(NNBCS,2),                     
     1                HSP(2),QQ0(2),HH0(2),QP0(2),LNP(2)                        
      COMMON /PARAM / AG(2),AGQ(2),DELE,HOLNG2(2),QCL,HFB,                      
     1                AX(3),AXH(3),XQQ(2),GQP(2)                                
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /UVUV  / ESP(NMD,6),EQP(NMD,6),U(NMD,6),V(NMD,6),                  
     1                KQN(NMD,6)                                                
      COMMON /DEF   / QEL(NMD,NMD,5,6),PEL(NMD,NMD,6),                          
     1                SEL(NMD,NMD,3,6),REL(NMD,NMD,3,6),                        
     2                HEL(NMD,NMD,6),AMP(NMD,NMD,6)                             
C                                                                               
      DATA IMAX/100/, ERR/1.D-12/                                               
C                                                                               
      INOUT=0                                                                   
      MODE=1                                                                    
      L=LNP(IT)                                                                 
      N=ANP(IT)+0.005                                                           
C                                                                               
      C=ALAMB(IT)                                                               
      D=DELTA(IT)                                                               
      G=AG(IT)                                                                  
      G2=0.5*G                                                                  
      PN=N                                                                      
      OSUM=0.                                                                   
      DO 1 I= 1,L                                                               
      OMEG= 1.                                                                  
C      OMEG= OMEGA(I,IT)                                                        
      OSUM= OSUM+OMEG                                                           
    1 CONTINUE                                                                  
C                                                                               
      IF (INOUT.GT.0) GOTO 6                                                    
      RHO= G2*OSUM                                                              
      CHI= 1.-PN/OSUM                                                           
C      IF (MODE.NE.0) GOTO 3                                                    
C ***  FIND TRIAL VALUES FOR SPHERICAL NUCLEUS                                  
C      EPAV= 0.                                                                 
C      EPSQ= 0.                                                                 
C      DO 2 I =1,L                                                              
C      OMEP= OMEGA(I,IT)*EP(I,IT)                                               
C      EPAV= EPAV+OMEP                                                          
C    2 EPSQ= EPSQ+OMEP*EP(I,IT)                                                 
C      EPAV= EPAV/OSUM                                                          
C      EPSQ= EPSQ/OSUM-EPAV*EPAV                                                
C      COSHTH= 1.+EPSQ/RHO**2                                                   
C      GOTO 5                                                                   
C ***  FIND TRIAL VALUES FOR DEFORMED NUCLEUS                                   
C    3 CONTINUE                                                                 
      EPMIN= EP(1,IT)                                                           
      EPMAX= EP(L,IT)                                                           
      DO 4 I= 1,L                                                               
      IF (EPMIN.GT.EP(I,IT)) EPMIN=EP(I,IT)                                     
      IF (EPMAX.LT.EP(I,IT)) EPMAX=EP(I,IT)                                     
    4 CONTINUE                                                                  
      EPAV= 0.5*(EPMAX+EPMIN)                                                   
      EPSQ= 0.5*(EPMAX-EPMIN)                                                   
      THETA= EPSQ/RHO                                                           
      RHO= EPSQ/SINH(THETA)                                                     
      COSHTH= COSH(THETA)                                                       
C ***  START ITERATION  ***                                                     
    5 CONTINUE                                                                  
      C= EPAV-RHO*CHI*COSHTH                                                    
      D= RHO*SQRT(1.-CHI**2)                                                    
    6 CONTINUE                                                                  
      DO 8 ITER=1,IMAX                                                          
      FA= 0.                                                                    
      FB= 0.                                                                    
      FC= 0.                                                                    
      FD= 0.                                                                    
      DO 7 I=1,L                                                                
      E= EP(I,IT)-C                                                             
      ESQ= E*E+D*D                                                              
      OMEG= 1.                                                                  
C      OMEG= OMEGA(I,IT)                                                        
      E1= OMEG/SQRT(ESQ)                                                        
      E3= E1/ESQ                                                                
      FA= FA+E1                                                                 
      FB= FB+E1*E                                                               
      FC= FC+E3*E                                                               
      FD= FD+E3                                                                 
    7 CONTINUE                                                                  
      FB= OSUM-FB                                                               
      FD= FD*D                                                                  
      A= 1.-G2*FA                                                               
      B= 1.-FB/PN                                                               
      ABSA= ABS(A)                                                              
      ABSB= ABS(B)                                                              
      IF (ABSA+ABSB.LT.ERR) GOTO 9                                              
      A= A*(D/G2)                                                               
      B= B*PN                                                                   
      DET= (FC*FC+FD*FD)*D                                                      
      CER= (A*FC+B*FD)/DET                                                      
      DER= (B*FC-A*FD)/DET                                                      
      C= C+CER                                                                  
      D= D+DER                                                                  
      IF (D.LT.0.) D = -D                                                       
    8 CONTINUE                                                                  
      GOTO 10                                                                   
    9 INOUT=ITER                                                                
      ALAMB(IT) = C                                                             
      DELTA(IT) = D                                                             
      SN = 0.                                                                   
      HSP(IT) = 0.                                                              
      QQ0(IT) = 0.                                                              
      HH0(IT) = 0.                                                              
      QP0(IT) = 0.                                                              
      MSHL=MXSH(IT)                                                             
      DO 16 NSHL=1,MSHL                                                         
      MIX = NSHL + (IT-1)*MXSH(1)                                               
      NL  = NILS(NSHL,IT)                                                       
      DO 15 I=1,NL                                                              
      E         = ESP(I,MIX)-C                                                  
      EQP(I,MIX)= SQRT(E*E+D*D)                                                 
      VV        = 1.-E/EQP(I,MIX)                                               
      VISQ      = 0.5*VV                                                        
      U(I,MIX)  = SQRT(1.-VISQ)                                                 
      V(I,MIX)  = SQRT(VISQ)                                                    
      VU        = V(I,MIX)*U(I,MIX)                                             
      SN        = SN+VV                                                         
      HSP(IT)   = HSP(IT)+VV*HEL(I,I,MIX)                                       
      QQ0(IT)   = QQ0(IT)+VV*QEL(I,I,3,MIX)                                     
      HH0(IT)   = HH0(IT)+VV*PEL(I,I,MIX)                                       
      QP0(IT)   = QP0(IT)+VU*QEL(I,I,3,MIX)                                     
   15 CONTINUE                                                                  
   16 CONTINUE                                                                  
      WRITE(1,101) IT,ALAMB(IT),DELTA(IT),SN                                    
  101 FORMAT(/' BCS: IT =',I2,'  LAMBDA =',F8.4,'  DELTA =',F8.4,               
     1        '  N (SUM 2.0*V(I)**2) =',F8.4)                                   
      RETURN                                                                    
C                                                                               
   10 INOUT=0                                                                   
      PRINT*, 'STOP: BCS FAILED',ABSA,ABSB                                      
      STOP                                                                      
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE NILBCS                                                         
C                                                                               
C CALLED BY MAIN                                                                
C                                                                               
C******************************************************************             
C        CALCULATION OF THE NILSON-BCS-BASIS                      *             
C        FOR NEUTRONS (IT=1) AND PROTONS (IT=2)                   *             
C******************************************************************             
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (NJMAX1=NJMAX+1)                                                
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /BCSBCS/ ANP(2),ALAMB(2),DELTA(2),EP(NNBCS,2),                     
     1                HSP(2),QQ0(2),HH0(2),QP0(2),LNP(2)                        
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /PARAM / AG(2),AGQ(2),DELE,HOLNG2(2),QCL,HFB,                      
     1                AX(3),AXH(3),XQQ(2),GQP(2)                                
      COMMON /NIL   / VAL(NMD),VEC(NMD,NMD),RSQ(NJMAX1,NJMAX1),                 
     1                W0(2),EP1,EP2,NS,ND                                       
      COMMON /NUCDAT/ AMAS,EPS(2),NPNUMB(2),NPCORE(2),LSF                       
C                                                                               
      DATA W/41.4678/                                                           
C                                                                               
c      WRITE(*,*)                                                               
c      WRITE(*,*) '----------------------------------------------'              
c      WRITE(*,*) 'LS FORCE: 1.[NILSSON ET AL] OR 2.[ZHANG ET AL]'              
c    4 READ*, LSF                                                               
c      IF (LSF.LT.1.OR.LSF.GT.2) GOTO 4                                         
c      WRITE(*,*) '----------------------------------------------'              
C                                                                               
      DO 1 IT=1,2                                                               
C                                                                               
      WRITE(1,*)                                                                
      IF (LSF.EQ.1) THEN                                                        
       WRITE(1,*) 'LS FORCE: NILSSON ET AL'                                     
      ELSE                                                                      
       WRITE(1,*) 'LS FORCE: ZHANG ET AL'                                       
      ENDIF                                                                     
C                                                                               
      MS=MXSH(IT)                                                               
      II=0                                                                      
      DO 2 N=1,MS                                                               
      NSHELL=NPSH(N,IT)                                                         
      CALL NILSON(NSHELL,IT)                                                    
      NL=ND                                                                     
      NILS(N,IT)=NL                                                             
      DO 3 I = 1,NL                                                             
      II = II + 1                                                               
    3 EP(II,IT) = VAL(I)                                                        
    2 CALL DEFSET(IT,N)                                                         
C                                                                               
      HOLNG2(IT)=W/W0(IT)                                                       
      LNP(IT)=II                                                                
      ANP(IT)=NPNUMB(IT)-NPCORE(IT)                                             
      CALL BCS(IT)                                                              
C                                                                               
    1 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE AINVBG(NDIM,A,B,ND,MD,INV,DET)                                 
C                                                                               
C CALLED BY SUB. CONTRA                                                         
c                                                                               
C*************************************************************                  
C       B=(1/A)*B       GAUSS ELIMINATION METHOD             *                  
C*************************************************************                  
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM)                                       
C                                                                               
      DATA F/1.D-32/                                                            
C                                                                               
      INV=1                                                                     
      E=0.D00                                                                   
      DO 1 I=1,ND                                                               
      DO 1 J=1,ND                                                               
      AIK=A(I,J)                                                                
C991  AIK=CONJ(AIK)                                                             
      E=E+A(I,J)*AIK                                                            
    1 CONTINUE                                                                  
      D=ND*ND                                                                   
      E=E*F/D                                                                   
      IF (ND.GT.1) GOTO 10                                                      
      DET=A(1,1)                                                                
      DO 5 L=1,MD                                                               
      B(1,L)=B(1,L)/DET                                                         
    5 CONTINUE                                                                  
      RETURN                                                                    
   10 DET=1.D00                                                                 
      N1=ND-1                                                                   
      DO 40 K=1,N1                                                              
      KK=K                                                                      
      K1=K+1                                                                    
      AIK=A(K,K)                                                                
C992  AIK=CONJ(AIK)                                                             
      C=A(K,K)*AIK                                                              
      DO 15 I=K1,ND                                                             
      AIK=A(I,K)                                                                
C993  AIK=CONJ(AIK)                                                             
      D=A(I,K)*AIK                                                              
      IF (C.GE.D) GOTO 15                                                       
      C=D                                                                       
      KK=I                                                                      
   15 CONTINUE                                                                  
      IF (KK.NE.K) INV=INV+1                                                    
      IF (C.LT.E) GOTO 60                                                       
      SUM=A(KK,K)                                                               
      A(KK,K)=A(K,K)                                                            
      A(K,K)=SUM                                                                
      DET=DET*SUM                                                               
      DO 20 I=K1,ND                                                             
      AIK=A(KK,I)/SUM                                                           
      A(KK,I)=A(K,I)                                                            
      A(K,I)=AIK                                                                
   20 CONTINUE                                                                  
      DO 25 L=1,MD                                                              
      AIK=B(KK,L)/SUM                                                           
      B(KK,L)=B(K,L)                                                            
      B(K,L)=AIK                                                                
   25 CONTINUE                                                                  
      DO 35 I=K1,ND                                                             
      AIK=A(I,K)                                                                
      DO 30 J=K1,ND                                                             
      A(I,J)=A(I,J)-AIK*A(K,J)                                                  
   30 CONTINUE                                                                  
      DO 35 L=1,MD                                                              
      B(I,L)=B(I,L)-AIK*B(K,L)                                                  
   35 CONTINUE                                                                  
   40 CONTINUE                                                                  
      SUM=A(ND,ND)                                                              
      AIK=SUM                                                                   
C994  AIK=CONJ(AIK)                                                             
      D=SUM*AIK                                                                 
      IF (D.LT.E) GOTO 60                                                       
      DET=DET*SUM                                                               
      IF (MOD(INV,2).EQ.0) DET=-DET                                             
      DO 45 L=1,MD                                                              
      B(ND,L)=B(ND,L)/SUM                                                       
   45 CONTINUE                                                                  
      DO 55 L=1,MD                                                              
      DO 55 K=1,N1                                                              
      I=ND-K                                                                    
      I1=I+1                                                                    
      SUM=B(I,L)                                                                
      DO 50 J=I1,ND                                                             
      SUM=SUM-A(I,J)*B(J,L)                                                     
   50 CONTINUE                                                                  
      B(I,L)=SUM                                                                
   55 CONTINUE                                                                  
      RETURN                                                                    
   60 INV=0                                                                     
      DET=0.D00                                                                 
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE RSDIAG(NMAX,N,H,E,Z,IS)                                        
C                                                                               
C CALLED BY SUBR. HNDIAG, NILSON                                                
C                                                                               
C****************************************************************               
C     DIE MATRIX H (INPUT LOWER TRIANGLE) WIRD DIAGONALISIERT   *               
C     E (OUTPUT) SIND DIE EIGENWERTE, Z IST HILFS-ARRAY         *               
C     H (OUTPUT) SIND DIE EIGENVEKTOREN                         *               
C     IS=+1     MIT (INK.) ORDNEN DER EIGENWERTE UND SIGNUM     *               
C        +2     WIE +1 MIT TESTEN                               *               
C        -1     MIT (DEK.) ORDNEN DER EIGENWERTE UND SIGNUM     *               
C        -2     WIE -1 MIT TESTEN                               *               
C         0     OHNE ORDNEN, OHNE TESTEN                        *               
C****************************************************************               
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      DIMENSION H(NMAX,NMAX),E(NMAX),Z(NMAX)                                    
C                                                                               
      DATA ERR,EPS/1.D-32,1.D-14/,ZERO,ONE,HALF/0.D0,1.D0,0.5D0/                
C                                                                               
      IF (N.GT.1) GOTO 1                                                        
      E(1)=H(1,1)                                                               
      H(1,1)=ONE                                                                
      RETURN                                                                    
C                                                                               
    1 TR1=ZERO                                                                  
      DO 3 J=1,N                                                                
      DO 2 I=1,N                                                                
      K=MIN(I,J)                                                                
      L=MAX(I,J)                                                                
      TR1=TR1+H(L,K)*H(L,K)                                                     
    2 CONTINUE                                                                  
    3 CONTINUE                                                                  
      TOL=ERR*TR1                                                               
CCC   HOUSEHOLDER-REDUKTION                                                     
   10 I=N                                                                       
   15 IF (I.LT.2) GOTO 200                                                      
   20 L=I-2                                                                     
      F=H(I,I-1)                                                                
      G=F                                                                       
      Q=ZERO                                                                    
      IF (L.LE.0) GOTO 31                                                       
   32 DO 30 K=1,L                                                               
      Q=Q+H(I,K)*H(I,K)                                                         
   30 CONTINUE                                                                  
   31 S=Q+F*F                                                                   
      IF (S.GT.TOL) GOTO 34                                                     
   33 Q=ZERO                                                                    
      GOTO 100                                                                  
   34 IF (Q.LE.TOL) GOTO 33                                                     
   40 L=L+1                                                                     
      G=SQRT(S)                                                                 
      IF (F.GE.ZERO) G=-G                                                       
      Q=S-F*G                                                                   
      T=ONE/Q                                                                   
      H(I,I-1)=F-G                                                              
      F=ZERO                                                                    
      IF (L.LE.0) GOTO 51                                                       
   52 DO 50 J=1,L                                                               
      H(J,I)=H(I,J)*T                                                           
      S=ZERO                                                                    
      DO 55 K=1,J                                                               
   55 S=S+H(J,K)*H(I,K)                                                         
      J1=J+1                                                                    
      IF (L.LT.J1) GOTO 57                                                      
   58 DO 59 K=J1,L                                                              
      S=S+H(K,J)*H(I,K)                                                         
   59 CONTINUE                                                                  
   57 Z(J)=S*T                                                                  
   50 F=F+S*H(J,I)                                                              
   51 F=HALF*F*T                                                                
C                                                                               
      IF (L.LE.0) GOTO 100                                                      
   62 DO 60 J=1,L                                                               
      S=H(I,J)                                                                  
      Z(J)=Z(J)-F*S                                                             
      P=Z(J)                                                                    
      DO 65 K=1,J                                                               
      H(J,K)=H(J,K)-S*Z(K)-H(I,K)*P                                             
   65 CONTINUE                                                                  
   60 CONTINUE                                                                  
  100 CONTINUE                                                                  
      E(I)=Q                                                                    
      Z(I-1)=G                                                                  
      I=I-1                                                                     
      GOTO 15                                                                   
CCC   BEREITSTELLEN DER TRANSFORMATIONMATRIX                                    
  200 E(1)=ZERO                                                                 
      Z(N)=ZERO                                                                 
      B=ZERO                                                                    
      F=ZERO                                                                    
      DO 210 I=1,N                                                              
      L=I-1                                                                     
      IF (E(I).EQ.ZERO) GOTO 221                                                
      IF (L.LE.0) GOTO 221                                                      
  222 DO 220 J=1,L                                                              
      S=ZERO                                                                    
      DO 225 K=1,L                                                              
      S=S+H(I,K)*H(K,J)                                                         
  225 CONTINUE                                                                  
      DO 226 K=1,L                                                              
      H(K,J)=H(K,J)-S*H(K,I)                                                    
  226 CONTINUE                                                                  
  220 CONTINUE                                                                  
C                                                                               
  221 E(I)=H(I,I)                                                               
      H(I,I)=ONE                                                                
      IF (L.LE.0) GOTO 210                                                      
  232 DO 230 J=1,L                                                              
      H(I,J)=ZERO                                                               
      H(J,I)=ZERO                                                               
  230 CONTINUE                                                                  
  210 CONTINUE                                                                  
CCC   DIAGONALISIEREN DER DREIECKSMATRIX                                        
      DO 300 L=1,N                                                              
      Q=EPS*(ABS(E(L))+ABS(Z(L)))                                               
      IF (Q.GT.B) B=Q                                                           
CCC   TEST FUER SPLITTING                                                       
      DO 310 J=L,N                                                              
      IF (ABS(Z(J)).LE.B) GOTO 320                                              
  310 CONTINUE                                                                  
CCC   TEST FUER KONVERGENZ                                                      
  320 IF (J.EQ.L) GOTO 300                                                      
  340 P=HALF*(E(L+1)-E(L))/Z(L)                                                 
      R=SQRT(P*P+ONE)                                                           
      T=P+R                                                                     
      IF (P.LT.ZERO) T=P-R                                                      
      Q=E(L)-Z(L)/T                                                             
      DO 350 I=L,N                                                              
      E(I)=E(I)-Q                                                               
  350 CONTINUE                                                                  
      F=F+Q                                                                     
CCC   QR-TRANSFORMATION                                                         
      P=E(J)                                                                    
      C=ONE                                                                     
      S=ZERO                                                                    
      I=J                                                                       
  360 I=I-1                                                                     
      IF (I.LT.L) GOTO 362                                                      
      G=C*Z(I)                                                                  
      Q=C*P                                                                     
      IF (ABS(P).LT.ABS(Z(I))) GOTO 363                                         
  364 C=Z(I)/P                                                                  
      R=SQRT(C*C+ONE)                                                           
      Z(I+1)=S*P*R                                                              
      S=C/R                                                                     
      C=ONE/R                                                                   
      GOTO 365                                                                  
  363 C=P/Z(I)                                                                  
      R=SQRT(C*C+ONE)                                                           
      Z(I+1)=S*Z(I)*R                                                           
      S=ONE/R                                                                   
      C=C/R                                                                     
  365 P=C*E(I)-S*G                                                              
      E(I+1)=Q+S*(C*G+S*E(I))                                                   
      DO 368 K=1,N                                                              
      Q=H(K,I+1)                                                                
      H(K,I+1)=H(K,I)*S+Q*C                                                     
      H(K,I)=H(K,I)*C-Q*S                                                       
  368 CONTINUE                                                                  
      GOTO 360                                                                  
  362 Z(L)=S*P                                                                  
      E(L)=C*P                                                                  
      IF (ABS(Z(L)).GT.B) GOTO 340                                              
  300 E(L)=E(L)+F                                                               
CCC   KONVERGENZ                                                                
      IF (IABS(IS).LT.2) GOTO 400                                               
      TR2=ZERO                                                                  
      DO 370 I=1,N                                                              
  370 TR2=TR2+E(I)*E(I)                                                         
      T=ZERO                                                                    
      DO 372 I=1,N                                                              
      DO 372 J=1,N                                                              
      S=ZERO                                                                    
      DO 371 K=1,N                                                              
  371 S=S+H(K,I)*H(K,J)                                                         
  372 T=T+S                                                                     
      ER1=TR2/TR1-ONE                                                           
      ER2=T/N-ONE                                                               
      IF (ABS(ER1).LT.EPS.AND.ABS(ER2).LT.EPS) GOTO 400                         
      PRINT*,'RSDISG TEST FAILED'                                               
      PRINT*,' TRACE CHECK OF H**2:',ER1                                        
      PRINT*,' UNITARITY CHECK:',ER2                                            
      STOP                                                                      
CCC   ORDNEN DER EIGENWERTE                                                     
  400 IF (IS.EQ.0) RETURN                                                       
      DO 450 I=1,N                                                              
      K=I                                                                       
      P=E(I)                                                                    
      J1=I+1                                                                    
      IF (J1.GT.N) GOTO 450                                                     
  401 DO 410 J=J1,N                                                             
      IF (IS.GT.0) THEN                                                         
       IF (E(J).GE.P) GOTO 410                                                  
      ELSE                                                                      
       IF (E(J).LE.P) GOTO 410                                                  
      ENDIF                                                                     
      K=J                                                                       
      P=E(J)                                                                    
  410 CONTINUE                                                                  
  420 IF (K.EQ.I) GOTO 450                                                      
      E(K)=E(I)                                                                 
      E(I)=P                                                                    
      DO 425 J=1,N                                                              
      P=H(J,I)                                                                  
      H(J,I)=H(J,K)                                                             
      H(J,K)=P                                                                  
  425 CONTINUE                                                                  
  450 CONTINUE                                                                  
CCC   SIGNUM (MAKE THE LARGEST AMPLITUDE POSITIVE)                              
      DO 71 K=1,N                                                               
      S=ZERO                                                                    
      DO 72 I=1,N                                                               
      Q=ABS(H(I,K))                                                             
      IF (Q.LT.S) GOTO 72                                                       
      S=Q                                                                       
      L=I                                                                       
   72 CONTINUE                                                                  
      IF (H(L,K).GE.ZERO) GOTO 71                                               
      DO 73 I=1,N                                                               
      H(I,K)=-H(I,K)                                                            
   73 CONTINUE                                                                  
   71 CONTINUE                                                                  
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE HNDIAG(N,M,EPS,IS)                                             
C                                                                               
C CALLED BY SUB. INTEG                                                          
C                                                                               
C*********************************************************************          
C     DIAGONALISIERT  HH * FF = EH * HN * FF                         *          
C       DURCH    HN * UU = EN * UU                                   *          
C         UND    HN**(-1/2) * HH * HN**(-1/2) * VV = EH * VV         *          
C       NB: EIGENWERTE DER NORM EN=<EPS WERDEN ABGESCHNITTEN         *          
C           DIE MATRIZEN HH UND HN BLEIBEN UNVERAENDERT              *          
C                                                                    *          
C     IS = 0   NUR EIGENWERTE EH,EN UND EIGENVEKTOREN VV,UU          *          
C     IS = 1   AUCH EIGENVEKTOREN FF = UU * HN**(-1/2) * VV          *          
C                                                                    *          
C     EH(N) EIGENWERTE DES VERALLGEMEINERTEN EIGENWERTPROBLEMS       *          
C     EN(N) EIGENWERTE DER NORMMATRIX HN                             *          
C     Z(N)  HILFSFELD                                                *          
C*********************************************************************          
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /HND   / HH(NNRANT,NNRANT),HN(NNRANT,NNRANT),                      
     1                UU(NNRANT,NNRANT),VV(NNRANT,NNRANT),                      
     2                FF(NNRANT,NNRANT),EH(NNRANT),EN(NNRANT),                  
     3                Z(NNRANT)                                                 
C                                                                               
      DO 1 I=1,N                                                                
      DO 1 K=1,I                                                                
      HH(K,I)=HH(I,K)                                                           
    1 UU(I,K)=HN(I,K)                                                           
C                                                                               
      CALL RSDIAG(NNRANT,N,UU,EN,Z,-1)                                          
      M=0                                                                       
      IF (EN(N).LT.-EPS) GOTO 4                                                 
      DO 2 I=1,N                                                                
      IF (EN(I).LE.EPS) GOTO 3                                                  
      M=M+1                                                                     
    2 CONTINUE                                                                  
    3 CONTINUE                                                                  
      IF (M.GT.0) GOTO 5                                                        
      WRITE(*,200) EPS                                                          
      WRITE(1,200) EPS                                                          
  200 FORMAT(/' NORM EIGENVALUES TOO SMALL: TOLERANCE',E12.4)                   
      RETURN                                                                    
C                                                                               
    4 CONTINUE                                                                  
      WRITE(*,300) EN(N)                                                        
      WRITE(1,300) EN(N)                                                        
  300 FORMAT(/' NEGATIVE NORM EIGENVALUE (',E12.4,') EXCEEDS TOLERANCE')        
      RETURN                                                                    
C                                                                               
    5 CONTINUE                                                                  
      DO 9 K=1,M                                                                
      DO 7 I=1,N                                                                
      S=0.D0                                                                    
      DO 6 J=1,N                                                                
    6 S=S+HH(I,J)*UU(J,K)                                                       
    7 Z(I)=S                                                                    
      DO 9 L=K,M                                                                
      S=0.D0                                                                    
      DO 8 I=1,N                                                                
    8 S=S+UU(I,L)*Z(I)                                                          
    9 VV(L,K)=S/SQRT(EN(L)*EN(K))                                               
C                                                                               
      CALL RSDIAG(NNRANT,M,VV,EH,Z,1)                                           
      IF (IS.EQ.0) RETURN                                                       
C                                                                               
      DO 11 I=1,N                                                               
      DO 11 L=1,M                                                               
      S=0.D0                                                                    
      DO 10 K=1,M                                                               
      S=S+UU(I,K)*VV(K,L)/SQRT(EN(K))                                           
   10 CONTINUE                                                                  
      FF(I,L)=S                                                                 
   11 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE ROTATE(IT,NSHL)                                                
C                                                                               
C CALLED BY SUB. INTEG                                                          
C                                                                               
C****************************************************************               
C     CALCULATES THE OVERLAPS BETWEEN                           *               
C     NILSSON WAVEFUNCTIONS AND ROTATED NILSSON WAVEFUNCTIONS   *               
C     FOR THE SHELL NSHL AND THE ISOSPIN IT                     *               
C****************************************************************               
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /SDF   / SDF1(NJMAX,NJMAX,NJMAX),SDF2(NJMAX,NJMAX,NJMAX)           
      COMMON /DEF   / QEL(NMD,NMD,5,6),PEL(NMD,NMD,6),                          
     1                SEL(NMD,NMD,3,6),REL(NMD,NMD,3,6),                        
     2                HEL(NMD,NMD,6),AMP(NMD,NMD,6)                             
      COMMON /WORK1 / W1(NMD,NMD),W2(NMD,NMD)                                   
      COMMON /WORK2 / A1(NMD,NMD),A2(NMD,NMD),                                  
     1                B1(NMD,NMD),B2(NMD,NMD),                                  
     2                C1(NMD,NMD),C2(NMD,NMD)                                   
C                                                                               
      MIX = NSHL+(IT-1)*MXSH(1)                                                 
      NS  = NPSH(NSHL,IT)+1                                                     
      NL  = NILS(NSHL,IT)                                                       
C                                                                               
      DO 10 N1=1,NL                                                             
      DO 10 N2=1,NL                                                             
      W1(N2,N1)=0.D0                                                            
      W2(N2,N1)=0.D0                                                            
   10 CONTINUE                                                                  
C                                                                               
      DO 5 J=1,NS                                                               
      JJ=(J*(J-1))/2                                                            
      DO 2 N2=1,NL                                                              
      DO 2 M=1,J                                                                
      A=0.D0                                                                    
      B=0.D0                                                                    
      DO 1 K=1,J                                                                
      L2=JJ+K                                                                   
      C=AMP(L2,N2,MIX)                                                          
      A=A+C*SDF1(K,M,J)                                                         
      B=B+C*SDF2(K,M,J)                                                         
    1 CONTINUE                                                                  
      A1(M,N2)=A                                                                
      A2(M,N2)=B                                                                
    2 CONTINUE                                                                  
      DO 4 N1=1,NL                                                              
      DO 4 N2=1,NL                                                              
      A=0.D0                                                                    
      B=0.D0                                                                    
      DO 3 M=1,J                                                                
      L1=JJ+M                                                                   
      C=AMP(L1,N1,MIX)                                                          
      A=A+C*A1(M,N2)                                                            
      B=B+C*A2(M,N2)                                                            
    3 CONTINUE                                                                  
      W1(N2,N1)=W1(N2,N1)+A                                                     
      W2(N2,N1)=W2(N2,N1)+B                                                     
    4 CONTINUE                                                                  
    5 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE CONTRA(IT,NSHL,OLAP)                                           
C                                                                               
C CALLED BY SUB. INTEG                                                          
C                                                                               
C*************************************************************                  
C       EVALUATION OF THREE CONTRACTIONS AND OVERLAP                            
C*************************************************************                  
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CON   / CA1(NMD,NMD,6),CA2(NMD,NMD,6),CB1(NMD,NMD,6),             
     1                CB2(NMD,NMD,6),CC1(NMD,NMD,6),CC2(NMD,NMD,6)              
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /UVUV  / ESP(NMD,6),EQP(NMD,6),U(NMD,6),V(NMD,6),                  
     1                KQN(NMD,6)                                                
      COMMON /WORK1 / W1(NMD,NMD),W2(NMD,NMD)                                   
      COMMON /WORK2 / A1(NMD,NMD),A2(NMD,NMD),                                  
     1                B1(NMD,NMD),B2(NMD,NMD),                                  
     2                C1(NMD,NMD),C2(NMD,NMD)                                   
C                                                                               
      MIX = NSHL + (IT-1)*MXSH(1)                                               
      NL  = NILS(NSHL,IT)                                                       
C                                                                               
      DO 1 N1=1,NL                                                              
      R=U(N1,MIX)                                                               
      S=V(N1,MIX)                                                               
      DO 2 N2=1,NL                                                              
      SS=U(N2,MIX)                                                              
      TT=V(N2,MIX)                                                              
      P=R*SS+S*TT                                                               
      Q=R*TT-S*SS                                                               
      SS=W1(N1,N2)                                                              
      TT=W2(N1,N2)                                                              
      A1(N1,N2)=P*SS                                                            
      A2(N1,N2)=P*TT                                                            
      B1(N1,N2)=Q*TT                                                            
      B2(N1,N2)=-Q*SS                                                           
    2 C1(N2,N1)=0.D00                                                           
    1 C1(N1,N1)=1.D00                                                           
C                                                                               
      CALL AINVBG(NMD,A1,C1,NL,NL,INV,TT)                                       
      IF (INV.EQ.0) GOTO 8                                                      
      IF (TT.LT.0.D0) TT=-TT                                                    
      DO 4 N2 = 1,NL                                                            
      DO 4 N1 = 1,NL                                                            
      SS = 0.D00                                                                
      DO 5 N3 = 1,NL                                                            
    5 SS = SS + C1(N1,N3) * A2(N3,N2)                                           
    4 A1(N1,N2) = SS                                                            
      DO 46 N2=1,NL                                                             
      DO 45 N1=1,NL                                                             
      C2(N1,N2) = 0.D00                                                         
      SS = 0.D00                                                                
      DO 44 N3 = 1,NL                                                           
   44 SS = SS + A1(N1,N3)*A1(N3,N2)                                             
   45 A2(N1,N2) = SS                                                            
      C2(N2,N2) = 1.D00                                                         
   46 A2(N2,N2) = 1.D00 + A2(N2,N2)                                             
C                                                                               
      CALL AINVBG(NMD,A2,C2,NL,NL,INV,SS)                                       
      IF (INV.EQ.0) GOTO 8                                                      
      IF (SS.LT.0.D0) SS=-SS                                                    
      OLAP=DSQRT(SS)*TT                                                         
      DO 52 N1=1,NL                                                             
      DO 52 N2=1,NL                                                             
      SS=0.D00                                                                  
      DO 51 N3=1,NL                                                             
   51 SS=SS+C2(N1,N3)*C1(N3,N2)                                                 
   52 CC1(N1,N2,MIX)=SS                                                         
      DO 55 N1=1,NL                                                             
      DO 55 N2=1,NL                                                             
      SS=0.D00                                                                  
      DO 54 N3=1,NL                                                             
      TT=CC1(N3,N2,MIX)                                                         
   54 SS=SS-A1(N1,N3)*TT                                                        
   55 CC2(N1,N2,MIX)=SS                                                         
      DO 68 N1=1,NL                                                             
      DO 68 N2=1,NL                                                             
      P=0.D00                                                                   
      Q=0.D00                                                                   
      R=0.D00                                                                   
      S=0.D00                                                                   
      DO 63 N3=1,NL                                                             
      SS=CC1(N3,N2,MIX)                                                         
      TT=CC2(N3,N2,MIX)                                                         
      UU=CC1(N1,N3,MIX)                                                         
      VV=CC2(N1,N3,MIX)                                                         
      P=P+B1(N1,N3)*SS-B2(N1,N3)*TT                                             
      Q=Q+UU*B1(N3,N2)-VV*B2(N3,N2)                                             
      R=R+B1(N1,N3)*TT+B2(N1,N3)*SS                                             
      S=S+UU*B2(N3,N2)+VV*B1(N3,N2)                                             
   63 CONTINUE                                                                  
      CA1(N1,N2,MIX)=P                                                          
      CB1(N1,N2,MIX)=Q                                                          
      CA2(N1,N2,MIX)=R                                                          
      CB2(N1,N2,MIX)=S                                                          
   68 CONTINUE                                                                  
      RETURN                                                                    
C                                                                               
    8 WRITE(*,100) MIX, NL                                                      
  100 FORMAT (' STOP: TRANSFORMATION MATRIX IN CONTRA SINGULAR',                
     1        ' AT MIX WITH NL =', I2,I3)                                       
      STOP                                                                      
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE SPHSET                                                         
C                                                                               
C CALLED BY MAIN                                                                
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (NIE=127)                                                       
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      DIMENSION MQN(12),MQ(2)                                                   
      DIMENSION CC(NJMAX),CG(NJMAX,NJMAX,2)                                     
      DIMENSION CD(NJMAX),CH(NJMAX,NJMAX,2)                                     
C                                                                               
      COMMON /GFVSQ / SQ(0:NIE)                                                 
      COMMON /MAXMAX/ EMAX(3),KMAX(3),NMAX(3),JMAX                              
      COMMON /SPHJLM/ J1(NMD),L1(NMD),M1(NMD)                                   
      COMMON /SPH   / ES(NMD),Q(NMD,NMD,5),R(NMD,NMD,3),P(NMD,NMD),             
     1                SP(NMD,NMD,3)                                             
C                                                                               
      DATA MQN/1,-1,3,-3,5,-5,7,-7,9,-9,11,-11/                                 
      DATA ZERO,HALF,TWO,FOUR/0.0,0.5,2.0,4.0/                                  
C                                                                               
      DO 5 I=1,NMD                                                              
      DO 5 J=1,NMD                                                              
      P(I,J)=ZERO                                                               
    5 CONTINUE                                                                  
      DO 6 K=1,3                                                                
      DO 6 I=1,NMD                                                              
      DO 6 J=1,NMD                                                              
      R(I,J,K)=ZERO                                                             
      SP(I,J,K)=ZERO                                                            
    6 CONTINUE                                                                  
      DO 7 K=1,5                                                                
      DO 7 I=1,NMD                                                              
      DO 7 J=1,NMD                                                              
      Q(I,J,K)=ZERO                                                             
    7 CONTINUE                                                                  
C                                                                               
      K=0                                                                       
      DO 4 JS=1,JMAX                                                            
      SJ=JS                                                                     
      SJ=SJ-HALF                                                                
      JMIN=MAX0(JS-4,1)                                                         
      S=1                                                                       
      IF (MOD(JS+JMIN,2).EQ.0) S=-1                                             
      DO 1 JR=JMIN,JS                                                           
      RJ=JR                                                                     
      RJ=RJ-HALF                                                                
      S=-S                                                                      
      CD(JR)=CLEBG(SJ,HALF,FOUR,ZERO,RJ,HALF)*S                                 
    1 CC(JR)=CLEBG(SJ,HALF, TWO,ZERO,RJ,HALF)*S                                 
C                                                                               
      DO 4 MJS=1,JS                                                             
      MS=MQN(MJS)                                                               
      SM=MS                                                                     
      SM=SM-HALF                                                                
      K=K+1                                                                     
      J1(K)=JS                                                                  
      M1(K)=MS                                                                  
C     WRITE(1,100) K,J1(K),M1(K)                                                
  100 FORMAT(I5,'    J1 =',I3,'   M1 =',I3)                                     
C                                                                               
      DO 20 JR=JMIN,JS                                                          
      RJ=JR                                                                     
      RJ=RJ-HALF                                                                
      S=1                                                                       
      IF (MOD(JR,2).EQ.0) S=-1                                                  
      DO 20 MJR=1,JR                                                            
      MR=MQN(MJR)                                                               
      RM=MR                                                                     
      RM=RM-HALF                                                                
      CH(JR,MJR,1)=CLEBG(RJ,-RM,FOUR,SM+RM,SJ,SM)*S                             
      CH(JR,MJR,2)=CLEBG(RJ,+RM,FOUR,SM-RM,SJ,SM)                               
      CG(JR,MJR,1)=CLEBG(RJ,-RM, TWO,SM+RM,SJ,SM)*S                             
   20 CG(JR,MJR,2)=CLEBG(RJ,+RM, TWO,SM-RM,SJ,SM)                               
C                                                                               
      DO 4 MV=1,2                                                               
      DO 3 L=1,K                                                                
      JJ=J1(K)-J1(L)                                                            
      IF (JJ.GT.4) GOTO 3                                                       
      JR=J1(L)                                                                  
      MR=M1(L)                                                                  
      MQ(1)=MS+MR-1                                                             
      MQ(2)=MS-MR                                                               
      MJR=MR                                                                    
      IF (MR.LT.0) MJR=1-MR                                                     
C     Y4-0                                                                      
      MU=MQ(MV)+5                                                               
      IF (MU.GT.9.OR.MU.LT.1) GOTO 3                                            
      IF (MU.EQ.5) P(K,L)=CD(JR)*CH(JR,MJR,MV)                                  
      IF (JJ.GT.2) GOTO 3                                                       
C     Y2-MU                                                                     
      MU=MQ(MV)+3                                                               
      IF (MU.GT.5.OR.MU.LT.1) GOTO 3                                            
      Q(K,L,MU)=CC(JR)*CG(JR,MJR,MV)                                            
      IF (MOD(MU,2).EQ.0) Q(L,K,MU)=-Q(K,L,MU)                                  
      IF (JJ.NE.0) GOTO 3                                                       
C     J-MU                                                                      
      MU=MQ(MV)+2                                                               
      IF (MU.GT.3.OR.MU.LT.1) GOTO 3                                            
      IF (MU.EQ.2) GOTO 2                                                       
      MT=(1+MQ(MV))/2                                                           
      MT=MQ(MV)*MR-MT                                                           
      S=SQ(JR+MT)*SQ(JR-MT)                                                     
      IF (MOD(JR,2).EQ.0) S=-S                                                  
      R(K,L,MU)=S                                                               
      R(L,K,MU)=S                                                               
      GOTO 3                                                                    
    2 RM=MR                                                                     
      RM=RM-HALF                                                                
      R(K,L,MU)=RM                                                              
      R(L,K,MU)=RM                                                              
    3 CONTINUE                                                                  
C                                                                               
C     END OF THE JS LOOP                                                        
C                                                                               
    4 CONTINUE                                                                  
C                                                                               
      DO 10 MU=1,5,2                                                            
      MUU=6-MU                                                                  
      DO 10 I=1,K                                                               
      DO 10 L=1,I                                                               
   10 Q(L,I,MU)=Q(I,L,MUU)                                                      
C     DO 11 MU=1,9,2                                                            
      DO 11 I=1,K                                                               
      DO 11 L=1,I                                                               
   11 P(L,I)=P(I,L)                                                             
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE DEFSET(IT,N)                                                   
C                                                                               
C CALLED BY SUB. NILBCS                                                         
C                                                                               
C****************************************************************               
C                                                               *               
C   COMMONBLOCK SPH/DEF: MATRICES IN SPHERICAL/DEFORMED BASIS   *               
C   ---------------------------------------------------------   *               
C     IT : 1 FOR NEUTRONS AND 2 FOR PROTONS                     *               
C     MIX: NUMBERING OF MAJOR SHELLS (THROUGHOUT IT=1 AND 2)    *               
C       SPHERICAL BASIS: (FOR A GIVEN SHELL N AND IT)           *               
C     ES(N1)         SINGLE-PARTICLE ENERGY                     *               
C     RSQ(N1,N2)     MATRIX ELEMENTS <N1|R**2|N2>               *               
C     Q(N1,N2,MU)    Y2-MU MATRIXELEMENTS IN SPHERICAL BASIS    *               
C     P(N1,N2)       Y4-0  MATRIXELEMENTS IN SPHERICAL BASIS    *               
C     R(N1,N2,MU)    J-MU                                       *               
C     SP(N1,M2,MU)   SPIN                                       *               
C       DEFORMED BASIS: (THROUGHOUT SHELLS AND IT=1,2)          *               
C     ESP(N1,MIX)       NILSSON ENERGY                          *               
C     AMP(N1,N2,MIX)    NILSSON WAVEFUNCTION                    *               
C     HEL(N1,N2,MIX)    SINGLE-PARTICLE HAMILTONIAN             *               
C     QEL(N1,N2,MU,MIX) R**2 * Y2-MU                            *               
C     PEL(N1,N2,MIX)    R**2 * Y4-0                             *               
C     REL(N1,N2,MU,MIX) ANGULAR MOMENTUM J-MU                   *               
C     SEL(N1,N2,MU,MIX) SPIN MATRIX ELEMENTS                    *               
C                                                               *               
C****************************************************************               
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (NJMAX1=NJMAX+1)                                                
      PARAMETER (NIE=127)                                                       
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /GFVSQ / SQ(0:NIE)                                                 
      COMMON /MAXMAX/ EMAX(3),KMAX(3),NMAX(3),JMAX                              
      COMMON /SPHJLM/ J1(NMD),L1(NMD),M1(NMD)                                   
      COMMON /UVUV  / ESP(NMD,6),EQP(NMD,6),U(NMD,6),V(NMD,6),                  
     1                KQN(NMD,6)                                                
      COMMON /NIL   / VAL(NMD),VEC(NMD,NMD),RSQ(NJMAX1,NJMAX1),                 
     1                W0(2),EP1,EP2,NS,ND                                       
      COMMON /SPH   / ES(NMD),Q(NMD,NMD,5),R(NMD,NMD,3),P(NMD,NMD),             
     1                SP(NMD,NMD,3)                                             
      COMMON /DEF   / QEL(NMD,NMD,5,6),PEL(NMD,NMD,6),                          
     1                SEL(NMD,NMD,3,6),REL(NMD,NMD,3,6),                        
     2                HEL(NMD,NMD,6),AMP(NMD,NMD,6)                             
C                                                                               
      EPS1=EP1*W0(IT)                                                           
      EPS2=EP2*W0(IT)                                                           
C                                                                               
C     STORE NILSSON ENERGY AND WAVEFUNCTION                                     
C                                                                               
      NL=ND                                                                     
      MIX=N+(IT-1)*MXSH(1)                                                      
      DO 1 I=1,NL                                                               
      ESP(I,MIX)=VAL(I)                                                         
      DO 1 K=1,ND                                                               
    1 AMP(K,I,MIX)=VEC(K,I)                                                     
C                                                                               
C     SPHERICAL SINGLE-PARTICLE HAMILTONIAN MATRIX ELEMENTS                     
C                                                                               
      DO 88 N1=1,NL                                                             
      DO 88 N2=1,N1                                                             
      SS=0.                                                                     
      DO 87 N3=1,NL                                                             
   87 SS=SS+VEC(N3,N1)*ES(N3)*VEC(N3,N2)                                        
      HEL(N1,N2,MIX)=SS                                                         
   88 HEL(N2,N1,MIX)=SS                                                         
C                                                                               
C     QUADRUPOLE MATRIX ELEMENTS                                                
C                                                                               
      DO 50 MU=1,5                                                              
      DO 2 N1=1,ND                                                              
      DO 2 N2=1,ND                                                              
    2 VEC(N2,N1)=Q(N2,N1,MU)*RSQ(L1(N2),L1(N1))                                 
      DO 40 N1=1,NL                                                             
      DO 19 N2=1,ND                                                             
      SS=0.                                                                     
      DO 18 N3=1,ND                                                             
   18 SS=SS+AMP(N3,N1,MIX)*VEC(N3,N2)                                           
      VAL(N2)=SS                                                                
   19 CONTINUE                                                                  
      DO 40 N2=1,ND                                                             
      SS=0.                                                                     
      DO 25 N3=1,ND                                                             
   25 SS=SS+VAL(N3)*AMP(N3,N2,MIX)                                              
   40 QEL(N1,N2,MU,MIX) = SS                                                    
   50 CONTINUE                                                                  
C                                                                               
C     HEXADECUPOLE MATRIX ELEMENTS                                              
C                                                                               
      DO 8 N1=1,ND                                                              
      DO 8 N2=1,ND                                                              
    8 VEC(N2,N1)=P(N2,N1)*RSQ(L1(N2),L1(N1))                                    
      DO 12 N1=1,NL                                                             
      DO 10 N2=1,ND                                                             
      SS=0.                                                                     
      DO 9 N3=1,ND                                                              
      SS=SS+AMP(N3,N1,MIX)*VEC(N3,N2)                                           
    9 CONTINUE                                                                  
      VAL(N2)=SS                                                                
   10 CONTINUE                                                                  
      DO 12 N2=1,ND                                                             
      SS=0.                                                                     
      DO 11 N3=1,ND                                                             
   11 SS=SS+VAL(N3)*AMP(N3,N2,MIX)                                              
   12 PEL(N1,N2,MIX) = SS                                                       
C                                                                               
C     ANGULAR MOMENTUM MATRIX ELEMENTS                                          
C                                                                               
      DO 70 MU=1,3                                                              
      DO 60 N1=1,NL                                                             
      DO 54 N2=1,ND                                                             
      SS=0.                                                                     
      DO 53 N3=1,ND                                                             
   53 SS=SS+AMP(N3,N1,MIX)*R(N3,N2,MU)                                          
      VAL(N2)=SS                                                                
   54 CONTINUE                                                                  
      DO 60 N2=1,N1                                                             
      SS=0.                                                                     
      DO 56 N3=1,ND                                                             
   56 SS=SS+VAL(N3)*AMP(N3,N2,MIX)                                              
      REL(N1,N2,MU,MIX) = SS                                                    
   60 REL(N2,N1,MU,MIX) = SS                                                    
   70 CONTINUE                                                                  
C                                                                               
C     K-QUANTUM NUMBER (2*J_0)                                                  
C                                                                               
      DO 3 I=1,NL                                                               
    3 KQN(I,MIX)=NINT(2.*REL(I,I,2,MIX))                                        
C                                                                               
C     SPIN MATRIX ELEMENTS                                                      
C                                                                               
      DO 200 N1=1,ND                                                            
      KC1=(J1(N1)-L1(N1)+1)*2-1                                                 
      SL=L1(N1)-0.5                                                             
      DL=2.*SL                                                                  
      AM1=M1(N1)-0.5                                                            
      DO 180 N2=1,N1                                                            
      KC2=(J1(N2)-L1(N2)+1)*2-1                                                 
      AM2=M1(N2)-0.5                                                            
      DO 110 MU=1,3                                                             
      SP(N2,N1,MU)=0.                                                           
      SP(N1,N2,MU)=0.                                                           
  110 CONTINUE                                                                  
      IF (L1(N1).NE.L1(N2)) GO TO 180                                           
      STM=1.                                                                    
      IF (M1(N1).EQ.M1(N2)) GO TO 120                                           
      AM2=-AM2                                                                  
      ITM=J1(N2)-M1(N2)                                                         
      IF (ITM.NE.(ITM/2)*2) STM=-1.                                             
      IF (M1(N1)+M1(N2)-1.EQ.1) GO TO 140                                       
      IF (M1(N1)+M1(N2)-1.EQ.-1) GO TO 160                                      
      GOTO 180                                                                  
  120 MU=2                                                                      
      IL=NINT(SL+AM1)                                                           
      IR=NINT(SL-AM1)                                                           
      IF (KC1.EQ.KC2) SPM=KC1*AM1/DL                                            
      IF (KC1.NE.KC2) SPM=SQ(IR)*SQ(IL)/DL                                      
      GOTO 170                                                                  
  140 MU=3                                                                      
      IL=NINT(SL+KC1*AM1)                                                       
      IR=NINT(SL-KC2*AM2)                                                       
      IF (KC1.EQ.KC2) SPM=KC2*SQ(IL)*SQ(IR)/DL                                  
      IF (KC1.NE.KC2) SPM=KC2*SQ(IL)*SQ(IR)/DL                                  
      GOTO 170                                                                  
  160 MU=1                                                                      
      IL=NINT(SL-KC1*AM1)                                                       
      IR=NINT(SL+KC2*AM2)                                                       
      IF (KC1.EQ.KC2) SPM=KC1*SQ(IL)*SQ(IR)/DL                                  
      IF (KC1.NE.KC2) SPM=KC1*SQ(IL)*SQ(IR)/DL                                  
  170 SP(N2,N1,MU)=SPM*STM                                                      
      SP(N1,N2,MU)=SP(N2,N1,MU)                                                 
      IF (M1(N1).NE.1.OR.M1(N2).NE.1) GO TO 180                                 
      ITM=J1(N2)-M1(N2)                                                         
      IF (ITM.NE.(ITM/2)*2) STM=-1.                                             
      IL=NINT(SL+KC1*AM1)                                                       
      IR=NINT(SL+KC2*AM1)                                                       
      SPM=KC2*SQ(IL)*SQ(IR)/DL                                                  
      SP(N2,N1,3)=SPM*STM                                                       
      SP(N1,N2,3)=SPM*STM                                                       
  180 CONTINUE                                                                  
  200 CONTINUE                                                                  
C                                                                               
      DO 270 MU=1,3                                                             
      DO 260 N1=1,NL                                                            
      DO 254 N2=1,ND                                                            
      SS=0.                                                                     
      DO 253 N3=1,ND                                                            
  253 SS=SS+AMP(N3,N1,MIX)*SP(N3,N2,MU)                                         
      VAL(N2)=SS                                                                
  254 CONTINUE                                                                  
      DO 258 N2=1,N1                                                            
      SS=0.                                                                     
      DO 256 N3=1,ND                                                            
  256 SS=SS+VAL(N3)*AMP(N3,N2,MIX)                                              
      SEL(N1,N2,MU,MIX)=SS                                                      
  258 CONTINUE                                                                  
  260 CONTINUE                                                                  
  270 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE SDFSET(JMAX,COSBET)                                            
C                                                                               
C CALLED BY SUB. INTEG                                                          
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      DIMENSION MQN(12)                                                         
C                                                                               
      COMMON /SDF   / SDF1(NJMAX,NJMAX,NJMAX),SDF2(NJMAX,NJMAX,NJMAX)           
C                                                                               
      DATA MQN/1,-1,3,-3,5,-5,7,-7,9,-9,11,-11/                                 
C                                                                               
      SN=-1                                                                     
      DO 1 J=1,JMAX                                                             
      SN=-SN                                                                    
      SJ=J-0.5                                                                  
      DO 1 M=1,J                                                                
      SM=MQN(M)-0.5                                                             
      DO 1 K=1,J                                                                
      SK=MQN(K)-0.5                                                             
      SDF1(K,M,J)=SDFUNC(SJ,SM,SK,COSBET)                                       
    1 SDF2(K,M,J)=SN*SDFUNC(SJ,-SM,SK,COSBET)                                   
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE GFV(IE)                                                        
C                                                                               
C CALLED BY MAIN                                                                
C                                                                               
      PARAMETER (NIE=127)                                                       
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /GFVIV / IV(-NIE:NIE)                                              
      COMMON /GFVSQ / SQ(0:NIE)                                                 
      COMMON /GFVLOG/ G(NIE),LOGTAB                                             
C                                                                               
      IF (IE.GT.NIE) THEN                                                       
       WRITE(*,*) ' SET THE PARAMETER NIE IN GFV, NILSON, SPHSET,',             
     +            ' DEFSET, INTEG, CLEBG AND SDFUNC TO',IE                      
       STOP                                                                     
      ENDIF                                                                     
      LOGTAB=IE                                                                 
      IV(0)=1                                                                   
      SQ(0)=0.                                                                  
      DO 10 I=1,IE                                                              
      IV(I)=-IV(I-1)                                                            
      IV(-I)=IV(I)                                                              
      X=I                                                                       
      SQ(I)=SQRT(X)                                                             
   10 CONTINUE                                                                  
      G(1)=0                                                                    
      G(2)=0                                                                    
      DO 20 I=3,IE                                                              
      X=I-1                                                                     
   20 G(I)=G(I-1)+DLOG(X)                                                       
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE LEGEND                                                         
C                                                                               
C CALLED BY MAIN                                                                
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /GAUSSL/ X(NNDEG),W(NNDEG),WEIT(NNDEG),N                           
      DIMENSION GAMMA(NNDEG)                                                    
C                                                                               
      H=0                                                                       
      DO 1 J=1,N                                                                
      H=H+2.                                                                    
      H1=J-1                                                                    
    1 GAMMA(J)=-H1*H1/(H-1.)/(H-3.)                                             
      NP1=N+1                                                                   
      NM1=N-1                                                                   
      Q=1                                                                       
      I=N/2                                                                     
      TEST=N/2.-I-0.1                                                           
      IF (TEST) 2,2,3                                                           
    2 J=I+1                                                                     
      GO TO 4                                                                   
    3 X(I+1)=0                                                                  
      J=I+2                                                                     
    4 H=0                                                                       
      DO 14 K=J,N                                                               
      DOWN=H                                                                    
      UP=1                                                                      
    5 M=0                                                                       
      EWN=(UP+DOWN)/2.                                                          
      DO 9 II=1,N                                                               
      IF (Q) 7,6,7                                                              
    6 Q=1.D-15                                                                  
    7 Q=GAMMA(II)/Q-EWN                                                         
      IF (Q) 8,9,9                                                              
    8 M=M+1                                                                     
    9 CONTINUE                                                                  
      IF (K-M) 10,10,11                                                         
   10 UP=EWN                                                                    
      GO TO 12                                                                  
   11 DOWN=EWN                                                                  
   12 Q=1.D-12-(UP-DOWN)                                                        
      IF (Q) 5,13,13                                                            
   13 H=(UP+DOWN)/2.                                                            
      X(K)=-H                                                                   
   14 X(NP1-K)=H                                                                
      J=I+1                                                                     
      DO 16 K=J,N                                                               
      SUM=1                                                                     
      UP=1                                                                      
      DOWN=0                                                                    
      H1=X(K)                                                                   
      DO 15 II=1,NM1                                                            
      H=2*II                                                                    
      Q=((H-1.)*H1*UP-(II-1)*DOWN)/II                                           
      SUM=SUM+Q*Q*(H+1.)                                                        
      DOWN=UP                                                                   
   15 UP=Q                                                                      
      W(K)=2./SUM                                                               
   16 W(NP1-K)=W(K)                                                             
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE LEGPOL(MAXDEG,COSBET)                                          
C                                                                               
C CALLED BY SUB. NILSON                                                         
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /LEG   / PLEG(NNDEG)                                               
C                                                                               
      PLEG(1)=1.D0                                                              
      PLEG(2)=COSBET                                                            
      M = MAXDEG - 1                                                            
      DO 1 L=2,M                                                                
      A=2*L-1                                                                   
      B=L-1                                                                     
      C=L                                                                       
    1 PLEG(L+1)=(A*COSBET*PLEG(L)-B*PLEG(L-1))/C                                
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE SPUR(JBET,IT,NSHL,NPQP,MU,IS)                                  
C                                                                               
C CALLED BY SUB. INTEG                                                          
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /CON   / CA1(NMD,NMD,6),CA2(NMD,NMD,6),CB1(NMD,NMD,6),             
     1                CB2(NMD,NMD,6),CC1(NMD,NMD,6),CC2(NMD,NMD,6)              
      COMMON /TPQAS / T11PQ(MD2,MD2,6),T20PQ(MD2,MD2,6),                        
     1                T02PQ(MD2,MD2,6)                                          
      COMMON /TQQAS / T11QQ(MD2,MD2,-2:2,6),T20QQ(MD2,MD2,-2:2,6),              
     1                T02QQ(MD2,MD2,-2:2,6)                                     
      COMMON /TPMUQA/ T11PMU(MD2,MD2,-2:2,6),T20PMU(MD2,MD2,-2:2,6),            
     1                T02PMU(MD2,MD2,-2:2,6)                                    
      COMMON /THQAS/ T11HQ(MD2,MD2,6),T20HQ(MD2,MD2,6),T02HQ(MD2,MD2,6)         
      COMMON /TRACE1/ ST02B(2,12),ST20B(2,12),                                  
     1                AAT1(NNRANK,12,2),AAT2(NNRANK,12,2),                      
     2                T1AA(NNRANK,12,2),T2AA(NNRANK,12,2)                       
C                                                                               
      DIMENSION T20(MD2,MD2),T02(MD2,MD2),CC(MD2,MD2)                           
C                                                                               
      MIX = NSHL+(IT-1)*MXSH(1)                                                 
      NL  = NILS(NSHL,IT)                                                       
      NL2 = NL * 2                                                              
      IF (JBET.EQ.1) THEN                                                       
        CALL LARGE(MD2,NMD,NL,CB1(1,1,MIX),CB2(1,1,MIX),CC,+1)                  
      ELSE                                                                      
        CALL LARGE(MD2,NMD,NL,CB1(1,1,MIX),CB2(1,1,MIX),CC, 0)                  
      ENDIF                                                                     
C                                                                               
C     CHOOSE THE OPERATORS ACCORDING TO THE NUMBER II                           
C                                                                               
      IF(NPQP.EQ.0) THEN                                                        
         DO 101 I=1,NL2                                                         
         DO 101 J=1,NL2                                                         
         T02(I,J) = T02PQ(I,J,MIX)                                              
         T20(I,J) = T20PQ(I,J,MIX)                                              
  101    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      IF(NPQP.EQ.1) THEN                                                        
         DO 102 I=1,NL2                                                         
         DO 102 J=1,NL2                                                         
         T02(I,J) = T02QQ(I,J,MU,MIX)                                           
         T20(I,J) = T20QQ(I,J,MU,MIX)                                           
  102    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      IF(NPQP.EQ.2) THEN                                                        
         DO 103 I=1,NL2                                                         
         DO 103 J=1,NL2                                                         
         T02(I,J) = T02PMU(I,J,MU,MIX)                                          
         T20(I,J) = T20PMU(I,J,MU,MIX)                                          
  103    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      IF(NPQP.EQ.3) THEN                                                        
         DO 104 I=1,NL2                                                         
         DO 104 J=1,NL2                                                         
         T02(I,J) = T02HQ(I,J,MIX)                                              
         T20(I,J) = T20HQ(I,J,MIX)                                              
  104    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      SUM1=0.                                                                   
      SUM2=0.                                                                   
      DO 1 I=1,NL2                                                              
      DO 1 J=1,NL2                                                              
      SUM1=SUM1+T20(I,J)*CC(J,I)                                                
      SUM2=SUM2+T02(I,J)*CC(J,I)                                                
    1 CONTINUE                                                                  
      ST20B(IT,IS)=ST20B(IT,IS)+0.5*SUM1                                        
      ST02B(IT,IS)=ST02B(IT,IS)+0.5*SUM2                                        
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE LARGE(MD2,MD,NL,B1,B2,B,K123)                                  
C                                                                               
C CALLED BY SUBS. SPUR; ELEM00                                                  
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      DIMENSION B(MD2,MD2),B1(MD,MD),B2(MD,MD)                                  
C                                                                               
      IF (K123) 5,6,7                                                           
    5 DO 3 I=1,NL                                                               
      DO 3 J=1,NL                                                               
      B(I,J)      =-B2(I,J)                                                     
      B(I,J+NL)   =-B1(I,J)                                                     
      B(I+NL,J)   = B1(I,J)                                                     
      B(I+NL,J+NL)=-B2(I,J)                                                     
    3 CONTINUE                                                                  
      GOTO 10                                                                   
C                                                                               
    6 DO 2 I=1,NL                                                               
      DO 2 J=1,NL                                                               
      B(I,J)      =-B1(I,J)                                                     
      B(I,J+NL)   = B2(I,J)                                                     
      B(I+NL,J)   =-B2(I,J)                                                     
      B(I+NL,J+NL)=-B1(I,J)                                                     
    2 CONTINUE                                                                  
      GOTO 10                                                                   
C                                                                               
    7 DO 1 I=1,NL                                                               
      DO 1 J=1,NL                                                               
      B(I,J)      = B1(I,J)                                                     
      B(I,J+NL)   = B2(I,J)                                                     
      B(I+NL,J)   =-B2(I,J)                                                     
      B(I+NL,J+NL)= B1(I,J)                                                     
    1 CONTINUE                                                                  
   10 CONTINUE                                                                  
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE INTERA                                                         
C                                                                               
C CALLED BY MAIN                                                                
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (NJMAX1=NJMAX+1)                                                
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /BCSBCS/ ANP(2),ALAMB(2),DELTA(2),EP(NNBCS,2),                     
     1                HSP(2),QQ0(2),HH0(2),QP0(2),LNP(2)                        
      COMMON /NIL   / VAL(NMD),VEC(NMD,NMD),RSQ(NJMAX1,NJMAX1),                 
     1                W0(2),EP1,EP2,NS,ND                                       
      COMMON /PARAM / AG(2),AGQ(2),DELE,HOLNG2(2),QCL,HFB,                      
     1                AX(3),AXH(3),XQQ(2),GQP(2)                                
C                                                                               
      QCL = DELE*QQ0(1)*HOLNG2(1) + (1.+DELE)*QQ0(2)*HOLNG2(2)                  
C                                                                               
      QUAD=EP1/(QQ0(1)*W0(1)+QQ0(2)*W0(2))                                      
      AX(1)=QUAD*W0(1)*W0(1)                                                    
      AX(2)=QUAD*W0(2)*W0(2)                                                    
      AX(3)=QUAD*W0(1)*W0(2)                                                    
      XQQ(1)=AX(1)*QQ0(1)+AX(3)*QQ0(2)                                          
      XQQ(2)=AX(2)*QQ0(2)+AX(3)*QQ0(1)                                          
      GQP(1)=AGQ(1)*QP0(1)                                                      
      GQP(2)=AGQ(2)*QP0(2)                                                      
C                                                                               
      HEXA=EP2/(HH0(1)*W0(1)+HH0(2)*W0(2))                                      
      AXH(1)=HEXA*W0(1)*W0(1)                                                   
      AXH(2)=HEXA*W0(2)*W0(2)                                                   
      AXH(3)=HEXA*W0(1)*W0(2)                                                   
C                                                                               
      AAA=SQRT(AX(1))*QQ0(1)+SQRT(AX(2))*QQ0(2)                                 
      AAA=0.5*AAA*AAA                                                           
      HFB=HSP(1)+HSP(2)-AAA-DELTA(1)**2/AG(1)-DELTA(2)**2/AG(2)                 
     1   -ALAMB(1)*ANP(1)-ALAMB(2)*ANP(2)-GQP(1)*QP0(1)-GQP(2)*QP0(2)           
C                                                                               
      WRITE(1,101) QUAD,AX                                                      
      WRITE(1,102) AG,AGQ                                                       
      WRITE(1,103) QQ0,QP0                                                      
      WRITE(1,104) XQQ,GQP                                                      
      WRITE(1,105) W0,(SQRT(HOLNG2(I)),I=1,2)                                   
      WRITE(1,106) HFB,QCL,DELE                                                 
  101 FORMAT(/2X,'QUAD =', E14.6,3X,' X_Q =',3E14.6)                            
  102 FORMAT( 2X,' G_M =',2E14.6,3X,' G_Q =',2E14.6)                            
  103 FORMAT( 2X,' Q_0 =',2E14.6,3X,' P_0 =',2E14.6)                            
  104 FORMAT( 2X,'XQ*Q =',2E14.6,3X,'GQ*P =',2E14.6)                            
  105 FORMAT( 2X,'W_HO =',2E14.6,3X,'B_HO =',2E14.6)                            
  106 FORMAT( 2X,'EHFB =', E14.6,3X,'Q_CL =', E14.6,3X,'DELE =',F5.2)           
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE QPBASE                                                         
C                                                                               
C CALLED BY MAIN                                                                
C                                                                               
C*************************************************************                  
C       TRANSFORMATION TO THE QUASIPARTICLE BASIS            *                  
C             NPQP = 0   FOR  P OPERATOR                     *                  
C             NPQP = 1   FOR  Q OPERATOR                     *                  
C             NPQP = 2   FOR PMU OPERATOR                    *                  
C             NPQP = 3   FOR  H  OPERATOR                    *                  
C               MU = -2,  -1,  0,  +1,  +2                   *                  
C*************************************************************                  
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
C                                                                               
      DO 12 IT = 1,2                                                            
      MISH=MXSH(IT)                                                             
      IF (MISH.LT.1) GOTO 12                                                    
      DO 11 NSHL=1,MISH                                                         
C     OPERATOR P                                                                
      CALL PQASI(IT,NSHL,0)                                                     
C     OPERATOR Q                                                                
      CALL QMQASI(IT,NSHL,1)                                                    
C     OPERATOR PMU                                                              
      CALL PMQASI(IT,NSHL,2)                                                    
C     OPERATOR H                                                                
      CALL HQASI(IT,NSHL,3)                                                     
   11 CONTINUE                                                                  
   12 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE HQASI(IT,NSHL,NPQP)                                            
C                                                                               
C CALLED BY SUB. QPBASE                                                         
C                                                                               
C************************************************************                   
C             NPQP = 0   FOR  P OPERATOR                                        
C             NPQP = 1   FOR  Q OPERATOR                                        
C             NPQP = 2   FOR PMU OPERATOR                                       
C             NPQP = 3   FOR  H  OPERATOR                                       
C               MU = -2,  -1,  0,  +1,  +2                                      
C*************************************************************                  
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /BCSBCS/ ANP(2),ALAMB(2),DELTA(2),EP(NNBCS,2),                     
     1                HSP(2),QQ0(2),HH0(2),QP0(2),LNP(2)                        
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /THQAS / T11HQ(MD2,MD2,6),T20HQ(MD2,MD2,6),                        
     1                T02HQ(MD2,MD2,6)                                          
      COMMON /UVUV  / ESP(NMD,6),EQP(NMD,6),U(NMD,6),V(NMD,6),                  
     1                KQN(NMD,6)                                                
      COMMON /WORK1 / G1(NMD,NMD),G2(NMD,NMD)                                   
      COMMON /DEF   / QEL(NMD,NMD,5,6),PEL(NMD,NMD,6),                          
     1                SEL(NMD,NMD,3,6),REL(NMD,NMD,3,6),                        
     2                HEL(NMD,NMD,6),AMP(NMD,NMD,6)                             
C                                                                               
      IF (NPQP.NE.3) RETURN                                                     
      MU = 0                                                                    
      MIX = NSHL + (IT-1) * MXSH(1)                                             
      NL  = NILS(NSHL,IT)                                                       
      NL2 = NL*2                                                                
C                                                                               
C===================================================================            
C         H  MATRIX  IN  QUASITEILCHENZUSTAND                                   
C===================================================================            
C                                                                               
      DO 111 I=1,NL                                                             
      DO 113 J=1,NL                                                             
      G1(I,J)=U(I,MIX)*HEL(I,J,MIX)                                             
  113 G2(I,J)=V(I,MIX)*HEL(I,J,MIX)                                             
      G1(I,I)=U(I,MIX)*(HEL(I,I,MIX)-ALAMB(IT))                                 
  111 G2(I,I)=V(I,MIX)*(HEL(I,I,MIX)-ALAMB(IT))                                 
      DO 112 I=1,NL                                                             
      DO 112 J=1,NL                                                             
      T11HQ(I,J,MIX)=G1(I,J)*U(J,MIX)-G2(I,J)*V(J,MIX)                          
      T11HQ(I+NL,J+NL,MIX)=T11HQ(I,J,MIX)                                       
      T20HQ(I,J+NL,MIX)=-G1(I,J)*V(J,MIX)-G2(I,J)*U(J,MIX)                      
      T20HQ(I+NL,J,MIX)=-T20HQ(I,J+NL,MIX)                                      
      T02HQ(I,J+NL,MIX)=T20HQ(I,J+NL,MIX)                                       
  112 T02HQ(I+NL,J,MIX)=-T02HQ(I,J+NL,MIX)                                      
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE QMQASI(IT,NSHL,NPQP)                                           
C                                                                               
C CALLED BY SUB. QPBASE                                                         
C                                                                               
C*************************************************************                  
C             NPQP = 0   FOR  P OPERATOR                     *                  
C             NPQP = 1   FOR  Q OPERATOR                     *                  
C             NPQP = 2   FOR PMU OPERATOR                    *                  
C             NPQP = 3   FOR  H  OPERATOR                    *                  
C               MU = -2,  -1,  0,  +1,  +2                   *                  
C*************************************************************                  
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /TQQAS / T11QQ(MD2,MD2,-2:2,6),T20QQ(MD2,MD2,-2:2,6),              
     1                T02QQ(MD2,MD2,-2:2,6)                                     
      COMMON /UVUV  / ESP(NMD,6),EQP(NMD,6),U(NMD,6),V(NMD,6),                  
     1                KQN(NMD,6)                                                
      COMMON /WORK1 / G1(NMD,NMD),G2(NMD,NMD)                                   
      COMMON /DEF   / QEL(NMD,NMD,5,6),PEL(NMD,NMD,6),                          
     1                SEL(NMD,NMD,3,6),REL(NMD,NMD,3,6),                        
     2                HEL(NMD,NMD,6),AMP(NMD,NMD,6)                             
C                                                                               
      IF (NPQP.NE.1) RETURN                                                     
      MIX = NSHL + (IT-1) * MXSH(1)                                             
      NL  = NILS(NSHL,IT)                                                       
      NL2 = NL*2                                                                
C                                                                               
      MU = 0                                                                    
      DO 99 IMU = 0,4                                                           
      IF (IMU.EQ.1.OR.IMU.EQ.3) MU = MU + IMU                                   
      IF (IMU.EQ.2.OR.IMU.EQ.4) MU = MU - IMU                                   
C                                                                               
      IF (MU.EQ.0) GOTO 21                                                      
      IF (MU.EQ.1) GOTO 22                                                      
      IF (MU.EQ.-1) GOTO 23                                                     
      IF (MU.EQ.2) GOTO 24                                                      
      IF (MU.EQ.-2) GOTO 25                                                     
C                                                                               
C===================================================================            
C         Q  MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = 0                      
C===================================================================            
C                                                                               
   21 DO 111 I=1,NL                                                             
      DO 111 J=1,NL                                                             
      G1(I,J)=U(I,MIX)*QEL(I,J,3,MIX)                                           
  111 G2(I,J)=V(I,MIX)*QEL(I,J,3,MIX)                                           
      DO 112 I=1,NL                                                             
      DO 112 J=1,NL                                                             
      T11QQ(I,J,0,MIX)=G1(I,J)*U(J,MIX)-G2(I,J)*V(J,MIX)                        
      T11QQ(I+NL,J+NL,0,MIX)=T11QQ(I,J,0,MIX)                                   
      T20QQ(I,J+NL,0,MIX)=-G1(I,J)*V(J,MIX)-G2(I,J)*U(J,MIX)                    
      T20QQ(I+NL,J,0,MIX)=-T20QQ(I,J+NL,0,MIX)                                  
      T02QQ(I,J+NL,0,MIX)=T20QQ(I,J+NL,0,MIX)                                   
  112 T02QQ(I+NL,J,0,MIX)=-T02QQ(I,J+NL,0,MIX)                                  
      GOTO 98                                                                   
C                                                                               
C===================================================================            
C         Q  MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = +1                     
C===================================================================            
C                                                                               
   22 DO 75 I=1,NL                                                              
      DO 75 J=1,NL                                                              
      G1(I,J)=U(I,MIX)*QEL(I,J,2,MIX)                                           
   75 G2(I,J)=V(I,MIX)*QEL(I,J,2,MIX)                                           
      DO 76 I=1,NL                                                              
      DO 76 J=1,NL                                                              
      T11QQ(I+NL,J,1,MIX)=G1(I,J)*U(J,MIX)-G2(I,J)*V(J,MIX)                     
   76 T20QQ(I+NL,J+NL,1,MIX)=-G1(I,J)*V(J,MIX)-G2(I,J)*U(J,MIX)                 
C                                                                               
      DO 77 I=1,NL                                                              
      DO 77 J=1,NL                                                              
      G1(I,J)=U(I,MIX)*QEL(I,J,4,MIX)                                           
   77 G2(I,J)=V(I,MIX)*QEL(I,J,4,MIX)                                           
      DO 78 I=1,NL                                                              
      DO 78 J=1,NL                                                              
      T11QQ(I,J+NL,1,MIX)=G1(I,J)*U(J,MIX)-G2(I,J)*V(J,MIX)                     
      T20QQ(I,J,1,MIX)=G1(I,J)*V(J,MIX)+G2(I,J)*U(J,MIX)                        
      T02QQ(I,J,1,MIX)=T20QQ(I+NL,J+NL,1,MIX)                                   
   78 T02QQ(I+NL,J+NL,1,MIX)=T20QQ(I,J,1,MIX)                                   
      GOTO 98                                                                   
C                                                                               
C===================================================================            
C         Q  MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = -1                     
C===================================================================            
C                                                                               
   23 DO 96 I=1,NL                                                              
      DO 96 J=1,NL                                                              
      T11QQ(I+NL,J,-1,MIX)=T11QQ(I,J+NL,1,MIX)                                  
   96 T11QQ(I,J+NL,-1,MIX)=T11QQ(I+NL,J,1,MIX)                                  
C                                                                               
      DO 97 I=1,NL                                                              
      DO 97 J=1,NL                                                              
      T20QQ(I,J,-1,MIX)=-T20QQ(I+NL,J+NL,1,MIX)                                 
   97 T20QQ(I+NL,J+NL,-1,MIX)=-T20QQ(I,J,1,MIX)                                 
C                                                                               
      DO 91 I=1,NL                                                              
      DO 91 J=1,NL                                                              
      T02QQ(I,J,-1,MIX)=-T02QQ(I+NL,J+NL,1,MIX)                                 
   91 T02QQ(I+NL,J+NL,-1,MIX)=-T02QQ(I,J,1,MIX)                                 
      GOTO 98                                                                   
C                                                                               
C===================================================================            
C         Q  MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = +2                     
C===================================================================            
C                                                                               
   24 DO 35 I=1,NL                                                              
      DO 35 J=1,NL                                                              
      G1(I,J)=U(I,MIX)*QEL(I,J,5,MIX)                                           
   35 G2(I,J)=V(I,MIX)*QEL(I,J,5,MIX)                                           
      DO 36 I=1,NL                                                              
      DO 36 J=1,NL                                                              
      T11QQ(I,J,2,MIX)=G1(I,J)*U(J,MIX)-G2(I,J)*V(J,MIX)                        
   36 T20QQ(I,J+NL,2,MIX)=-G1(I,J)*V(J,MIX)-G2(I,J)*U(J,MIX)                    
C                                                                               
      DO 39 I=1,NL                                                              
      DO 39 J=1,NL                                                              
      G1(I,J)=U(I,MIX)*QEL(I,J,1,MIX)                                           
   39 G2(I,J)=V(I,MIX)*QEL(I,J,1,MIX)                                           
      DO 38 I=1,NL                                                              
      DO 38 J=1,NL                                                              
      T11QQ(I+NL,J+NL,2,MIX)=G1(I,J)*U(J,MIX)-G2(I,J)*V(J,MIX)                  
      T20QQ(I+NL,J,2,MIX)=G1(I,J)*V(J,MIX)+G2(I,J)*U(J,MIX)                     
      T02QQ(I+NL,J,2,MIX)=-T20QQ(I,J+NL,2,MIX)                                  
   38 T02QQ(I,J+NL,2,MIX)=-T20QQ(I+NL,J,2,MIX)                                  
      GOTO 98                                                                   
C                                                                               
C===================================================================            
C         Q  MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = -2                     
C===================================================================            
C                                                                               
   25 DO 106 I=1,NL                                                             
      DO 106 J=1,NL                                                             
      T11QQ(I,J,-2,MIX)=T11QQ(I+NL,J+NL,2,MIX)                                  
  106 T11QQ(I+NL,J+NL,-2,MIX)=T11QQ(I,J,2,MIX)                                  
C                                                                               
      DO 107 I=1,NL                                                             
      DO 107 J=1,NL                                                             
      T20QQ(I+NL,J,-2,MIX)=-T20QQ(I,J+NL,2,MIX)                                 
  107 T20QQ(I,J+NL,-2,MIX)=-T20QQ(I+NL,J,2,MIX)                                 
C                                                                               
      DO 108 I=1,NL                                                             
      DO 108 J=1,NL                                                             
      T02QQ(I+NL,J,-2,MIX)=-T02QQ(I,J+NL,2,MIX)                                 
  108 T02QQ(I,J+NL,-2,MIX)=-T02QQ(I+NL,J,2,MIX)                                 
C                                                                               
   98 CONTINUE                                                                  
C                                                                               
   99 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE PQASI(IT,NSHL,NPQP)                                            
C                                                                               
C CALLED BY SUB. QPBASE                                                         
C                                                                               
C*************************************************************                  
C             NPQP = 0   FOR  P OPERATOR                     *                  
C             NPQP = 1   FOR  Q OPERATOR                     *                  
C             NPQP = 2   FOR PMU OPERATOR                    *                  
C               MU = -2,  -1,  0,  +1,  +2                   *                  
C*************************************************************                  
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /TPQAS / T11PQ(MD2,MD2,6),T20PQ(MD2,MD2,6),                        
     1                T02PQ(MD2,MD2,6)                                          
      COMMON /UVUV  / ESP(NMD,6),EQP(NMD,6),U(NMD,6),V(NMD,6),                  
     1                KQN(NMD,6)                                                
C                                                                               
      IF (NPQP.NE.0) RETURN                                                     
      MU = 0                                                                    
      MIX = NSHL + (IT-1) * MXSH(1)                                             
      NL  = NILS(NSHL,IT)                                                       
      NL2 = NL*2                                                                
C                                                                               
C===================================================================            
C         P  MATRIX  IN  QUASITEILCHENZUSTAND                                   
C===================================================================            
C                                                                               
      DO 12 I=1,NL                                                              
      T11PQ(I,I,MIX)=-U(I,MIX)*V(I,MIX)                                         
      T11PQ(I+NL,I+NL,MIX)=T11PQ(I,I,MIX)                                       
      T20PQ(I,I+NL,MIX)=V(I,MIX)*V(I,MIX)                                       
      T20PQ(I+NL,I,MIX)=-T20PQ(I,I+NL,MIX)                                      
      T02PQ(I,I+NL,MIX)=-U(I,MIX)*U(I,MIX)                                      
   12 T02PQ(I+NL,I,MIX)=-T02PQ(I,I+NL,MIX)                                      
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE PMQASI(IT,NSHL,NPQP)                                           
C                                                                               
C CALLED BY SUB. QPBASE                                                         
C                                                                               
C*************************************************************                  
C             NPQP = 0   FOR  P OPERATOR                     *                  
C             NPQP = 1   FOR  Q OPERATOR                     *                  
C             NPQP = 2   FOR PMU OPERATOR                    *                  
C               MU = -2,  -1,  0,  +1,  +2                   *                  
C*************************************************************                  
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /TPMUQA/ T11PMU(MD2,MD2,-2:2,6),T20PMU(MD2,MD2,-2:2,6),            
     1                T02PMU(MD2,MD2,-2:2,6)                                    
      COMMON /UVUV  / ESP(NMD,6),EQP(NMD,6),U(NMD,6),V(NMD,6),                  
     1                KQN(NMD,6)                                                
      COMMON /WORK1 / G1(NMD,NMD),G2(NMD,NMD)                                   
      COMMON /DEF   / QEL(NMD,NMD,5,6),PEL(NMD,NMD,6),                          
     1                SEL(NMD,NMD,3,6),REL(NMD,NMD,3,6),                        
     2                HEL(NMD,NMD,6),AMP(NMD,NMD,6)                             
C                                                                               
      IF (NPQP.NE.2) RETURN                                                     
      MIX = NSHL + (IT-1) * MXSH(1)                                             
      NL  = NILS(NSHL,IT)                                                       
      NL2 = NL*2                                                                
C                                                                               
      MU = 0                                                                    
      DO 99 IMU = 0,4                                                           
      IF (IMU.EQ.1.OR.IMU.EQ.3) MU = MU + IMU                                   
      IF (IMU.EQ.2.OR.IMU.EQ.4) MU = MU - IMU                                   
C                                                                               
      IF (MU.EQ.0) GOTO 21                                                      
      IF (MU.EQ.1) GOTO 22                                                      
      IF (MU.EQ.-1) GOTO 23                                                     
      IF (MU.EQ.2) GOTO 24                                                      
      IF (MU.EQ.-2) GOTO 25                                                     
C                                                                               
C=============================================================                  
C       PMU  MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = 0                      
C=============================================================                  
C                                                                               
   21 DO 311 I=1,NL                                                             
      DO 311 J=1,NL                                                             
      G1(I,J)=U(I,MIX)*QEL(I,J,3,MIX)                                           
  311 G2(I,J)=V(I,MIX)*QEL(I,J,3,MIX)                                           
      DO 312 I=1,NL                                                             
      DO 312 J=1,NL                                                             
      T11PMU(I,J,0,MIX)=-G2(I,J)*U(J,MIX)                                       
      T11PMU(I+NL,J+NL,0,MIX)=T11PMU(I,J,0,MIX)                                 
      T20PMU(I,J+NL,0,MIX)=G2(I,J)*V(J,MIX)                                     
      T20PMU(I+NL,J,0,MIX)=-T20PMU(I,J+NL,0,MIX)                                
      T02PMU(I,J+NL,0,MIX)=-G1(I,J)*U(J,MIX)                                    
  312 T02PMU(I+NL,J,0,MIX)=-T02PMU(I,J+NL,0,MIX)                                
      GOTO 98                                                                   
C                                                                               
C=============================================================                  
C        PMU MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = +1                     
C=============================================================                  
C                                                                               
   22 DO 375 I=1,NL                                                             
      DO 375 J=1,NL                                                             
      G1(I,J)=U(I,MIX)*QEL(I,J,2,MIX)                                           
  375 G2(I,J)=V(I,MIX)*QEL(I,J,2,MIX)                                           
      DO 376 I=1,NL                                                             
      DO 376 J=1,NL                                                             
      T11PMU(I,J+NL,1,MIX)=G2(I,J)*U(J,MIX)                                     
      T20PMU(I,J,1,MIX)=G2(I,J)*V(J,MIX)                                        
  376 T02PMU(I+NL,J+NL,1,MIX)=-G1(I,J)*U(J,MIX)                                 
C                                                                               
      DO 377 I=1,NL                                                             
      DO 377 J=1,NL                                                             
      G1(I,J)=U(I,MIX)*QEL(I,J,4,MIX)                                           
  377 G2(I,J)=V(I,MIX)*QEL(I,J,4,MIX)                                           
      DO 378 I=1,NL                                                             
      DO 378 J=1,NL                                                             
      T11PMU(I+NL,J,1,MIX)=G2(I,J)*U(J,MIX)                                     
      T20PMU(I+NL,J+NL,1,MIX)=-G2(I,J)*V(J,MIX)                                 
  378 T02PMU(I,J,1,MIX)=G1(I,J)*U(J,MIX)                                        
      GOTO 98                                                                   
C                                                                               
C=============================================================                  
C        PMU MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = -1                     
C=============================================================                  
C                                                                               
   23 DO 476 I=1,NL                                                             
      DO 476 J=1,NL                                                             
      T11PMU(I,J+NL,-1,MIX)=T11PMU(I+NL,J,1,MIX)                                
  476 T11PMU(I+NL,J,-1,MIX)=T11PMU(I,J+NL,1,MIX)                                
C                                                                               
      DO 477 I=1,NL                                                             
      DO 477 J=1,NL                                                             
      T20PMU(I,J,-1,MIX)=-T20PMU(I+NL,J+NL,1,MIX)                               
  477 T20PMU(I+NL,J+NL,-1,MIX)=-T20PMU(I,J,1,MIX)                               
C                                                                               
      DO 478 I=1,NL                                                             
      DO 478 J=1,NL                                                             
      T02PMU(I,J,-1,MIX)=-T02PMU(I+NL,J+NL,1,MIX)                               
  478 T02PMU(I+NL,J+NL,-1,MIX)=-T02PMU(I,J,1,MIX)                               
      GOTO 98                                                                   
C                                                                               
C=============================================================                  
C        PMU MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = +2                     
C=============================================================                  
C                                                                               
   24 DO 335 I=1,NL                                                             
      DO 335 J=1,NL                                                             
      G1(I,J)=U(I,MIX)*QEL(I,J,5,MIX)                                           
  335 G2(I,J)=V(I,MIX)*QEL(I,J,5,MIX)                                           
      DO 336 I=1,NL                                                             
      DO 336 J=1,NL                                                             
      T11PMU(I+NL,J+NL,2,MIX)=-G2(I,J)*U(J,MIX)                                 
      T20PMU(I+NL,J,2,MIX)=-G2(I,J)*V(J,MIX)                                    
  336 T02PMU(I,J+NL,2,MIX)=-G1(I,J)*U(J,MIX)                                    
C                                                                               
      DO 359 I=1,NL                                                             
      DO 359 J=1,NL                                                             
      G1(I,J)=U(I,MIX)*QEL(I,J,1,MIX)                                           
  359 G2(I,J)=V(I,MIX)*QEL(I,J,1,MIX)                                           
      DO 340 I=1,NL                                                             
      DO 340 J=1,NL                                                             
      T11PMU(I,J,2,MIX)=-G2(I,J)*U(J,MIX)                                       
      T20PMU(I,J+NL,2,MIX)=G2(I,J)*V(J,MIX)                                     
  340 T02PMU(I+NL,J,2,MIX)=G1(I,J)*U(J,MIX)                                     
      GOTO 98                                                                   
C                                                                               
C=============================================================                  
C        PMU MATRIX  IN  QUASITEILCHENZUSTAND  MIT  MU = -2                     
C=============================================================                  
C                                                                               
   25 DO 436 I=1,NL                                                             
      DO 436 J=1,NL                                                             
      T11PMU(I,J,-2,MIX)=T11PMU(I+NL,J+NL,2,MIX)                                
  436 T11PMU(I+NL,J+NL,-2,MIX)=T11PMU(I,J,2,MIX)                                
C                                                                               
      DO 437 I=1,NL                                                             
      DO 437 J=1,NL                                                             
      T20PMU(I,J+NL,-2,MIX)=-T20PMU(I+NL,J,2,MIX)                               
  437 T20PMU(I+NL,J,-2,MIX)=-T20PMU(I,J+NL,2,MIX)                               
C                                                                               
      DO 438 I=1,NL                                                             
      DO 438 J=1,NL                                                             
      T02PMU(I,J+NL,-2,MIX)=-T02PMU(I+NL,J,2,MIX)                               
  438 T02PMU(I+NL,J,-2,MIX)=-T02PMU(I,J+NL,2,MIX)                               
C                                                                               
   98 CONTINUE                                                                  
   99 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE KERNEL                                                         
C                                                                               
C CALLED BY MAIN                                                                
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /GAUSSL/ X(NNDEG),W(NNDEG),WEIT(NNDEG),NDEG                        
      COMMON /MAXMAX/ EMAX(3),KMAX(3),NMAX(3),JMAX                              
      COMMON /TRACE1/ ST02B(2,12),ST20B(2,12),                                  
     1                AAT1(NNRANK,12,2),AAT2(NNRANK,12,2),                      
     2                T1AA(NNRANK,12,2),T2AA(NNRANK,12,2)                       
C                                                                               
      WRITE(1,300)                                                              
  300 FORMAT(/' IBET',5X,'BETA',7X,'<H[R]>',9X,'<R>')                           
C                                                                               
C==================================                                             
C     CALCULATION OF THE KERNELS                                                
C==================================                                             
C                                                                               
      IBETI=1                                                                   
      IBETF=NDEG/2                                                              
      DO 10 IBET=IBETI,IBETF                                                    
C                                                                               
      WRITE(*,100) IBET                                                         
  100 FORMAT(/,20('*'),' BETA LOOP :',I3)                                       
C                                                                               
      COSBET=X(IBET)                                                            
      WEIT(IBET)=W(IBET)                                                        
      CALL SDFSET(JMAX,COSBET)                                                  
C                                                                               
      DO 3 IT=1,2                                                               
      MS=MXSH(IT)                                                               
      DO 3 NSHL=1,MS                                                            
      CALL ROTATE(IT,NSHL)                                                      
      CALL CONTRA(IT,NSHL,OLAP)                                                 
      WEIT(IBET) = OLAP*WEIT(IBET)                                              
    3 CONTINUE                                                                  
C                                                                               
      DO 13 JBET = 1,2                                                          
C                                                                               
      DO 5 IT= 1,2                                                              
      DO 5 IS= 1,12                                                             
      ST02B(IT,IS) = 0.                                                         
    5 ST20B(IT,IS) = 0.                                                         
C                                                                               
      DO 11 IT=1,2                                                              
      MS=MXSH(IT)                                                               
      DO 11 NSHL=1,MS                                                           
C                                                                               
      IS = 0                                                                    
      DO 1 NPQP = 0,3                                                           
      DO 1 MU  = -2,2                                                           
      IF (NPQP.EQ.0.AND.MU.NE.0) GOTO 1                                         
      IF (NPQP.EQ.3.AND.MU.NE.0) GOTO 1                                         
      IS = IS + 1                                                               
      CALL SPUR(JBET,IT,NSHL,NPQP,MU,IS)                                        
    1 CONTINUE                                                                  
   11 CONTINUE                                                                  
C                                                                               
      CALL HAMVAC(HR00)                                                         
C---------------------------------------------------                            
C     BERECHNUNG VON MATRIZEN IM INTRUDER ORBIT                                 
C---------------------------------------------------                            
      CALL GETCON(JBET)                                                         
      CALL INNER1                                                               
      CALL INNER2                                                               
      CALL ELEM00(IBET,JBET,HR00)                                               
   13 CONTINUE                                                                  
C                                                                               
      BETA=DACOS(COSBET)                                                        
      OLAP=WEIT(IBET)/W(IBET)                                                   
      WRITE(1,200) IBET,BETA,HR00,OLAP                                          
  200 FORMAT(I4,F11.5,2E15.5)                                                   
C                                                                               
   10 CONTINUE                                                                  
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE INTEG                                                          
C                                                                               
C CALLED BY MAIN                                                                
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
      PARAMETER (NNDEH=NNDEG/2)                                                 
      PARAMETER (NIE=127)                                                       
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      DIMENSION WDF(-NKMAX:NKMAX,-NKMAX:NKMAX,NNDEH),INDQP(NNRANT)              
      DIMENSION EYRAST(NJOTMX),EBAND(NJOTMX,NNRANT),NBAND(NJOTMX)               
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /GAUSSL/ X(NNDEG),W(NNDEG),WEIT(NNDEG),NDEG                        
      COMMON /GFVIV / IV(-NIE:NIE)                                              
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /JJJJJJ/ JJ(NJOTMX),JOTMAX                                         
      COMMON /MAXMAX/ EMAX(3),KMAX(3),NMAX(3),JMAX                              
      COMMON /NORHAM/ WK1(NNRANT,NNRANT,2,NNDEH),                               
     1                WK2(NNRANT,NNRANT,2,NNDEH)                                
      COMMON /TQB1  / NOQ(2),NOQR(2),IOQ(2,MD2),IOQR(2,MD2),                    
     1                NTQ(3),NTQH(2)                                            
      COMMON /TQB2  / ETQ(NNRANT),KTQ(NNRANT),ITQ(NNRANT,3),NRANK               
      COMMON /TQB3  / IOQB(2,MD2),ITQB(2,NNATA,2),NATA(2)                       
      COMMON /TRACE1/ ST02B(2,12),ST20B(2,12),                                  
     1                AAT1(NNRANK,12,2),AAT2(NNRANK,12,2),                      
     2                T1AA(NNRANK,12,2),T2AA(NNRANK,12,2)                       
      COMMON /TRACE2/ AT1A(NNATA,12,2),AT2A(NNATA,12,2)                         
      COMMON /HND   / HH(NNRANT,NNRANT),HN(NNRANT,NNRANT),                      
     1                UU(NNRANT,NNRANT),VV(NNRANT,NNRANT),                      
     2                FF(NNRANT,NNRANT),EH(NNRANT),EN(NNRANT),                  
     3                XX(NNRANT)                                                
      COMMON/NHELEM / NHOUT                                                     
C                                                                               
      DATA EPS/0.001/                                                           
C                                                                               
      IBETI=1                                                                   
      IBETF=NDEG/2                                                              
      MAXJOT=0                                                                  
      DO 30 JOT = 1,JOTMAX                                                      
      NBAND(JOT)=0                                                              
      IJOT=JJ(JOT)                                                              
      AJOT=IJOT                                                                 
      KM=MAX(KMAX(1),KMAX(2))                                                   
      KM=MAX(KM,KMAX(3))                                                        
      KM=MIN(KM,IJOT)                                                           
      WRITE(*,201) IJOT                                                         
      WRITE(1,201) IJOT                                                         
  201 FORMAT(/,'******** SPIN =',I3,' ********')                                
C------------------------------------------------                               
C     WIGNER-FUNKTION                                                           
C------------------------------------------------                               
      DO 31 IBET=IBETI,IBETF                                                    
      COSBET=X(IBET)                                                            
      DO 32 K1=-KM,KM                                                           
      AK1=K1                                                                    
      DO 32 K2=-KM,K1                                                           
      AK2=K2                                                                    
      SN1=IV(K2-K1)                                                             
      WDF(K1,K2,IBET)=SDFUNC(AJOT,AK1,AK2,COSBET)                               
      WDF(K2,K1,IBET)=SN1*WDF(K1,K2,IBET)                                       
   32 CONTINUE                                                                  
   31 CONTINUE                                                                  
C------------------------------------------------                               
C     LOOP OVER STATES                                                          
C------------------------------------------------                               
      I1 = 0                                                                    
      DO 33 IP1 = 1,NRANK                                                       
      K1=KTQ(IP1)                                                               
      IF (IABS(K1).GT.IJOT) GOTO 33                                             
      SIK=IV(IJOT+K1)                                                           
      I1 = I1 + 1                                                               
      INDQP(I1) = IP1                                                           
      I2 = 0                                                                    
      DO 34 IP2 = 1,IP1                                                         
      K2=KTQ(IP2)                                                               
      IF (IABS(K2).GT.IJOT) GOTO 34                                             
      I2 = I2 + 1                                                               
C                                                                               
      SUMH=0.D0                                                                 
      SUMN=0.D0                                                                 
      DO 35 IBET=IBETI,IBETF                                                    
      WW1 = WEIT(IBET) * WDF(K1, K2,IBET)                                       
      WW2 = WEIT(IBET) * WDF(K1,-K2,IBET) * SIK                                 
      SUMN=SUMN+WW1*WK1(IP1,IP2,1,IBET)+WW2*WK1(IP1,IP2,2,IBET)                 
      SUMH=SUMH+WW1*WK2(IP1,IP2,1,IBET)+WW2*WK2(IP1,IP2,2,IBET)                 
   35 CONTINUE                                                                  
      HN(I1,I2)=SUMN                                                            
      HH(I1,I2)=SUMH                                                            
C                                                                               
   34 CONTINUE                                                                  
   33 CONTINUE                                                                  
      N=I1                                                                      
      IF (N.EQ.0) THEN                                                          
       PRINT*, ' STOP: DIMENSION OF THE PROJECTED MATRIX IS 0'                  
       STOP                                                                     
      ENDIF                                                                     
C----------------------------------------------------                           
C     DIAGONALIZATION                                                           
C----------------------------------------------------                           
      AVHN=0.                                                                   
      DO 40 I=1,N                                                               
      IF (HN(I,I).LE.0.) THEN                                                   
       PRINT*, ' STOP: NON-POSITIVE NORM FOR BAND =',INDQP(I),HN(I,I)           
       STOP                                                                     
      ENDIF                                                                     
      AVHN=AVHN+HN(I,I)                                                         
   40 CONTINUE                                                                  
      AVHN=AVHN/N                                                               
      TOL=EPS*AVHN                                                              
C                                                                               
      CALL HNDIAG(N,M,TOL,0)                                                    
      IF (M.LT.1) GOTO 36                                                       
C                                                                               
      IF (M.LT.N) THEN                                                          
       M1=M+1                                                                   
       DO 41 I=M1,N                                                             
       EH(I)=0.                                                                 
   41  CONTINUE                                                                 
      ENDIF                                                                     
C                                                                               
      WRITE(1,101)                                                              
  101 FORMAT(/' BAND KQP=',3X,'EQP=',4X,'HH/HN=',4X,'HH=',9X,'HN=',             
     +     6X,' EIGEN NORM  STATE  ENERGY')                                     
      DO 43 I = 1,N                                                             
      K=INDQP(I)                                                                
      KK=KTQ(K)                                                                 
      EB=HH(I,I)/HN(I,I)                                                        
      EBAND(JOT,K)=EB                                                           
      WRITE(1,102) K,KK,ETQ(K),EB,HH(I,I),HN(I,I),EN(I),I,EH(I)                 
  102 FORMAT(I4,I5,2F9.4,3E12.4,I6,F9.4)                                        
   43 CONTINUE                                                                  
C                                                                               
      IF (NHOUT.NE.0) THEN                                                      
       WRITE(4,201) IJOT                                                        
       WRITE(4,*) 'BAND  NORM MATRIX ELEMENTS'                                  
       DO 45 I=1,N                                                              
       WRITE(4,500) INDQP(I),(HN(I,J),J=1,I)                                    
   45  CONTINUE                                                                 
       WRITE(4,*) 'BAND  HAMILTONIAN ELEMENTS'                                  
       DO 46 I=1,N                                                              
       WRITE(4,500) INDQP(I),(HH(I,J),J=1,I)                                    
   46  CONTINUE                                                                 
  500  FORMAT(I4,99(1PD12.3))                                                   
      ENDIF                                                                     
C                                                                               
      IF (JOT.EQ.1) EGS=EH(1)                                                   
      EYRAST(JOT) = EH(1)                                                       
      NBAND(JOT)=N                                                              
      MM=MIN(M,6)                                                               
      WRITE(1,204) MM,M,(EH(I),I=1,MM)                                          
  204 FORMAT(/' ENERGIES (LOWEST',I2,' OUT OF ',I2,'):',6F7.3)                  
      WRITE(1,205) TOL,EN(N),AVHN,EN(1)                                         
  205 FORMAT(' CUTOFF OF THE NORM EIGENVALUE:',E12.4,' (',3E11.4,')')           
C                                                                               
      IF (EN(N).LT.0.) THEN                                                     
        NN=0                                                                    
        DO 47 I=N,1,-1                                                          
        IF (EN(I).GT.0.) GOTO 49                                                
        NN=NN+1                                                                 
   47   CONTINUE                                                                
   49   WRITE(1,206) NN,EN(N),EN(N-NN+1)                                        
  206   FORMAT(' WARNING:',I3,' NEGATIVE NORM EIGENVALUES BETWEEN',             
     +         E12.4,' AND',E12.4)                                              
      ENDIF                                                                     
C                                                                               
      MAXJOT=MAXJOT+1                                                           
   30 CONTINUE                                                                  
C------------------------                                                       
C     END LOOP OVER J                                                           
C------------------------                                                       
   36 JOTMAX=MAXJOT                                                             
      IF (JOTMAX.EQ.0) STOP                                                     
C                                                                               
      WRITE(1,108)                                                              
      DO 37 JOT = 1,JOTMAX                                                      
      IF (NBAND(JOT).EQ.0) GOTO 37                                              
      IJOT = JJ(JOT)                                                            
      AJOT = IJOT                                                               
      IF (JOT.EQ.1) THEN                                                        
        WRITE(1,109) IJOT,EGS                                                   
      ELSE                                                                      
        DE = EYRAST(JOT)-EYRAST(JOT-1)                                          
        WRITE(1,109) IJOT,EYRAST(JOT),DE                                        
      ENDIF                                                                     
   37 CONTINUE                                                                  
C                                                                               
C===============================================================                
C     EE_OUT2 (WRITE 2) AND EE_OUT3 (WRITE 3) FOR XMGR GRAPHICS                 
C===============================================================                
C                                                                               
      WRITE(2,105) EGS                                                          
      DO 60 J=1,JOTMAX                                                          
      AJOT=JJ(J)                                                                
      WRITE(2,103) AJOT,EYRAST(J)-EGS                                           
   60 CONTINUE                                                                  
C                                                                               
c      NN=MIN(NRANK-1,29)                                                       
      NN=NRANK-1                                                                
      DO 62 I=0,NN                                                              
      IF (I.EQ.0) THEN                                                          
        K=NRANK                                                                 
      ELSE                                                                      
        K=INDQP(I)                                                              
      ENDIF                                                                     
      WRITE(2,106) K,KTQ(K)                                                     
      DO 61 J=1,JOTMAX                                                          
      IF (NBAND(J).EQ.0) GOTO 61                                                
      IJOT=JJ(J)                                                                
      AJOT=IJOT                                                                 
      IF (IABS(KTQ(K)).LE.IJOT) WRITE(2,103) AJOT,EBAND(J,K)-EGS                
   61 CONTINUE                                                                  
   62 CONTINUE                                                                  
C                                                                               
      IF (JJ(2)-JJ(1).EQ.2) THEN                                                
       I1=2                                                                     
       I2=1                                                                     
      ELSE                                                                      
       I1=3                                                                     
       I2=2                                                                     
      ENDIF                                                                     
      WRITE(3,107)                                                              
      DO 70 J=I1,JOTMAX,I2                                                      
      AJOT=JJ(J)                                                                
      DE=EYRAST(J)-EYRAST(J-I2)                                                 
      OM=0.5*DE                                                                 
      AM=(2.*AJOT-1.)/OM                                                        
      OM=OM*OM                                                                  
      WRITE(3,104) OM,AM                                                        
   70 CONTINUE                                                                  
C                                                                               
  103 FORMAT(F5.1,2F10.4)                                                       
  104 FORMAT(3F10.4)                                                            
  105 FORMAT('# ENRGIES RELATIVE TO',F8.4,/,'# YRAST')                          
  106 FORMAT('&',/,'# BAND:',I3,3X,'K =',I3)                                    
  107 FORMAT('# BACKBENDING PLOT')                                              
  108 FORMAT(//' SPIN   YRAST-E   DELTA-E'/                                     
     +         ' ~~~~   ~~~~~~~   ~~~~~~~')                                     
  109 FORMAT(I5,2F10.4)                                                         
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE HAMVAC(HR00)                                                   
C                                                                               
C CALLED BY SUB. INTEG                                                          
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /BCSBCS/ ANP(2),ALAMB(2),DELTA(2),EP(NNBCS,2),                     
     1                HSP(2),QQ0(2),HH0(2),QP0(2),LNP(2)                        
      COMMON /PARAM / AG(2),AGQ(2),DELE,HOLNG2(2),QCL,HFB,                      
     1                AX(3),AXH(3),XQQ(2),GQP(2)                                
      COMMON /TRACE1/ ST02B(2,12),ST20B(2,12),                                  
     1                AAT1(NNRANK,12,2),AAT2(NNRANK,12,2),                      
     2                T1AA(NNRANK,12,2),T2AA(NNRANK,12,2)                       
C                                                                               
      WRITE(*,100)                                                              
  100 FORMAT(20('*'),' HAMVAC')                                                 
C                                                                               
      QMU1 = 0.                                                                 
      PMU1 = 0.                                                                 
      DO 1 MU = -2,2                                                            
      IS = 4 + MU                                                               
      QMU1 = QMU1 + AX(1)*ST20B(1,IS)*ST02B(1,IS)                               
     1            + AX(2)*ST20B(2,IS)*ST02B(2,IS)                               
     2            + AX(3)*(ST20B(1,IS)*ST02B(2,IS)                              
     3                    +ST20B(2,IS)*ST02B(1,IS))                             
      IS = 9 + MU                                                               
      PMU1 = PMU1 + AGQ(1)*ST20B(1,IS)*ST02B(1,IS)                              
     +            + AGQ(2)*ST20B(2,IS)*ST02B(2,IS)                              
    1 CONTINUE                                                                  
      QMU1 = .5*QMU1                                                            
      QMU2 = XQQ(1)*ST02B(1,4)+XQQ(2)*ST02B(2,4)                                
      PMU2 = AGQ(1)*QP0(1)*(ST20B(1,9)+ST02B(1,9))                              
     +      +AGQ(2)*QP0(2)*(ST20B(2,9)+ST02B(2,9))                              
      PNP1 = AG(1)*ST20B(1,1)*ST02B(1,1)                                        
     +      +AG(2)*ST20B(2,1)*ST02B(2,1)                                        
      PNP2 = DELTA(1)*(ST20B(1,1)+ST02B(1,1))                                   
     +      +DELTA(2)*(ST20B(2,1)+ST02B(2,1))                                   
      HHHH = ST02B(1,12)+ST02B(2,12)                                            
C                                                                               
      HR00 = HHHH-QMU1-PMU1-PNP1-QMU2-PMU2-PNP2                                 
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE GETCON(JBET)                                                   
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /TINQAS/ AA(MD2,MD2,2),BB(MD2,MD2,2),CC(MD2,MD2,2)                 
      COMMON /CON   / CA1(NMD,NMD,6),CA2(NMD,NMD,6),CB1(NMD,NMD,6),             
     1                CB2(NMD,NMD,6),CC1(NMD,NMD,6),CC2(NMD,NMD,6)              
C                                                                               
      DO 10 IT=1,2                                                              
      MIX = MIXINT(IT)                                                          
      NL  = NLINT(IT)                                                           
      IF (JBET.EQ.1) THEN                                                       
        CALL LARGE(MD2,NMD,NL,CA1(1,1,MIX),CA2(1,1,MIX),AA(1,1,IT),+1)          
        CALL LARGE(MD2,NMD,NL,CB1(1,1,MIX),CB2(1,1,MIX),BB(1,1,IT),+1)          
        CALL LARGE(MD2,NMD,NL,CC1(1,1,MIX),CC2(1,1,MIX),CC(1,1,IT),+1)          
      ELSE                                                                      
        CALL LARGE(MD2,NMD,NL,CA1(1,1,MIX),CA2(1,1,MIX),AA(1,1,IT), 0)          
        CALL LARGE(MD2,NMD,NL,CB1(1,1,MIX),CB2(1,1,MIX),BB(1,1,IT), 0)          
        CALL LARGE(MD2,NMD,NL,CC1(1,1,MIX),CC2(1,1,MIX),CC(1,1,IT),-1)          
      ENDIF                                                                     
   10 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE INNER1                                                         
C                                                                               
C CALLED BY SUB. INTEG                                                          
C                                                                               
C**************************************************************                 
C     IT: 1 -- NEUTRON ;      2 -- PROTON                     *                 
C     II: 1 -- P ; 2 -- Q0 ; 3 -- Q1 ; 4 -- Q2 ;              *                 
C                  5 -- P0 ; 6 -- P1 ; 7 -- P2                *                 
C     IS: 1 -- T11B  ; 2 -- T11C  ; 3 -- BT20B ; 4 -- BT20C;  *                 
C         5 -- CT20C ; 6 -- BT02B ; 7 -- BT02C ; 8 -- CT02C   *                 
C**************************************************************                 
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /PARAM / AG(2),AGQ(2),DELE,HOLNG2(2),QCL,HFB,                      
     1                AX(3),AXH(3),XQQ(2),GQP(2)                                
      COMMON /TQB1  / NOQ(2),NOQR(2),IOQ(2,MD2),IOQR(2,MD2),                    
     1                NTQ(3),NTQH(2)                                            
      COMMON /TQB2  / ETQ(NNRANT),KTQ(NNRANT),ITQ(NNRANT,3),NRANK               
      COMMON /TQB3  / IOQB(2,MD2),ITQB(2,NNATA,2),NATA(2)                       
      COMMON /TPQAS / T11PQ(MD2,MD2,6),T20PQ(MD2,MD2,6),                        
     1                T02PQ(MD2,MD2,6)                                          
      COMMON /TQQAS / T11QQ(MD2,MD2,-2:2,6),T20QQ(MD2,MD2,-2:2,6),              
     1                T02QQ(MD2,MD2,-2:2,6)                                     
      COMMON /TPMUQA/ T11PMU(MD2,MD2,-2:2,6),T20PMU(MD2,MD2,-2:2,6),            
     1                T02PMU(MD2,MD2,-2:2,6)                                    
      COMMON /THQAS / T11HQ(MD2,MD2,6),T20HQ(MD2,MD2,6),                        
     1                T02HQ(MD2,MD2,6)                                          
      COMMON /TINQAS/ AA(MD2,MD2,2),BB(MD2,MD2,2),CC(MD2,MD2,2)                 
      COMMON /TRACE1/ ST02B(2,12),ST20B(2,12),                                  
     1                AAT1(NNRANK,12,2),AAT2(NNRANK,12,2),                      
     2                T1AA(NNRANK,12,2),T2AA(NNRANK,12,2)                       
      COMMON /TRACE2/ AT1A(NNATA,12,2),AT2A(NNATA,12,2)                         
C                                                                               
      DIMENSION W1(MD2,MD2),W2(MD2,MD2),W3(MD2,MD2),                            
     1          W4(MD2,MD2),W5(MD2,MD2)                                         
      DIMENSION T11(MD2,MD2),T20(MD2,MD2),T02(MD2,MD2)                          
      DIMENSION IOQ1(NNRANK),IOQ2(NNRANK)                                       
C                                                                               
      WRITE(*,101)                                                              
  101 FORMAT(20('*'),' INNER1')                                                 
C                                                                               
      DO 92 IT=1,2                                                              
      MIX = MIXINT(IT)                                                          
      NL  = NLINT(IT)                                                           
      NL2 = NL * 2                                                              
      NOQR2 = NOQR(IT)*2                                                        
C                                                                               
      DO 90 I = 1,NTQ(IT)                                                       
      IF (IT.EQ.1) THEN                                                         
        I1 = ITQ(I,1)                                                           
        I2 = ITQ(I,2)                                                           
      ELSE                                                                      
        I1 = ITQ(NTQ(1)+I,1)                                                    
        I2 = ITQ(NTQ(1)+I,2)                                                    
      ENDIF                                                                     
      IOQ1(I)=I1                                                                
      IOQ2(I)=I2                                                                
   90 CONTINUE                                                                  
C                                                                               
C     CHOOSE THE OPERATORS ACCORDING TO THE NUMBER II                           
C                                                                               
      IS = 0                                                                    
      DO 91 II = 0,3                                                            
      DO 91 MU = -2,2                                                           
      IF (II.EQ.0.AND.MU.NE.0) GOTO 91                                          
      IF (II.EQ.3.AND.MU.NE.0) GOTO 91                                          
      IS = IS + 1                                                               
C                                                                               
      IF (II.EQ.0) THEN                                                         
         DO 201 I=1,NL2                                                         
         DO 201 J=1,NL2                                                         
         T11(I,J) = T11PQ(I,J,MIX)                                              
         T02(I,J) = T02PQ(I,J,MIX)                                              
         T20(I,J) = T20PQ(I,J,MIX)                                              
  201    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      IF (II.EQ.1) THEN                                                         
         DO 202 I=1,NL2                                                         
         DO 202 J=1,NL2                                                         
         T11(I,J) = T11QQ(I,J,MU,MIX)                                           
         T02(I,J) = T02QQ(I,J,MU,MIX)                                           
         T20(I,J) = T20QQ(I,J,MU,MIX)                                           
  202    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      IF (II.EQ.2) THEN                                                         
         DO  203 I=1,NL2                                                        
         DO  203 J=1,NL2                                                        
         T11(I,J) = T11PMU(I,J,MU,MIX)                                          
         T02(I,J) = T02PMU(I,J,MU,MIX)                                          
         T20(I,J) = T20PMU(I,J,MU,MIX)                                          
  203    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      IF (II.EQ.3) THEN                                                         
         DO  204 I=1,NL2                                                        
         DO  204 J=1,NL2                                                        
         T11(I,J) = T11HQ(I,J,MIX)                                              
         T02(I,J) = T02HQ(I,J,MIX)                                              
         T20(I,J) = T20HQ(I,J,MIX)                                              
  204    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
      DO 11 I1=1,NOQR2                                                          
      I = IOQB(IT,I1)                                                           
      DO 11 J=1,NL2                                                             
      SUM1=0.D0                                                                 
      SUM2=0.D0                                                                 
      DO 110 K=1,NL2                                                            
      SUM1=SUM1+BB(I,K,IT)*T02(K,J)                                             
  110 SUM2=SUM2+BB(I,K,IT)*T20(K,J)                                             
      W3(I,J)=SUM1                                                              
  11  W4(I,J)=SUM2                                                              
C                                                                               
C     ( A A T )                                                                 
C                                                                               
      DO 13 I = 1,NTQ(IT)                                                       
      I1=IOQ1(I)                                                                
      I2=IOQ2(I)                                                                
      SUM1=0.D0                                                                 
      SUM2=0.D0                                                                 
      DO 130 K=1,NL2                                                            
      SUM1=SUM1-BB(I1,K,IT)*T11(I2,K)+BB(I2,K,IT)*T11(I1,K)                     
  130 SUM2=SUM2-W3(I1,K)*BB(I2,K,IT)+W3(I2,K)*BB(I1,K,IT)                       
  13  AAT1(I,IS,IT)=SUM1+0.5D0*(SUM2+T20(I1,I2)-T20(I2,I1))                     
C                                                                               
      IF (II.EQ.3) GOTO 111                                                     
C                                                                               
C     ( A A T+ )                                                                
C                                                                               
      DO 14 I = 1,NTQ(IT)                                                       
      I1=IOQ1(I)                                                                
      I2=IOQ2(I)                                                                
      SUM1=0.D0                                                                 
      SUM2=0.D0                                                                 
      DO 140 K=1,NL2                                                            
      SUM1=SUM1-BB(I1,K,IT)*T11(K,I2)+BB(I2,K,IT)*T11(K,I1)                     
  140 SUM2=SUM2-W4(I1,K)*BB(I2,K,IT)+W4(I2,K)*BB(I1,K,IT)                       
  14  AAT2(I,IS,IT)=SUM1+0.5D0*(SUM2+T02(I1,I2)-T02(I2,I1))                     
C                                                                               
  111 CONTINUE                                                                  
C                                                                               
C     ( T A+ A+ )                                                               
C                                                                               
      DO 17 I = 1,NTQ(IT)                                                       
      I1=IOQ1(I)                                                                
      I2=IOQ2(I)                                                                
      DO 16 J=1,NL2                                                             
      SUM1=0.D0                                                                 
      SUM2=0.D0                                                                 
      DO 160 K=1,NL2                                                            
      SUM1=SUM1+CC(K,I1,IT)*T02(K,J)                                            
  160 SUM2=SUM2+CC(K,I2,IT)*T02(K,J)                                            
      W1(I1,J)=SUM1                                                             
  16  W1(I2,J)=SUM2                                                             
C                                                                               
      SUM=0.D0                                                                  
      DO 170 K=1,NL2                                                            
  170 SUM=SUM+W1(I1,K)*CC(K,I2,IT)-W1(I2,K)*CC(K,I1,IT)                         
  17  T1AA(I,IS,IT)=0.5D0*SUM                                                   
C                                                                               
      IF (II.EQ.3) GOTO 112                                                     
C                                                                               
C     ( T+ A+ A+ )                                                              
C                                                                               
      DO 19 I = 1,NTQ(IT)                                                       
      I1=IOQ1(I)                                                                
      I2=IOQ2(I)                                                                
      DO 18 J=1,NL2                                                             
      SUM1=0.D0                                                                 
      SUM2=0.D0                                                                 
      DO 180 K=1,NL2                                                            
      SUM1=SUM1+CC(K,I1,IT)*T20(K,J)                                            
  180 SUM2=SUM2+CC(K,I2,IT)*T20(K,J)                                            
      W1(I1,J)=SUM1                                                             
  18  W1(I2,J)=SUM2                                                             
C                                                                               
      SUM=0.D0                                                                  
      DO 190 K=1,NL2                                                            
  190 SUM=SUM+W1(I1,K)*CC(K,I2,IT)-W1(I2,K)*CC(K,I1,IT)                         
  19  T2AA(I,IS,IT)=0.5D0*SUM                                                   
C                                                                               
  112 CONTINUE                                                                  
C                                                                               
C     ( A T A+ )                                                                
C                                                                               
      DO 20 I1=1,NOQR2                                                          
      I = IOQB(IT,I1)                                                           
      DO 20 J=1,NL2                                                             
      SUM=0.D0                                                                  
      DO 200 K=1,NL2                                                            
  200 SUM=SUM+BB(I,K,IT)*T02(J,K)                                               
  20  W5(I,J)=SUM                                                               
C                                                                               
      DO 21 I1=1,NOQR2                                                          
      J1 = IOQB(IT,I1)                                                          
      DO 21 I2=1,NOQR2                                                          
      J2 = IOQB(IT,I2)                                                          
      SUM1=0.D0                                                                 
      SUM2=0.D0                                                                 
      DO 210 K=1,NL2                                                            
      SUM1=SUM1+W3(J1,K)*CC(K,J2,IT)                                            
  210 SUM2=SUM2+W5(J1,K)*CC(K,J2,IT)                                            
      W1(J1,J2)=SUM1                                                            
  21  W2(J1,J2)=SUM2                                                            
C                                                                               
      DO 22 I1=1,NOQR2                                                          
      J1 = IOQB(IT,I1)                                                          
      DO 22 I2=1,NOQR2                                                          
      J2 = IOQB(IT,I2)                                                          
      SUM=0.D0                                                                  
      DO 220 K=1,NL2                                                            
  220 SUM=SUM+T11(J1,K)*CC(K,J2,IT)                                             
  22  W3(J1,J2)=SUM                                                             
C                                                                               
      I = 0                                                                     
      DO 5 I1 = 1,NOQR2                                                         
      J1 = IOQB(IT,I1)                                                          
      DO 5 I2 = 1,NOQR2                                                         
      J2 = IOQB(IT,I2)                                                          
      I = I + 1                                                                 
      AT1A(I,IS,IT)=W3(J2,J1)+0.5D0*(W2(J2,J1)-W1(J2,J1))                       
    5 CONTINUE                                                                  
C                                                                               
      IF (II.EQ.3) GOTO 113                                                     
C                                                                               
C     ( A T+ A+ )                                                               
C                                                                               
      DO 25 I1=1,NOQR2                                                          
      I = IOQB(IT,I1)                                                           
      DO 25 J=1,NL2                                                             
      SUM=0.D0                                                                  
      DO 250 K=1,NL2                                                            
  250 SUM=SUM+BB(I,K,IT)*T20(J,K)                                               
  25  W5(I,J)=SUM                                                               
C                                                                               
      DO 26 I1=1,NOQR2                                                          
      J1 = IOQB(IT,I1)                                                          
      DO 26 I2=1,NOQR2                                                          
      J2 = IOQB(IT,I2)                                                          
      SUM1=0.D0                                                                 
      SUM2=0.D0                                                                 
      SUM3=0.D0                                                                 
      DO 260 K=1,NL2                                                            
      SUM1=SUM1+T11(K,J1)*CC(K,J2,IT)                                           
      SUM2=SUM2+W5(J1,K)*CC(K,J2,IT)                                            
  260 SUM3=SUM3+W4(J1,K)*CC(K,J2,IT)                                            
      W1(J1,J2)=SUM1                                                            
      W2(J1,J2)=SUM2                                                            
  26  W3(J1,J2)=SUM3                                                            
C                                                                               
      I = 0                                                                     
      DO 6 I1 = 1,NOQR2                                                         
      J1 = IOQB(IT,I1)                                                          
      DO 6 I2 = 1,NOQR2                                                         
      J2 = IOQB(IT,I2)                                                          
      I = I + 1                                                                 
      AT2A(I,IS,IT)=W1(J2,J1)+0.5D0*(W2(J2,J1)-W3(J2,J1))                       
    6 CONTINUE                                                                  
C                                                                               
  113 CONTINUE                                                                  
   91 CONTINUE                                                                  
   92 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE INNER2                                                         
C                                                                               
C CALLED BY SUB. INTEG                                                          
C                                                                               
C**************************************************************                 
C     IT: 1 -- NEUTRON ;      2 -- PROTON                     *                 
C     II: 1 -- P ; 2 -- Q0 ; 3 -- Q1 ; 4 -- Q2 ;              *                 
C                  5 -- P0 ; 6 -- P1 ; 7 -- P2                *                 
C     IS: 1 -- T11B  ; 2 -- T11C  ; 3 -- BT20B ; 4 -- BT20C;  *                 
C         5 -- CT20C ; 6 -- BT02B ; 7 -- BT02C ; 8 -- CT02C   *                 
C**************************************************************                 
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /BCSBCS/ ANP(2),ALAMB(2),DELTA(2),EP(NNBCS,2),                     
     1                HSP(2),QQ0(2),HH0(2),QP0(2),LNP(2)                        
      COMMON /HMILTN/ HM1(NNRANK,4),HM2(NNRANK,NNRANK,6),HM3(NNATA,2),          
     1                HM4(NNATA,NNRANK,2),HM5(NNRANK,NNATA,2),                  
     2                HM6(NNATA,NNATA)                                          
      COMMON /PARAM / AG(2),AGQ(2),DELE,HOLNG2(2),QCL,HFB,                      
     1                AX(3),AXH(3),XQQ(2),GQP(2)                                
      COMMON /TQB1  / NOQ(2),NOQR(2),IOQ(2,MD2),IOQR(2,MD2),                    
     1                NTQ(3),NTQH(2)                                            
      COMMON /TQB2  / ETQ(NNRANT),KTQ(NNRANT),ITQ(NNRANT,3),NRANK               
      COMMON /TQB3  / IOQB(2,MD2),ITQB(2,NNATA,2),NATA(2)                       
      COMMON /TRACE1/ ST02B(2,12),ST20B(2,12),                                  
     1                AAT1(NNRANK,12,2),AAT2(NNRANK,12,2),                      
     2                T1AA(NNRANK,12,2),T2AA(NNRANK,12,2)                       
      COMMON /TRACE2/ AT1A(NNATA,12,2),AT2A(NNATA,12,2)                         
C                                                                               
      WRITE(*,100)                                                              
  100 FORMAT(20('*'),' INNER2')                                                 
C                                                                               
C     ( H AN AN AP AP )                                                         
C                                                                               
      DO 11 I1 = 1,NTQ(1)                                                       
      DO 11 I2 = 1,NTQ(2)                                                       
      SUM = 0.D0                                                                
      DO 21 MU = -2,+2                                                          
      IS = MU+4                                                                 
   21 SUM=SUM+T2AA(I1,IS,1)*T1AA(I2,IS,2)+T2AA(I2,IS,2)*T1AA(I1,IS,1)           
      HM2(I1,I2,1)=-.5*AX(3)*SUM                                                
   11 CONTINUE                                                                  
C                                                                               
C     ( AP AP AN AN H )                                                         
C                                                                               
      DO 12 I1 = 1,NTQ(1)                                                       
      DO 12 I2 = 1,NTQ(2)                                                       
      SUM = 0.D0                                                                
      DO 22 MU = -2,+2                                                          
      IS = MU+4                                                                 
   22 SUM=SUM+AAT1(I2,IS,2)*AAT2(I1,IS,1)+AAT2(I2,IS,2)*AAT1(I1,IS,1)           
      HM2(I1,I2,2)=-.5*AX(3)*SUM                                                
   12 CONTINUE                                                                  
C                                                                               
C     ( AN H AN AP AP )                                                         
C                                                                               
      DO 13 I1 = 1,NATA(1)                                                      
      DO 13 I2 = 1,NTQ(2)                                                       
      SUM = 0.D0                                                                
      DO 23 MU = -2,+2                                                          
      IS = MU+4                                                                 
   23 SUM=SUM+AT2A(I1,IS,1)*T1AA(I2,IS,2)+T2AA(I2,IS,2)*AT1A(I1,IS,1)           
      HM4(I1,I2,1)=-.5*AX(3)*SUM                                                
   13 CONTINUE                                                                  
C                                                                               
C     ( AP H AN AN AP )                                                         
C                                                                               
      DO 14 I1 = 1,NTQ(1)                                                       
      DO 14 I2 = 1,NATA(2)                                                      
      SUM = 0.D0                                                                
      DO 24 MU = -2,+2                                                          
      IS = MU+4                                                                 
   24 SUM=SUM+T2AA(I1,IS,1)*AT1A(I2,IS,2)+AT2A(I2,IS,2)*T1AA(I1,IS,1)           
      HM5(I1,I2,1)=-.5*AX(3)*SUM                                                
   14 CONTINUE                                                                  
C                                                                               
C     ( AP AP AN H AN )                                                         
C                                                                               
      DO 15 I1 = 1,NATA(1)                                                      
      DO 15 I2 = 1,NTQ(2)                                                       
      SUM = 0.D0                                                                
      DO 25 MU = -2,+2                                                          
      IS = MU+4                                                                 
   25 SUM=SUM+AT2A(I1,IS,1)*AAT1(I2,IS,2)+AAT2(I2,IS,2)*AT1A(I1,IS,1)           
      HM4(I1,I2,2)=-.5*AX(3)*SUM                                                
   15 CONTINUE                                                                  
C                                                                               
C     ( AP AN AN H AP )                                                         
C                                                                               
      DO 16 I1 = 1,NTQ(1)                                                       
      DO 16 I2 = 1,NATA(2)                                                      
      SUM = 0.D0                                                                
      DO 26 MU = -2,+2                                                          
      IS = MU+4                                                                 
   26 SUM=SUM+AAT2(I1,IS,1)*AT1A(I2,IS,2)+AT2A(I2,IS,2)*AAT1(I1,IS,1)           
      HM5(I1,I2,2)=-.5*AX(3)*SUM                                                
   16 CONTINUE                                                                  
C                                                                               
C     ( AP AP H AN AN ), ( AN AN H AP AP )                                      
C                                                                               
      DO 18 IT = 1,2                                                            
      DO 17 I1 = 1,NTQ(IT)                                                      
      DO 17 I2 = 1,NTQ(3-IT)                                                    
      SUM = 0.D0                                                                
      DO 27 MU = -2,+2                                                          
      IS = MU+4                                                                 
   27 SUM=SUM+T2AA(I1,IS,IT)*AAT1(I2,IS,3-IT)                                   
     +       +AAT2(I2,IS,3-IT)*T1AA(I1,IS,IT)                                   
      HM2(I2,I1,IT+2)=-.5*AX(3)*SUM                                             
   17 CONTINUE                                                                  
   18 CONTINUE                                                                  
C                                                                               
C     ( AP AN H AN AP )                                                         
C                                                                               
      DO 19 I1 = 1,NATA(1)                                                      
      DO 19 I2 = 1,NATA(2)                                                      
      SUM = 0.D0                                                                
      DO 29 MU = -2,+2                                                          
      IS = MU+4                                                                 
   29 SUM=SUM+AT2A(I1,IS,1)*AT1A(I2,IS,2)+AT2A(I2,IS,2)*AT1A(I1,IS,1)           
      HM6(I1,I2)=-.5*AX(3)*SUM                                                  
   19 CONTINUE                                                                  
C                                                                               
C     ( AN AN H AN AN ), ( AP AP H AP AP )                                      
C                                                                               
      DO 32 IT = 1,2                                                            
      NOQR2 = NOQR(IT)*2                                                        
      DO 31 I1 = 1,NTQ(IT)                                                      
      DO 31 I2 = 1,NTQ(IT)                                                      
      IF (IT.EQ.1) THEN                                                         
        J11=ITQ(I1,1)                                                           
        J12=ITQ(I1,2)                                                           
        J21=ITQ(I2,1)                                                           
        J22=ITQ(I2,2)                                                           
      ELSE                                                                      
        J11=ITQ(NTQ(1)+I1,1)                                                    
        J12=ITQ(NTQ(1)+I1,2)                                                    
        J21=ITQ(NTQ(1)+I2,1)                                                    
        J22=ITQ(NTQ(1)+I2,2)                                                    
      ENDIF                                                                     
      DO 41 II = 0,2                                                            
      HMSS = 0.D0                                                               
      DO 51 MU = -2,+2                                                          
      IF (II.EQ.0.AND.MU.NE.0) GOTO 51                                          
      IS = MU+4+(II-1)*5                                                        
      IF (II.EQ.0) IS = 1                                                       
      HMSS=HMSS+AAT2(I2,IS,IT)*T1AA(I1,IS,IT)                                   
     +         +AAT1(I2,IS,IT)*T2AA(I1,IS,IT)                                   
C                                                                               
      N12 = 0                                                                   
      DO 101 N1 = 1,NOQR2                                                       
      DO 101 N2 = 1,NOQR2                                                       
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J11.AND.IOQB(IT,N2).EQ.J21) GOTO 102                   
  101 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 101'                                  
      STOP                                                                      
  102 N21 = 0                                                                   
      DO 103 N1 = 1,NOQR2                                                       
      DO 103 N2 = 1,NOQR2                                                       
      N21 = N21 + 1                                                             
      IF (IOQB(IT,N1).EQ.J12.AND.IOQB(IT,N2).EQ.J22) GOTO 104                   
  103 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 103'                                  
      STOP                                                                      
  104 CONTINUE                                                                  
      HMSS=HMSS+AT2A(N12,IS,IT)*AT1A(N21,IS,IT)                                 
     +         +AT2A(N21,IS,IT)*AT1A(N12,IS,IT)                                 
C                                                                               
      N12 = 0                                                                   
      DO 105 N1 = 1,NOQR2                                                       
      DO 105 N2 = 1,NOQR2                                                       
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J11.AND.IOQB(IT,N2).EQ.J22) GOTO 106                   
  105 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 105'                                  
      STOP                                                                      
  106 N21 = 0                                                                   
      DO 107 N1 = 1,NOQR2                                                       
      DO 107 N2 = 1,NOQR2                                                       
      N21 = N21 + 1                                                             
      IF (IOQB(IT,N1).EQ.J12.AND.IOQB(IT,N2).EQ.J21) GOTO 108                   
  107 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 107'                                  
      STOP                                                                      
  108 CONTINUE                                                                  
      HMSS=HMSS-AT2A(N12,IS,IT)*AT1A(N21,IS,IT)                                 
     +         -AT2A(N21,IS,IT)*AT1A(N12,IS,IT)                                 
   51 CONTINUE                                                                  
      IF (II.EQ.0) HMSS0=-AG(IT)*HMSS                                           
      IF (II.EQ.1) HMSS1=-.5*AX(IT)*HMSS                                        
      IF (II.EQ.2) HMSS2=-AGQ(IT)*HMSS                                          
   41 CONTINUE                                                                  
C                                                                               
      HM2(I2,I1,IT+4)=HMSS0+HMSS1+HMSS2                                         
   31 CONTINUE                                                                  
   32 CONTINUE                                                                  
C                                                                               
C     ( AN AN H ), ( AP AP H )                                                  
C                                                                               
      DO 34 IT = 1,2                                                            
      AXX = XQQ(IT)                                                             
      DO 33 I1 = 1,NTQ(IT)                                                      
      DO 43 II = 0,2                                                            
      HMSS = 0.D0                                                               
      DO 53 MU = -2,+2                                                          
      IF (II.EQ.0.AND.MU.NE.0) GOTO 53                                          
      IS = MU+4+(II-1)*5                                                        
      IF (II.EQ.0) IS = 1                                                       
      IF (II.EQ.1) THEN                                                         
        HMSS=HMSS+AX(IT)*                                                       
     +  (ST20B(IT,IS)*AAT1(I1,IS,IT)+ST02B(IT,IS)*AAT2(I1,IS,IT))               
     +           +AX(3)*                                                        
     +  (ST02B(3-IT,IS)*AAT2(I1,IS,IT)+ST20B(3-IT,IS)*AAT1(I1,IS,IT))           
      ELSE                                                                      
        HMSS=HMSS+                                                              
     +   ST20B(IT,IS)*AAT1(I1,IS,IT)+ST02B(IT,IS)*AAT2(I1,IS,IT)                
      ENDIF                                                                     
   53 CONTINUE                                                                  
      IF (II.EQ.0) HMSS0 = -AG(IT)*HMSS                                         
      IF (II.EQ.1) HMSS1 = -.5*HMSS                                             
      IF (II.EQ.2) HMSS2 = -AGQ(IT)*HMSS                                        
   43 CONTINUE                                                                  
      HM1(I1,IT)=-AGQ(IT)*QP0(IT)*(AAT2(I1,9,IT)+AAT1(I1,9,IT))                 
     +        +AAT1(I1,12,IT)-AXX*AAT1(I1,4,IT)                                 
     +        -DELTA(IT)*(AAT2(I1,1,IT)+AAT1(I1,1,IT))                          
      HM1(I1,IT)=HM1(I1,IT)+HMSS0+HMSS1+HMSS2                                   
   33 CONTINUE                                                                  
   34 CONTINUE                                                                  
C                                                                               
C     ( AN H AN ), ( AP H AP )                                                  
C                                                                               
      DO 36 IT = 1,2                                                            
      AXX = XQQ(IT)                                                             
      DO 35 I1 = 1,NATA(IT)                                                     
      DO 45 II = 0,2                                                            
      HMSS = 0.D0                                                               
      DO 55 MU = -2,+2                                                          
      IF (II.EQ.0.AND.MU.NE.0) GOTO 55                                          
      IS = MU+4+(II-1)*5                                                        
      IF (II.EQ.0) IS = 1                                                       
      IF (II.EQ.1) THEN                                                         
        HMSS=HMSS+AX(IT)*                                                       
     +  (ST20B(IT,IS)*AT1A(I1,IS,IT)+ST02B(IT,IS)*AT2A(I1,IS,IT))               
     +           +AX(3)*                                                        
     +  (ST02B(3-IT,IS)*AT2A(I1,IS,IT)+ST20B(3-IT,IS)*AT1A(I1,IS,IT))           
      ELSE                                                                      
        HMSS=HMSS+                                                              
     +       ST20B(IT,IS)*AT1A(I1,IS,IT)+ST02B(IT,IS)*AT2A(I1,IS,IT)            
      ENDIF                                                                     
   55 CONTINUE                                                                  
      IF (II.EQ.0) HMSS0 = -AG(IT)*HMSS                                         
      IF (II.EQ.1) HMSS1 = -.5*HMSS                                             
      IF (II.EQ.2) HMSS2 = -AGQ(IT)*HMSS                                        
   45 CONTINUE                                                                  
      HM3(I1,IT)=-AGQ(IT)*QP0(IT)*(AT2A(I1,9,IT)+AT1A(I1,9,IT))                 
     +        +AT1A(I1,12,IT)-AXX*AT1A(I1,4,IT)                                 
     +        -DELTA(IT)*(AT2A(I1,1,IT)+AT1A(I1,1,IT))                          
      HM3(I1,IT)=HM3(I1,IT)+HMSS0+HMSS1+HMSS2                                   
   35 CONTINUE                                                                  
   36 CONTINUE                                                                  
C                                                                               
C     ( H AN AN ), ( H AP AP )                                                  
C                                                                               
      DO 38 IT = 1,2                                                            
      AXX = XQQ(IT)                                                             
      DO 37 I1 = 1,NTQ(IT)                                                      
      DO 47 II = 0,2                                                            
      HMSS = 0.D0                                                               
      DO 57 MU = -2,+2                                                          
      IF (II.EQ.0.AND.MU.NE.0) GOTO 57                                          
      IS = MU+4+(II-1)*5                                                        
      IF (II.EQ.0) IS = 1                                                       
      IF (II.EQ.1) THEN                                                         
        HMSS=HMSS+AX(IT)*                                                       
     +  (ST20B(IT,IS)*T1AA(I1,IS,IT)+ST02B(IT,IS)*T2AA(I1,IS,IT))               
     +           +AX(3)*                                                        
     +  (ST02B(3-IT,IS)*T2AA(I1,IS,IT)+ST20B(3-IT,IS)*T1AA(I1,IS,IT))           
      ELSE                                                                      
        HMSS=HMSS+                                                              
     +       ST20B(IT,IS)*T1AA(I1,IS,IT)+ST02B(IT,IS)*T2AA(I1,IS,IT)            
      ENDIF                                                                     
   57 CONTINUE                                                                  
      IF (II.EQ.0) HMSS0 = -AG(IT)*HMSS                                         
      IF (II.EQ.1) HMSS1 = -.5*HMSS                                             
      IF (II.EQ.2) HMSS2 = -AGQ(IT)*HMSS                                        
   47 CONTINUE                                                                  
      HM1(I1,IT+2)=-AGQ(IT)*QP0(IT)*(T2AA(I1,9,IT)+T1AA(I1,9,IT))               
     +             +T1AA(I1,12,IT)-AXX*T1AA(I1,4,IT)                            
     +             -DELTA(IT)*(T2AA(I1,1,IT)+T1AA(I1,1,IT))                     
      HM1(I1,IT+2)=HM1(I1,IT+2)+HMSS0+HMSS1+HMSS2                               
   37 CONTINUE                                                                  
   38 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE ELEM00(IBET,JBET,HR00)                                         
C                                                                               
C CALLED BY SUB. INTEG                                                          
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
      PARAMETER (NNDEH=NNDEG/2)                                                 
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /HMILTN/ HM1(NNRANK,4),HM2(NNRANK,NNRANK,6),HM3(NNATA,2),          
     1                HM4(NNATA,NNRANK,2),HM5(NNRANK,NNATA,2),                  
     2                HM6(NNATA,NNATA)                                          
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /NORHAM/ WK1(NNRANT,NNRANT,2,NNDEH),                               
     1                WK2(NNRANT,NNRANT,2,NNDEH)                                
      COMMON /TINQAS/ AA(MD2,MD2,2),BB(MD2,MD2,2),CC(MD2,MD2,2)                 
      COMMON /TQB1  / NOQ(2),NOQR(2),IOQ(2,MD2),IOQR(2,MD2),                    
     1                NTQ(3),NTQH(2)                                            
      COMMON /TQB2  / ETQ(NNRANT),KTQ(NNRANT),ITQ(NNRANT,3),NRANK               
      COMMON /TQB3  / IOQB(2,MD2),ITQB(2,NNATA,2),NATA(2)                       
C                                                                               
      DIMENSION C4(2),NP(2,2)                                                   
      DIMENSION HL2A(NNRANK,NNRANK)                                             
C                                                                               
      WRITE(*,100)                                                              
  100 FORMAT(20('*'),' ELEM00')                                                 
C                                                                               
      WK1(NRANK,NRANK,JBET,IBET)=1.D0                                           
      WK2(NRANK,NRANK,JBET,IBET)=HR00                                           
C                                                                               
C     <AN AN H>, <H AN AN>, <AP AP H>, <H AP AP>                                
C                                                                               
      DO 10 IT=1,2                                                              
C                                                                               
      DO 1 I1 = 1,NTQ(IT)                                                       
      IF (IT.EQ.1) THEN                                                         
        J1=ITQ(I1,1)                                                            
        J2=ITQ(I1,2)                                                            
        WK1(I1,NRANK,JBET,IBET)=BB(J2,J1,IT)                                    
        WK2(I1,NRANK,JBET,IBET)=BB(J2,J1,IT)*HR00+HM1(I1,IT)                    
      ELSE                                                                      
        J1=ITQ(NTQ(1)+I1,1)                                                     
        J2=ITQ(NTQ(1)+I1,2)                                                     
        WK1(NTQ(1)+I1,NRANK,JBET,IBET)=BB(J2,J1,IT)                             
        WK2(NTQ(1)+I1,NRANK,JBET,IBET)=BB(J2,J1,IT)*HR00+HM1(I1,IT)             
      ENDIF                                                                     
    1 CONTINUE                                                                  
C                                                                               
      DO 2 I1 = 1,NTQ(IT)                                                       
      IF (IT.EQ.1) THEN                                                         
        J1=ITQ(I1,1)                                                            
        J2=ITQ(I1,2)                                                            
        WK1(NRANK,I1,JBET,IBET)=AA(J1,J2,IT)                                    
        WK2(NRANK,I1,JBET,IBET)=AA(J1,J2,IT)*HR00+HM1(I1,IT+2)                  
      ELSE                                                                      
        J1=ITQ(NTQ(1)+I1,1)                                                     
        J2=ITQ(NTQ(1)+I1,2)                                                     
        WK1(NRANK,NTQ(1)+I1,JBET,IBET)=AA(J1,J2,IT)                             
        WK2(NRANK,NTQ(1)+I1,JBET,IBET)=AA(J1,J2,IT)*HR00+HM1(I1,IT+2)           
      ENDIF                                                                     
    2 CONTINUE                                                                  
   10 CONTINUE                                                                  
C                                                                               
C     <AP AP AN AN H>                                                           
C                                                                               
      DO 3 I1 = 1,NTQ(3)                                                        
      J1=ITQ(NTQ(1)+NTQ(2)+I1,1)                                                
      J2=ITQ(NTQ(1)+NTQ(2)+I1,2)                                                
      J11=ITQ(J1,1)                                                             
      J12=ITQ(J1,2)                                                             
      J21=ITQ(NTQ(1)+J2,1)                                                      
      J22=ITQ(NTQ(1)+J2,2)                                                      
      WK1(NTQ(1)+NTQ(2)+I1,NRANK,JBET,IBET)=                                    
     +          BB(J22,J21,2)*BB(J12,J11,1)                                     
      WK2(NTQ(1)+NTQ(2)+I1,NRANK,JBET,IBET)=HM2(J1,J2,2)                        
     1   + BB(J22,J21,2)*HM1(J1,1) + BB(J12,J11,1)*HM1(J2,2)                    
     2   + BB(J22,J21,2)*BB(J12,J11,1)*HR00                                     
    3 CONTINUE                                                                  
C                                                                               
C     <H AN AN AP AP>                                                           
C                                                                               
      DO 4 I1 = 1,NTQ(3)                                                        
      J1=ITQ(NTQ(1)+NTQ(2)+I1,1)                                                
      J2=ITQ(NTQ(1)+NTQ(2)+I1,2)                                                
      J11=ITQ(J1,1)                                                             
      J12=ITQ(J1,2)                                                             
      J21=ITQ(NTQ(1)+J2,1)                                                      
      J22=ITQ(NTQ(1)+J2,2)                                                      
      WK1(NRANK,NTQ(1)+NTQ(2)+I1,JBET,IBET)=                                    
     1          AA(J11,J12,1)*AA(J21,J22,2)                                     
      WK2(NRANK,NTQ(1)+NTQ(2)+I1,JBET,IBET)=HM2(J1,J2,1)                        
     1   + AA(J11,J12,1)*HM1(J2,4) + AA(J21,J22,2)*HM1(J1,3)                    
     2   + AA(J11,J12,1)*AA(J21,J22,2)*HR00                                     
    4 CONTINUE                                                                  
C                                                                               
C     <AN AN H AN AN>, <AP AP H AP AP>                                          
C                                                                               
      DO 12 IT=1,2                                                              
      DO 5 I1 = 1,NTQ(IT)                                                       
      DO 5 I2 = 1,NTQ(IT)                                                       
      IF (IT.EQ.1) THEN                                                         
        J11=ITQ(I1,1)                                                           
        J12=ITQ(I1,2)                                                           
        J21=ITQ(I2,1)                                                           
        J22=ITQ(I2,2)                                                           
      ELSE                                                                      
        J11=ITQ(NTQ(1)+I1,1)                                                    
        J12=ITQ(NTQ(1)+I1,2)                                                    
        J21=ITQ(NTQ(1)+I2,1)                                                    
        J22=ITQ(NTQ(1)+I2,2)                                                    
      ENDIF                                                                     
      ANORM=BB(J22,J21,IT)*AA(J11,J12,IT)                                       
     1     +CC(J22,J12,IT)*CC(J21,J11,IT)                                       
     2     -CC(J22,J11,IT)*CC(J21,J12,IT)                                       
      HL2A(I2,I1)=BB(J22,J21,IT)*HM1(I1,IT+2)                                   
     1           +AA(J11,J12,IT)*HM1(I2,IT)                                     
     2           +ANORM*HR00+HM2(I2,I1,IT+4)                                    
C                                                                               
      N12 = 0                                                                   
      NOQR2 = NOQR(IT)*2                                                        
      DO 25 N1 = 1,NOQR2                                                        
      DO 25 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J12.AND.IOQB(IT,N2).EQ.J21) GOTO 26                    
   25 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 25'                                   
      WRITE(*,*) ' STOP: FEHLER IN ELEM00 25'                                   
      STOP                                                                      
   26 CONTINUE                                                                  
      HL2A(I2,I1)=HL2A(I2,I1)-CC(J22,J11,IT)*HM3(N12,IT)                        
C                                                                               
      N12 = 0                                                                   
      DO 35 N1 = 1,NOQR2                                                        
      DO 35 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J11.AND.IOQB(IT,N2).EQ.J21) GOTO 36                    
   35 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 35'                                   
      STOP                                                                      
   36 CONTINUE                                                                  
      HL2A(I2,I1)=HL2A(I2,I1)+CC(J22,J12,IT)*HM3(N12,IT)                        
C                                                                               
      N12 = 0                                                                   
      DO 45 N1 = 1,NOQR2                                                        
      DO 45 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J12.AND.IOQB(IT,N2).EQ.J22) GOTO 46                    
   45 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 45'                                   
      STOP                                                                      
   46 CONTINUE                                                                  
      HL2A(I2,I1)=HL2A(I2,I1)+CC(J21,J11,IT)*HM3(N12,IT)                        
C                                                                               
      N12 = 0                                                                   
      DO 55 N1 = 1,NOQR2                                                        
      DO 55 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J11.AND.IOQB(IT,N2).EQ.J22) GOTO 56                    
   55 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 55'                                   
      STOP                                                                      
   56 CONTINUE                                                                  
      HL2A(I2,I1)=HL2A(I2,I1)-CC(J21,J12,IT)*HM3(N12,IT)                        
C                                                                               
      IF (IT.EQ.1) THEN                                                         
        WK1(I2,I1,JBET,IBET)=ANORM                                              
        WK2(I2,I1,JBET,IBET)=HL2A(I2,I1)                                        
      ELSE                                                                      
        WK1(I2+NTQ(1),I1+NTQ(1),JBET,IBET)=ANORM                                
        WK2(I2+NTQ(1),I1+NTQ(1),JBET,IBET)=HL2A(I2,I1)                          
      ENDIF                                                                     
    5 CONTINUE                                                                  
   12 CONTINUE                                                                  
C                                                                               
C     <AN AN H AP AP>                                                           
C                                                                               
      DO 6 I1 = 1,NTQ(1)                                                        
      DO 6 I2 = 1,NTQ(2)                                                        
      J11=ITQ(I2+NTQ(1),1)                                                      
      J12=ITQ(I2+NTQ(1),2)                                                      
      J21=ITQ(I1,1)                                                             
      J22=ITQ(I1,2)                                                             
      WK1(I1,I2+NTQ(1),JBET,IBET)=BB(J22,J21,1)*AA(J11,J12,2)                   
      WK2(I1,I2+NTQ(1),JBET,IBET)=HM2(I1,I2,4)                                  
     1   + BB(J22,J21,1)*HM1(I2,4) + AA(J11,J12,2)*HM1(I1,1)                    
     2   + BB(J22,J21,1)*AA(J11,J12,2)*HR00                                     
    6 CONTINUE                                                                  
C                                                                               
C     <AP AP H AN AN>                                                           
C                                                                               
      DO 7 I1 = 1,NTQ(1)                                                        
      DO 7 I2 = 1,NTQ(2)                                                        
      J11=ITQ(I1,1)                                                             
      J12=ITQ(I1,2)                                                             
      J21=ITQ(I2+NTQ(1),1)                                                      
      J22=ITQ(I2+NTQ(1),2)                                                      
      WK1(I2+NTQ(1),I1,JBET,IBET)=BB(J22,J21,2)*AA(J11,J12,1)                   
      WK2(I2+NTQ(1),I1,JBET,IBET)=HM2(I2,I1,3)                                  
     1   + BB(J22,J21,2)*HM1(I1,3) + AA(J11,J12,1)*HM1(I2,2)                    
     2   + BB(J22,J21,2)*AA(J11,J12,1)*HR00                                     
    7 CONTINUE                                                                  
C                                                                               
C     <AN AN H AN AN AP AP>, <AP AP H AN AN AP AP>                              
C                                                                               
      DO 102 IT = 1,2                                                           
      DO 8 I1 = 1,NTQ(IT)                                                       
      DO 8 I2 = 1,NTQ(3)                                                        
      J11=ITQ(NTQ(1)+NTQ(2)+I2,1)                                               
      J12=ITQ(NTQ(1)+NTQ(2)+I2,2)                                               
      NP(1,1)=ITQ(J11,1)                                                        
      NP(2,1)=ITQ(J11,2)                                                        
      NP(1,2)=ITQ(J12+NTQ(1),1)                                                 
      NP(2,2)=ITQ(J12+NTQ(1),2)                                                 
      IF (IT.EQ.1) THEN                                                         
        J21=ITQ(I1,1)                                                           
        J22=ITQ(I1,2)                                                           
      ELSE                                                                      
        J21=ITQ(NTQ(1)+I1,1)                                                    
        J22=ITQ(NTQ(1)+I1,2)                                                    
      ENDIF                                                                     
C                                                                               
      C4(IT) = BB(J22,J21,IT)*AA(NP(1,IT),NP(2,IT),IT)                          
     1       + CC(J22,NP(2,IT),IT)*CC(J21,NP(1,IT),IT)                          
     2       - CC(J22,NP(1,IT),IT)*CC(J21,NP(2,IT),IT)                          
      ANORM = C4(IT)*AA(NP(1,3-IT),NP(2,3-IT),3-IT)                             
C                                                                               
      SUM=ANORM*HR00                                                            
      SUM=SUM+BB(J22,J21,IT)*HM2(J11,J12,1)                                     
      SUM=SUM+AA(NP(1,1),NP(2,1),1)*AA(NP(1,2),NP(2,2),2)*HM1(I1,IT)            
C                                                                               
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM + AA(NP(1,2),NP(2,2),2)*HM2(I1,J11,IT+4)                        
     1          + AA(NP(1,1),NP(2,1),1)*HM2(I1,J12,5-IT)                        
     2          + C4(1)*HM1(J12,5-IT)                                           
     3          + BB(J22,J21,1)*AA(NP(1,2),NP(2,2),2)*HM1(J11,IT+2)             
      ELSE                                                                      
        SUM=SUM + AA(NP(1,1),NP(2,1),1)*HM2(I1,J12,IT+4)                        
     1          + AA(NP(1,2),NP(2,2),2)*HM2(I1,J11,5-IT)                        
     2          + C4(2)*HM1(J11,5-IT)                                           
     3          + BB(J22,J21,2)*AA(NP(1,1),NP(2,1),1)*HM1(J12,IT+2)             
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      NOQR2 = NOQR(IT)*2                                                        
      DO 27 N1 = 1,NOQR2                                                        
      DO 27 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(2,IT).AND.IOQB(IT,N2).EQ.J21) GOTO 28               
   27 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 27'                                   
      STOP                                                                      
   28 SUM=SUM                                                                   
     1 -CC(J22,NP(1,IT),IT)*AA(NP(1,3-IT),NP(2,3-IT),3-IT)*HM3(N12,IT)          
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM-CC(J22,NP(1,1),1)*HM4(N12,J12,1)                                
      ELSE                                                                      
        SUM=SUM-CC(J22,NP(1,2),2)*HM5(J11,N12,1)                                
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      DO 37 N1 = 1,NOQR2                                                        
      DO 37 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(1,IT).AND.IOQB(IT,N2).EQ.J21) GOTO 38               
   37 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 37'                                   
      STOP                                                                      
   38 SUM=SUM                                                                   
     1 +CC(J22,NP(2,IT),IT)*AA(NP(1,3-IT),NP(2,3-IT),3-IT)*HM3(N12,IT)          
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM+CC(J22,NP(2,1),1)*HM4(N12,J12,1)                                
      ELSE                                                                      
        SUM=SUM+CC(J22,NP(2,2),2)*HM5(J11,N12,1)                                
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      DO 47 N1 = 1,NOQR2                                                        
      DO 47 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(2,IT).AND.IOQB(IT,N2).EQ.J22) GOTO 48               
   47 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 47'                                   
      STOP                                                                      
   48 SUM=SUM                                                                   
     1 +CC(J21,NP(1,IT),IT)*AA(NP(1,3-IT),NP(2,3-IT),3-IT)*HM3(N12,IT)          
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM+CC(J21,NP(1,1),1)*HM4(N12,J12,1)                                
      ELSE                                                                      
        SUM=SUM+CC(J21,NP(1,2),2)*HM5(J11,N12,1)                                
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      DO 57 N1 = 1,NOQR2                                                        
      DO 57 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(1,IT).AND.IOQB(IT,N2).EQ.J22) GOTO 58               
   57 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 57'                                   
      STOP                                                                      
   58 SUM=SUM                                                                   
     1 -CC(J21,NP(2,IT),IT)*AA(NP(1,3-IT),NP(2,3-IT),3-IT)*HM3(N12,IT)          
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM-CC(J21,NP(2,1),1)*HM4(N12,J12,1)                                
      ELSE                                                                      
        SUM=SUM-CC(J21,NP(2,2),2)*HM5(J11,N12,1)                                
      ENDIF                                                                     
C                                                                               
      IF (IT.EQ.1) THEN                                                         
        WK1(I1,NTQ(1)+NTQ(2)+I2,JBET,IBET)=ANORM                                
        WK2(I1,NTQ(1)+NTQ(2)+I2,JBET,IBET)=SUM                                  
      ELSE                                                                      
        WK1(NTQ(1)+I1,NTQ(1)+NTQ(2)+I2,JBET,IBET)=ANORM                         
        WK2(NTQ(1)+I1,NTQ(1)+NTQ(2)+I2,JBET,IBET)=SUM                           
      ENDIF                                                                     
    8 CONTINUE                                                                  
  102 CONTINUE                                                                  
C                                                                               
C     <AP AP AN AN H AN AN>, <AP AP AN AN H AP AP>                              
C                                                                               
      DO 103 IT = 1,2                                                           
      DO 9 I1 = 1,NTQ(3)                                                        
      DO 9 I2 = 1,NTQ(IT)                                                       
      IF (IT.EQ.1) THEN                                                         
        J11=ITQ(I2,1)                                                           
        J12=ITQ(I2,2)                                                           
      ELSE                                                                      
        J11=ITQ(NTQ(1)+I2,1)                                                    
        J12=ITQ(NTQ(1)+I2,2)                                                    
      ENDIF                                                                     
      J21=ITQ(NTQ(1)+NTQ(2)+I1,1)                                               
      J22=ITQ(NTQ(1)+NTQ(2)+I1,2)                                               
      NP(1,1)=ITQ(J21,1)                                                        
      NP(2,1)=ITQ(J21,2)                                                        
      NP(1,2)=ITQ(J22+NTQ(1),1)                                                 
      NP(2,2)=ITQ(J22+NTQ(1),2)                                                 
C                                                                               
      C4(IT) = BB(NP(2,IT),NP(1,IT),IT)*AA(J11,J12,IT)                          
     1       + CC(NP(2,IT),J12,IT)*CC(NP(1,IT),J11,IT)                          
     2       - CC(NP(2,IT),J11,IT)*CC(NP(1,IT),J12,IT)                          
      ANORM = C4(IT)*BB(NP(2,3-IT),NP(1,3-IT),3-IT)                             
C                                                                               
      SUM=ANORM*HR00                                                            
      SUM=SUM+AA(J11,J12,IT)*HM2(J21,J22,2)                                     
     1  +BB(NP(2,2),NP(1,2),2)*BB(NP(2,1),NP(1,1),1)*HM1(I2,IT+2)               
C                                                                               
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM + BB(NP(2,2),NP(1,2),2)*HM2(J21,I2,IT+4)                        
     1          + BB(NP(2,1),NP(1,1),1)*HM2(J22,I2,IT+2)                        
     2          + C4(1)*HM1(J22,3-IT)                                           
     3          + BB(NP(2,2),NP(1,2),2)*AA(J11,J12,1)*HM1(J21,IT)               
      ELSE                                                                      
        SUM=SUM + BB(NP(2,1),NP(1,1),1)*HM2(J22,I2,IT+4)                        
     1          + BB(NP(2,2),NP(1,2),2)*HM2(J21,I2,IT+2)                        
     2          + C4(2)*HM1(J21,3-IT)                                           
     3          + BB(NP(2,1),NP(1,1),1)*AA(J11,J12,2)*HM1(J22,IT)               
      ENDIF                                                                     
C                                                                               
      NOQR2 = NOQR(IT)*2                                                        
      N12 = 0                                                                   
      DO 29 N1 = 1,NOQR2                                                        
      DO 29 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J11.AND.IOQB(IT,N2).EQ.NP(2,IT)) GOTO 30               
   29 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 29'                                   
      STOP                                                                      
   30 SUM=SUM                                                                   
     1 -BB(NP(2,3-IT),NP(1,3-IT),3-IT)*CC(NP(1,IT),J12,IT)*HM3(N12,IT)          
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM-CC(NP(1,1),J12,1)*HM4(N12,J22,2)                                
      ELSE                                                                      
        SUM=SUM-CC(NP(1,2),J12,2)*HM5(J21,N12,2)                                
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      DO 39 N1 = 1,NOQR2                                                        
      DO 39 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J11.AND.IOQB(IT,N2).EQ.NP(1,IT)) GOTO 40               
   39 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 39'                                   
      STOP                                                                      
   40 SUM=SUM                                                                   
     1 +BB(NP(2,3-IT),NP(1,3-IT),3-IT)*CC(NP(2,IT),J12,IT)*HM3(N12,IT)          
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM+CC(NP(2,1),J12,1)*HM4(N12,J22,2)                                
      ELSE                                                                      
        SUM=SUM+CC(NP(2,2),J12,2)*HM5(J21,N12,2)                                
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      DO 49 N1 = 1,NOQR2                                                        
      DO 49 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J12.AND.IOQB(IT,N2).EQ.NP(2,IT)) GOTO 50               
   49 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 49'                                   
      STOP                                                                      
   50 SUM=SUM                                                                   
     1 +BB(NP(2,3-IT),NP(1,3-IT),3-IT)*CC(NP(1,IT),J11,IT)*HM3(N12,IT)          
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM+CC(NP(1,1),J11,1)*HM4(N12,J22,2)                                
      ELSE                                                                      
        SUM=SUM+CC(NP(1,2),J11,2)*HM5(J21,N12,2)                                
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      DO 59 N1 = 1,NOQR2                                                        
      DO 59 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.J12.AND.IOQB(IT,N2).EQ.NP(1,IT)) GOTO 60               
   59 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN ELEM00 59'                                   
      STOP                                                                      
   60 SUM=SUM                                                                   
     1 -BB(NP(2,3-IT),NP(1,3-IT),3-IT)*CC(NP(2,IT),J11,IT)*HM3(N12,IT)          
      IF (IT.EQ.1) THEN                                                         
        SUM=SUM-CC(NP(2,1),J11,1)*HM4(N12,J22,2)                                
      ELSE                                                                      
        SUM=SUM-CC(NP(2,2),J11,2)*HM5(J21,N12,2)                                
      ENDIF                                                                     
C                                                                               
      IF (IT.EQ.1) THEN                                                         
        WK1(NTQ(1)+NTQ(2)+I1,I2,JBET,IBET)=ANORM                                
        WK2(NTQ(1)+NTQ(2)+I1,I2,JBET,IBET)=SUM                                  
      ELSE                                                                      
        WK1(NTQ(1)+NTQ(2)+I1,NTQ(1)+I2,JBET,IBET)=ANORM                         
        WK2(NTQ(1)+NTQ(2)+I1,NTQ(1)+I2,JBET,IBET)=SUM                           
      ENDIF                                                                     
    9 CONTINUE                                                                  
  103 CONTINUE                                                                  
C                                                                               
C     < AP AP AN AN H AN AN AP AP >                                             
C                                                                               
      DO 70 I1 = 1,NTQ(3)                                                       
      DO 70 I2 = 1,NTQ(3)                                                       
      WK2(NTQ(1)+NTQ(2)+I2,NTQ(1)+NTQ(2)+I1,JBET,IBET) =                        
     1    H2N2P(I1,I2,HR00,ANORM)                                               
      WK1(NTQ(1)+NTQ(2)+I2,NTQ(1)+NTQ(2)+I1,JBET,IBET) = ANORM                  
   70 CONTINUE                                                                  
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE BASE00                                                         
C                                                                               
C CALLED BY MAIN                                                                
C                                                                               
C************************************************************                   
C        QUANTUM NUMBERS OF TWO-QUASIPARTICLE STATES        *                   
C     ITQ(IL,1)    INDEX OF A 2-NEUTRON STATE 'IL'          *                   
C     ITQ(IL,2)    INDEX OF A 2-PROTON  STATE 'IL'          *                   
C     ITQ(IL,3)    1/2 FOR K1 +/- K2; 0 FOR PAIRED 2-QP     *                   
C     KTQ(IL)      K-QUANTUM-NUMBER   K1 +/- K2             *                   
C     ETQ(IL)      QUASIPARTICLE-ENERGY E1 + E2             *                   
C************************************************************                   
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /CONFIG/ MXSH(2),NPSH(3,2),NILS(3,2)                               
      COMMON /INTRUD/ ITINT(2),NNINT(2),MIXINT(2),NLINT(2)                      
      COMMON /MAXMAX/ EMAX(3),KMAX(3),NMAX(3),JMAX                              
      COMMON /TQB1  / NOQ(2),NOQR(2),IOQ(2,MD2),IOQR(2,MD2),                    
     1                NTQ(3),NTQH(2)                                            
      COMMON /TQB2  / ETQ(NNRANT),KTQ(NNRANT),ITQ(NNRANT,3),NRANK               
      COMMON /TQB3  / IOQB(2,MD2),ITQB(2,NNATA,2),NATA(2)                       
      COMMON /TQBM  / ETQM(2,NNRANK),KTQM(2,NNRANK),ITQM(2,NNRANK,3)            
      COMMON /UVUV  / ESP(NMD,6),EQP(NMD,6),U(NMD,6),V(NMD,6),                  
     1                KQN(NMD,6)                                                
C                                                                               
      DIMENSION MI(2),MF(2)                                                     
C                                                                               
      CHARACTER*15 TITLE                                                        
C                                                                               
      NP=0                                                                      
      MP=0                                                                      
      DO 15 J=1,2                                                               
      MS=MXSH(J)                                                                
      DO 15 I=1,MS                                                              
      MP=MP+1                                                                   
   15 IF (NILS(I,J).GT.NP) NP=NILS(I,J)                                         
C                                                                               
      WRITE(1,*)                                                                
      WRITE(1,*) ' SINGLE-PARTICLE BASIS'                                       
      WRITE(1,*) ' ~~~~~~~~~~~~~~~~~~~~~'                                       
      WRITE(1,*)                                                                
C                                                                               
      TITLE='NILSSON LEVELS:'                                                   
      WRITE(1,115) TITLE,MXSH,((NPSH(I,J),I=1,MXSH(J)),J=1,2)                   
      DO 16 I=1,NP                                                              
   16 WRITE(1,116) I,(KQN(I,J),ESP(I,J),J=1,MP)                                 
      WRITE(1,*)                                                                
      WRITE(1,*)                                                                
      TITLE='QUASIPARTICLES:'                                                   
      WRITE(1,115) TITLE,MXSH,((NPSH(I,J),I=1,MXSH(J)),J=1,2)                   
      DO 17 I=1,NP                                                              
   17 WRITE(1,116) I,(KQN(I,J),EQP(I,J),J=1,MP)                                 
      WRITE(1,*)                                                                
C                                                                               
      WRITE(1,*)                                                                
      WRITE(1,*) ' MULTI-QUASIPARTICLE BASIS'                                   
      WRITE(1,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~'                                   
C                                                                               
      DO 3 IT=1,2                                                               
      IF (IT.EQ.1) THEN                                                         
        WRITE(1,111)                                                            
        WRITE(1,101) KMAX(IT),EMAX(IT),NMAX(IT)                                 
      ELSE                                                                      
        WRITE(1,112)                                                            
        WRITE(1,102) KMAX(IT),EMAX(IT),NMAX(IT)                                 
      ENDIF                                                                     
      WRITE(1,110) MIXINT(IT),NNINT(IT),NLINT(IT)                               
C                                                                               
C-----------------------------                                                  
C     SELECT 1QP STATES                                                         
C-----------------------------                                                  
C                                                                               
      I=0                                                                       
      DO 1 I1=1,NLINT(IT)                                                       
      I=I+1                                                                     
      IOQ(IT,I)=I1                                                              
    1 IOQR(IT,I1)=I                                                             
      NOQ(IT)=I                                                                 
C                                                                               
C----------------------------------------------------------                     
C     SELECT 2QP STATES                                                         
C     IT DEPENDENTS ON EMAX(IT),KMAX(IT) AND NMAX(IT).                          
C----------------------------------------------------------                     
C                                                                               
      LOW=1                                                                     
      ELOW=EQP(1,MIXINT(IT))                                                    
      DO 50 I=1,NLINT(IT)                                                       
      IF (EQP(I,MIXINT(IT)).GE.ELOW) GOTO 50                                    
      LOW=I                                                                     
      ELOW=EQP(I,MIXINT(IT))                                                    
   50 CONTINUE                                                                  
      NH=NMAX(IT)/2                                                             
      MI(IT)=MAX(1,LOW-NH/2)                                                    
      MF(IT)=MIN(NLINT(IT),MI(IT)+NH-1)                                         
      IF (MI(IT).EQ.1) MF(IT)=MI(IT)+NH-1                                       
      IF (MF(IT).EQ.NLINT(IT)) MI(IT)=MF(IT)-NH+1                               
C                                                                               
        IL = 0                                                                  
        DO 7 I1=MI(IT),MF(IT)                                                   
        DO 8 I2=I1,MF(IT)                                                       
        IF (IL.GE.NH) GOTO 9                                                    
        IF (I1.EQ.I2) GOTO 8                                                    
        ECHK = EQP(I1,MIXINT(IT))+EQP(I2,MIXINT(IT))                            
        IF (ECHK.GT.EMAX(IT)) GOTO 8                                            
        KCHK = (KQN(I1,MIXINT(IT))+KQN(I2,MIXINT(IT)))/2                        
        IF (IABS(KCHK).GT.KMAX(IT)) GOTO 8                                      
        IL = IL + 1                                                             
        ITQM(IT,IL,1)=I1                                                        
        ITQM(IT,IL,2)=I2                                                        
        ITQM(IT,IL,3)=1                                                         
        KTQM(IT,IL)=KCHK                                                        
        ETQM(IT,IL)=ECHK                                                        
    8   CONTINUE                                                                
    7   CONTINUE                                                                
    9   NTQH(IT)  = IL                                                          
C                                                                               
        DO 11 I1=MI(IT),MF(IT)                                                  
        DO 12 I2=I1,MF(IT)                                                      
        IF (IL.GE.NMAX(IT)) GOTO 10                                             
        ECHK = EQP(I1,MIXINT(IT))+EQP(I2,MIXINT(IT))                            
        IF (ECHK.GT.EMAX(IT)) GOTO 12                                           
        KCHK = (KQN(I1,MIXINT(IT))-KQN(I2,MIXINT(IT)))/2                        
        IF (IABS(KCHK).GT.KMAX(IT)) GOTO 12                                     
        IL = IL + 1                                                             
        ITQM(IT,IL,1)=I1                                                        
        ITQM(IT,IL,2)=I2                                                        
        ITQM(IT,IL,3)=2                                                         
        IF (I1.EQ.I2) ITQM(IT,IL,3)=0                                           
        KTQM(IT,IL)=KCHK                                                        
        ETQM(IT,IL)=ECHK                                                        
   12   CONTINUE                                                                
   11   CONTINUE                                                                
   10   NTQ(IT) = IL                                                            
C                                                                               
        CALL ORD2QP(IT)                                                         
C                                                                               
C--------------------------------------------------------                       
C     FIND SINGLE PARTICLES IN SELECTED 2QP STATES                              
C--------------------------------------------------------                       
C                                                                               
      II=0                                                                      
      DO 20 I1=1,NOQ(IT)                                                        
      DO 21 I2=1,NTQ(IT)                                                        
      IF (IOQ(IT,I1).EQ.ITQM(IT,I2,1).OR.                                       
     1    IOQ(IT,I1).EQ.ITQM(IT,I2,2)) THEN                                     
        II=II+1                                                                 
        IOQ(IT,II)=IOQ(IT,I1)                                                   
        IOQR(IT,IOQ(IT,II))=II                                                  
        GOTO 20                                                                 
      ENDIF                                                                     
   21 CONTINUE                                                                  
   20 CONTINUE                                                                  
      NOQR(IT)=II                                                               
C                                                                               
      IIN=NTQH(IT)+1                                                            
      DO 24 II = IIN,NTQ(IT)                                                    
   24 ITQM(IT,II,2)=ITQM(IT,II,2)+NLINT(IT)                                     
C                                                                               
      NOQR2=NOQR(IT)*2                                                          
      N12 = NOQR2*NOQR2                                                         
      IF (N12.GT.NNATA) THEN                                                    
        WRITE(*,*) 'SET NNATA IN EE_PARA FILE TO',N12                           
        STOP                                                                    
      ENDIF                                                                     
      NATA(IT) = N12                                                            
      N12 = 0                                                                   
      DO 27 I1= 1,NOQR2                                                         
      DO 27 I2= 1,NOQR2                                                         
      N12 = N12 + 1                                                             
      ITQB(IT,N12,1)=I1                                                         
   27 ITQB(IT,N12,2)=I2                                                         
      DO 28 I=1,NOQR(IT)                                                        
      IOQB(IT,I)=IOQ(IT,I)                                                      
   28 IOQB(IT,I+NOQR(IT))=IOQ(IT,I)+NLINT(IT)                                   
C                                                                               
C-----------------------------------------                                      
C     PRINT OUT                                                                 
C-----------------------------------------                                      
C                                                                               
      WRITE(1,104) NOQR(IT)                                                     
      IF (NOQR(IT).EQ.0) GOTO 32                                                
      DO 30 I1=1,NOQR(IT)                                                       
      II=IOQ(IT,I1)                                                             
      UMV = U(II,MIXINT(IT))**2-V(II,MIXINT(IT))**2                             
      UV  = U(II,MIXINT(IT))*V(II,MIXINT(IT))                                   
   30 WRITE(1,105) I1,II,KQN(II,MIXINT(IT)),EQP(II,MIXINT(IT)),                 
     +             UMV,UV,U(II,MIXINT(IT)),V(II,MIXINT(IT))                     
   32 CONTINUE                                                                  
C                                                                               
      WRITE(1,106) NTQ(IT)                                                      
      IF (NTQ(IT).EQ.0) GOTO 33                                                 
      DO 31 II = 1,NTQ(IT)                                                      
      I1 = ITQM(IT,II,1)                                                        
      IF (II.LE.NTQH(IT)) THEN                                                  
        I2 = ITQM(IT,II,2)                                                      
      ELSE                                                                      
        I2 = ITQM(IT,II,2)-NLINT(IT)                                            
      ENDIF                                                                     
      UMV = U(I1,MIXINT(IT))**2 - V(I1,MIXINT(IT))**2                           
     1    + U(I2,MIXINT(IT))**2 - V(I2,MIXINT(IT))**2                           
      UV  = U(I1,MIXINT(IT))*V(I1,MIXINT(IT))                                   
     1    + U(I2,MIXINT(IT))*V(I2,MIXINT(IT))                                   
   31 WRITE(1,107) II,(ITQM(IT,II,J),J=1,3),KTQM(IT,II),ETQM(IT,II),            
     1             UMV,UV                                                       
   33 CONTINUE                                                                  
C                                                                               
    3 CONTINUE                                                                  
C--------------------------------------------                                   
C     REARANGEMENT                                                              
C--------------------------------------------                                   
      II=0                                                                      
      DO 4 IT=1,2                                                               
      DO 4 I=1,NTQ(IT)                                                          
      II=II+1                                                                   
      KTQ(II) = KTQM(IT,I)                                                      
      ETQ(II) = ETQM(IT,I)                                                      
      DO 5 J=1,3                                                                
    5 ITQ(II,J) = ITQM(IT,I,J)                                                  
    4 CONTINUE                                                                  
C--------------------------------------------                                   
C     MAKE 4QP STATES [2N X 2P], T=1 only!?                                     
C--------------------------------------------                                   
      NH1=NTQ(1)                                                                
      NH2=NTQ(2)                                                                
      IL = 0                                                                    
      DO 41 I1 = 1,NH1                                                          
      DO 41 I2 = 1,NH2                                                          
      IF (IL.GE.NMAX(3)) GOTO 42                                                
      ECHK=ETQ(I1)+ETQ(NTQ(1)+I2)                                               
      KCHK=KTQ(I1)+KTQ(NTQ(1)+I2)                                               
      IF (ECHK.GT.EMAX(3)) GOTO 41                                              
      IF (IABS(KCHK).GT.KMAX(3)) GOTO 41                                        
      IL = IL + 1                                                               
      ITQM(1,IL,1)=I1                                                           
      ITQM(1,IL,2)=I2                                                           
      ITQM(1,IL,3)=1                                                            
      KTQM(1,IL)=KCHK                                                           
      ETQM(1,IL)=ECHK                                                           
   41 CONTINUE                                                                  
   42 NTQ(3) = IL                                                               
C                                                                               
      II=NTQ(1)+NTQ(2)                                                          
      DO 44 I=1,NTQ(3)                                                          
      II=II+1                                                                   
      KTQ(II) = KTQM(1,I)                                                       
      ETQ(II) = ETQM(1,I)                                                       
      DO 44 J=1,3                                                               
   44 ITQ(II,J) = ITQM(1,I,J)                                                   
      NRANK=II+1                                                                
      ETQ(NRANK) = 0.                                                           
      KTQ(NRANK) = 0                                                            
      DO 45 J=1,3                                                               
   45 ITQ(NRANK,J) = 0                                                          
C----------------------------------------------                                 
C     0QP [QP VACUUM] STATE STORED AT THE END                                   
C----------------------------------------------                                 
        IF (NRANK.EQ.0) THEN                                                    
         WRITE(*,100)                                                           
         WRITE(1,100)                                                           
  100    FORMAT(' STOP: CONFIG. SPACE EMPTY')                                   
         STOP                                                                   
        ENDIF                                                                   
C                                                                               
      IF (NRANK.GT.NNRANT) THEN                                                 
        WRITE(*,*) 'SET NNRANT IN EE_PARA FILE TO',NRANK                        
        STOP                                                                    
      ENDIF                                                                     
C                                                                               
      WRITE(1,108)                                                              
      WRITE(1,109) NRANK, (NTQ(I),I=1,3)                                        
      WRITE(1,*)                                                                
      WRITE(1,*) ' ----------   THE QP BASE CHOSEN:   ----------'               
      WRITE(1,*)                                                                
      JJ=3                                                                      
      IF (NTQ(3).EQ.0) JJ=2                                                     
      II=0                                                                      
      DO 72 J=1,JJ                                                              
      IF (J.EQ.1)                                                               
     + WRITE(1,*)' BAND  I1   I2    T   KQP   EQP : 2N STATES'                  
      IF (J.EQ.2)                                                               
     + WRITE(1,*)' BAND  I1   I2    T   KQP   EQP : 2P STATES'                  
      IF (J.EQ.3)                                                               
     + WRITE(1,*)' BAND  I1   I2    T   KQP   EQP : 2NX2P STATES'               
      NI=II+1                                                                   
      NF=II+NTQ(J)                                                              
      DO 71 I=NI,NF                                                             
      II=II+1                                                                   
      I1=ITQ(I,1)                                                               
      I2=ITQ(I,2)                                                               
      I3=ITQ(I,3)                                                               
      WRITE(1,114) II,I1,I2,I3,KTQ(I),ETQ(I)                                    
   71 CONTINUE                                                                  
   72 CONTINUE                                                                  
      II=II+1                                                                   
      I1=0                                                                      
      I2=0                                                                      
      I3=0                                                                      
      WRITE(1,*)' BAND  I1   I2    T   KQP   EQP : VACUUM STATE'                
      WRITE(1,114) II,I1,I2,I3,KTQ(II),ETQ(II)                                  
C                                                                               
  101 FORMAT(/5X,'KMAX =',I3,5X,'EMAX =',F5.2,5X,'NMAX =',I3,                   
     1        5X,'[2-QUASINEUTRON]')                                            
  102 FORMAT(/5X,'KMAX =',I3,5X,'EMAX =',F5.2,5X,'NMAX =',I3,                   
     1        5X,'[2-QUASIPROTON]')                                             
  103 FORMAT('NTQH =',2I3,'       NTQ =',2I3,'      NOQ =',2I3)                 
  104 FORMAT(/,' 1-QP BASIS:  NOQ =',I3,' LEVELS',/,30('-'),/,                  
     1      '  I',' IOQ(I)','    K',6X,'E',10X,'U**2-V**2',                     
     2      6X,'UV',10X,'U',11X,'V')                                            
  105 FORMAT(I3,I4,I7,'/2',5F12.5)                                              
  106 FORMAT(/,' 2-QP BASIS:  NTQ =',I3,' PAIRS',/,30('-'),/,                   
     + ' II   I1 I2   T   K',5X,'E',10X,'U**2-V**2',6X,'UV')                    
  107 FORMAT(I3,I5,I3,2I4,3F12.5)                                               
  108 FORMAT(/)                                                                 
  109 FORMAT(' DIMENSION OF THE CONFIG. SPACE:',I4,' (',3I3,' 1 )')             
  110 FORMAT(3X,'MIXINT =',I3,4X,'NNINT =',I3,6X,'NLINT =',I3)                  
  111 FORMAT(/,'  NEUTRON STATES :')                                            
  112 FORMAT(/,'  PROTON  STATES :')                                            
  113 FORMAT(4I10)                                                              
  114 FORMAT(4I5,I6,F7.3)                                                       
  115 FORMAT(1X,'NO. 2*K ENERGY OF ',A15,I2,' NEUTRON- AND',I2,                 
     1      ' PROTON-SHELLS'/1X,'--- --- ------'/4X,6(6X,'N=',I1,3X))           
  116 FORMAT(I3,6(I5,F7.3))                                                     
C                                                                               
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE ORD2QP(IT)                                                     
C                                                                               
      INCLUDE 'EE_PARA'                                                         
C                                                                               
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /TQB1  / NOQ(2),NOQR(2),IOQ(2,MD2),IOQR(2,MD2),                    
     1                NTQ(3),NTQH(2)                                            
      COMMON /TQBM  / ETQM(2,NNRANK),KTQM(2,NNRANK),ITQM(2,NNRANK,3)            
C                                                                               
      INI=1                                                                     
      IFI=NTQH(IT)                                                              
      DO 4 M=1,2                                                                
      DO 3 I=INI,IFI                                                            
      L=I                                                                       
      E=ETQM(IT,L)                                                              
      DO 1 J=I,IFI                                                              
      IF (ETQM(IT,J).GE.E) GOTO 1                                               
      L=J                                                                       
      E=ETQM(IT,L)                                                              
    1 CONTINUE                                                                  
      IF (L.EQ.I) GOTO 3                                                        
      ETQM(IT,L)=ETQM(IT,I)                                                     
      ETQM(IT,I)=E                                                              
      K=KTQM(IT,L)                                                              
      KTQM(IT,L)=KTQM(IT,I)                                                     
      KTQM(IT,I)=K                                                              
      DO 2 N=1,3                                                                
      K=ITQM(IT,L,N)                                                            
      ITQM(IT,L,N)=ITQM(IT,I,N)                                                 
      ITQM(IT,I,N)=K                                                            
    2 CONTINUE                                                                  
    3 CONTINUE                                                                  
      INI=NTQH(IT)+1                                                            
      IFI=NTQ(IT)                                                               
    4 CONTINUE                                                                  
      END                                                                       
C                                                                               
C                                                                               
      DOUBLE PRECISION FUNCTION H2N2P(I1,I2,HR00,ANORM)                         
C                                                                               
C CALLED BY SUB. ELEM00                                                         
C                                                                               
      INCLUDE 'EE_PARA'                                                         
      PARAMETER (MD2=2*NMD)                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /HMILTN/ HM1(NNRANK,4),HM2(NNRANK,NNRANK,6),HM3(NNATA,2),          
     1                HM4(NNATA,NNRANK,2),HM5(NNRANK,NNATA,2),                  
     2                HM6(NNATA,NNATA)                                          
      COMMON /TINQAS/ AA(MD2,MD2,2),BB(MD2,MD2,2),CC(MD2,MD2,2)                 
      COMMON /TQB1  / NOQ(2),NOQR(2),IOQ(2,MD2),IOQR(2,MD2),                    
     1                NTQ(3),NTQH(2)                                            
      COMMON /TQB2  / ETQ(NNRANT),KTQ(NNRANT),ITQ(NNRANT,3),NRANK               
      COMMON /TQB3  / IOQB(2,MD2),ITQB(2,NNATA,2),NATA(2)                       
C                                                                               
      DIMENSION NP(4,2),C4(2),C6(6,2)                                           
C                                                                               
      J11=ITQ(NTQ(1)+NTQ(2)+I1,1)                                               
      J12=ITQ(NTQ(1)+NTQ(2)+I1,2)                                               
      J21=ITQ(NTQ(1)+NTQ(2)+I2,1)                                               
      J22=ITQ(NTQ(1)+NTQ(2)+I2,2)                                               
      NP(1,1)=ITQ(J11,1)                                                        
      NP(2,1)=ITQ(J11,2)                                                        
      NP(3,1)=ITQ(J21,1)                                                        
      NP(4,1)=ITQ(J21,2)                                                        
      NP(1,2)=ITQ(NTQ(1)+J12,1)                                                 
      NP(2,2)=ITQ(NTQ(1)+J12,2)                                                 
      NP(3,2)=ITQ(NTQ(1)+J22,1)                                                 
      NP(4,2)=ITQ(NTQ(1)+J22,2)                                                 
C                                                                               
      SUM=BB(NP(4,2),NP(3,2),2)*BB(NP(4,1),NP(3,1),1)*HM2(J11,J12,1)            
     1   +AA(NP(1,1),NP(2,1),1)*AA(NP(1,2),NP(2,2),2)*HM2(J21,J22,2)            
C                                                                               
      DO 11 IT = 1,2                                                            
      NOQR2 = NOQR(IT)*2                                                        
      N12 = 0                                                                   
      DO 65 N1 = 1,NOQR2                                                        
      DO 65 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(1,IT).AND.IOQB(IT,N2).EQ.NP(3,IT)) GOTO 66          
   65 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 65'                                    
      STOP                                                                      
   66 CONTINUE                                                                  
      IF (IT.EQ.1) THEN                                                         
      SUM=SUM+CC(NP(4,1),NP(2,1),1)*(BB(NP(4,2),NP(3,2),2)                      
     1       *HM4(N12,J12,1)+AA(NP(1,2),NP(2,2),2)*HM4(N12,J22,2))              
      ELSE                                                                      
      SUM=SUM+CC(NP(4,2),NP(2,2),2)*(BB(NP(4,1),NP(3,1),1)                      
     1       *HM5(J11,N12,1)+AA(NP(1,1),NP(2,1),1)*HM5(J21,N12,2))              
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      DO 75 N1 = 1,NOQR2                                                        
      DO 75 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(2,IT).AND.IOQB(IT,N2).EQ.NP(3,IT)) GOTO 76          
   75 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 75'                                    
      STOP                                                                      
   76 CONTINUE                                                                  
      IF (IT.EQ.1) THEN                                                         
      SUM=SUM-CC(NP(4,1),NP(1,1),1)*(BB(NP(4,2),NP(3,2),2)                      
     1       *HM4(N12,J12,1)+AA(NP(1,2),NP(2,2),2)*HM4(N12,J22,2))              
      ELSE                                                                      
      SUM=SUM-CC(NP(4,2),NP(1,2),2)*(BB(NP(4,1),NP(3,1),1)                      
     1       *HM5(J11,N12,1)+AA(NP(1,1),NP(2,1),1)*HM5(J21,N12,2))              
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      DO 85 N1 = 1,NOQR2                                                        
      DO 85 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(2,IT).AND.IOQB(IT,N2).EQ.NP(4,IT)) GOTO 86          
   85 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 85'                                    
      STOP                                                                      
   86 CONTINUE                                                                  
      IF (IT.EQ.1) THEN                                                         
      SUM=SUM+CC(NP(3,1),NP(1,1),1)*(BB(NP(4,2),NP(3,2),2)                      
     1       *HM4(N12,J12,1)+AA(NP(1,2),NP(2,2),2)*HM4(N12,J22,2))              
      ELSE                                                                      
      SUM=SUM+CC(NP(3,2),NP(1,2),2)*(BB(NP(4,1),NP(3,1),1)                      
     1       *HM5(J11,N12,1)+AA(NP(1,1),NP(2,1),1)*HM5(J21,N12,2))              
      ENDIF                                                                     
C                                                                               
      N12 = 0                                                                   
      DO 95 N1 = 1,NOQR2                                                        
      DO 95 N2 = 1,NOQR2                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(1,IT).AND.IOQB(IT,N2).EQ.NP(4,IT)) GOTO 96          
   95 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 95'                                    
      STOP                                                                      
   96 CONTINUE                                                                  
      IF (IT.EQ.1) THEN                                                         
      SUM=SUM-CC(NP(3,1),NP(2,1),1)*(BB(NP(4,2),NP(3,2),2)                      
     1       *HM4(N12,J12,1)+AA(NP(1,2),NP(2,2),2)*HM4(N12,J22,2))              
      ELSE                                                                      
      SUM=SUM-CC(NP(3,2),NP(2,2),2)*(BB(NP(4,1),NP(3,1),1)                      
     1       *HM5(J11,N12,1)+AA(NP(1,1),NP(2,1),1)*HM5(J21,N12,2))              
      ENDIF                                                                     
   11 CONTINUE                                                                  
C                                                                               
      DO 10 IT = 1,2                                                            
      C4(IT) = BB(NP(4,IT),NP(3,IT),IT)*AA(NP(1,IT),NP(2,IT),IT)                
     1        +CC(NP(4,IT),NP(2,IT),IT)*CC(NP(3,IT),NP(1,IT),IT)                
     2        -CC(NP(4,IT),NP(1,IT),IT)*CC(NP(3,IT),NP(2,IT),IT)                
      C6(1,3-IT) = AA(NP(1,3-IT),NP(2,3-IT),3-IT)*C4(IT)                        
      C6(2,3-IT) = BB(NP(4,3-IT),NP(3,3-IT),3-IT)*C4(IT)                        
      C6(3,3-IT) = CC(NP(4,3-IT),NP(2,3-IT),3-IT)*C4(IT)                        
      C6(4,3-IT) = CC(NP(4,3-IT),NP(1,3-IT),3-IT)*C4(IT)                        
      C6(5,3-IT) = CC(NP(3,3-IT),NP(1,3-IT),3-IT)*C4(IT)                        
      C6(6,3-IT) = CC(NP(3,3-IT),NP(2,3-IT),3-IT)*C4(IT)                        
   10 CONTINUE                                                                  
C                                                                               
      SUM=SUM+C6(1,1)*HM1(J21,1)+C6(2,1)*HM1(J11,3)                             
     1       +C6(1,2)*HM1(J22,2)+C6(2,2)*HM1(J12,4)                             
C                                                                               
      DO 1 IT = 1,2                                                             
      NOQR2 = NOQR(IT)*2                                                        
      N12 = 0                                                                   
      DO 2 N1 = 1,NOQR2                                                         
      DO 2 N2 = 1,NOQR2                                                         
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(1,IT).AND.IOQB(IT,N2).EQ.NP(3,IT)) GOTO 3           
    2 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 2'                                     
      STOP                                                                      
    3 CONTINUE                                                                  
      SUM=SUM+C6(3,IT)*HM3(N12,IT)                                              
C                                                                               
      N12 = 0                                                                   
      DO 4 N1 = 1,NOQR2                                                         
      DO 4 N2 = 1,NOQR2                                                         
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(2,IT).AND.IOQB(IT,N2).EQ.NP(3,IT)) GOTO 5           
    4 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 4'                                     
      STOP                                                                      
    5 CONTINUE                                                                  
      SUM=SUM-C6(4,IT)*HM3(N12,IT)                                              
C                                                                               
      N12 = 0                                                                   
      DO 6 N1 = 1,NOQR2                                                         
      DO 6 N2 = 1,NOQR2                                                         
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(2,IT).AND.IOQB(IT,N2).EQ.NP(4,IT)) GOTO 7           
    6 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 6'                                     
      STOP                                                                      
    7 CONTINUE                                                                  
      SUM=SUM+C6(5,IT)*HM3(N12,IT)                                              
C                                                                               
      N12 = 0                                                                   
      DO 8 N1 = 1,NOQR2                                                         
      DO 8 N2 = 1,NOQR2                                                         
      N12 = N12 + 1                                                             
      IF (IOQB(IT,N1).EQ.NP(1,IT).AND.IOQB(IT,N2).EQ.NP(4,IT)) GOTO 9           
    8 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 8'                                     
      STOP                                                                      
    9 CONTINUE                                                                  
      SUM=SUM-C6(6,IT)*HM3(N12,IT)                                              
    1 CONTINUE                                                                  
C                                                                               
      ANORM=C4(1)*C4(2)                                                         
      SUM=SUM+ANORM*HR00                                                        
C                                                                               
      SUM=SUM+C4(1)*HM2(J22,J12,6)+C4(2)*HM2(J21,J11,5)                         
     1   +BB(NP(4,1),NP(3,1),1)*AA(NP(1,2),NP(2,2),2)*HM2(J22,J11,3)            
     2   +BB(NP(4,2),NP(3,2),2)*AA(NP(1,1),NP(2,1),1)*HM2(J21,J12,4)            
C                                                                               
      NOQR1 = NOQR(1)*2                                                         
      NOQR2 = NOQR(2)*2                                                         
      N12 = 0                                                                   
      DO 13 N1 = 1,NOQR1                                                        
      DO 13 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(2,1).AND.IOQB(1,N2).EQ.NP(4,1)) GOTO 14              
   13 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 13'                                    
      STOP                                                                      
   14 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 15 N1 = 1,NOQR2                                                        
      DO 15 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(2,2).AND.IOQB(2,N2).EQ.NP(4,2)) GOTO 16              
   15 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 15'                                    
      STOP                                                                      
   16 CONTINUE                                                                  
      SUM=SUM+CC(NP(3,1),NP(1,1),1)*CC(NP(3,2),NP(1,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 23 N1 = 1,NOQR1                                                        
      DO 23 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(2,1).AND.IOQB(1,N2).EQ.NP(4,1)) GOTO 24              
   23 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 23'                                    
      STOP                                                                      
   24 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 25 N1 = 1,NOQR2                                                        
      DO 25 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(1,2).AND.IOQB(2,N2).EQ.NP(4,2)) GOTO 26              
   25 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 25'                                    
      STOP                                                                      
   26 CONTINUE                                                                  
      SUM=SUM-CC(NP(3,1),NP(1,1),1)*CC(NP(3,2),NP(2,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 33 N1 = 1,NOQR1                                                        
      DO 33 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(2,1).AND.IOQB(1,N2).EQ.NP(4,1)) GOTO 34              
   33 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 33'                                    
      STOP                                                                      
   34 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 35 N1 = 1,NOQR2                                                        
      DO 35 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(1,2).AND.IOQB(2,N2).EQ.NP(3,2)) GOTO 36              
   35 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 35'                                    
      STOP                                                                      
   36 CONTINUE                                                                  
      SUM=SUM+CC(NP(3,1),NP(1,1),1)*CC(NP(4,2),NP(2,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 43 N1 = 1,NOQR1                                                        
      DO 43 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(2,1).AND.IOQB(1,N2).EQ.NP(4,1)) GOTO 44              
   43 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 43'                                    
      STOP                                                                      
   44 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 45 N1 = 1,NOQR2                                                        
      DO 45 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(2,2).AND.IOQB(2,N2).EQ.NP(3,2)) GOTO 46              
   45 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 45'                                    
      STOP                                                                      
   46 CONTINUE                                                                  
      SUM=SUM-CC(NP(3,1),NP(1,1),1)*CC(NP(4,2),NP(1,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 17 N1 = 1,NOQR1                                                        
      DO 17 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(1,1).AND.IOQB(1,N2).EQ.NP(4,1)) GOTO 18              
   17 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 17'                                    
      STOP                                                                      
   18 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 19 N1 = 1,NOQR2                                                        
      DO 19 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(2,2).AND.IOQB(2,N2).EQ.NP(4,2)) GOTO 20              
   19 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 19'                                    
      STOP                                                                      
   20 CONTINUE                                                                  
      SUM=SUM-CC(NP(3,1),NP(2,1),1)*CC(NP(3,2),NP(1,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 27 N1 = 1,NOQR1                                                        
      DO 27 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(1,1).AND.IOQB(1,N2).EQ.NP(4,1)) GOTO 28              
   27 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 27'                                    
      STOP                                                                      
   28 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 29 N1 = 1,NOQR2                                                        
      DO 29 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(1,2).AND.IOQB(2,N2).EQ.NP(4,2)) GOTO 30              
   29 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 29'                                    
      STOP                                                                      
   30 CONTINUE                                                                  
      SUM=SUM+CC(NP(3,1),NP(2,1),1)*CC(NP(3,2),NP(2,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 37 N1 = 1,NOQR1                                                        
      DO 37 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(1,1).AND.IOQB(1,N2).EQ.NP(4,1)) GOTO 40              
   37 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 37'                                    
      STOP                                                                      
   40 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 38 N1 = 1,NOQR2                                                        
      DO 38 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(1,2).AND.IOQB(2,N2).EQ.NP(3,2)) GOTO 39              
   38 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 38'                                    
      STOP                                                                      
   39 CONTINUE                                                                  
      SUM=SUM-CC(NP(3,1),NP(2,1),1)*CC(NP(4,2),NP(2,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 47 N1 = 1,NOQR1                                                        
      DO 47 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(1,1).AND.IOQB(1,N2).EQ.NP(4,1)) GOTO 48              
   47 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 47'                                    
      STOP                                                                      
   48 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 49 N1 = 1,NOQR2                                                        
      DO 49 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(2,2).AND.IOQB(2,N2).EQ.NP(3,2)) GOTO 50              
   49 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 49'                                    
      STOP                                                                      
   50 CONTINUE                                                                  
      SUM=SUM+CC(NP(3,1),NP(2,1),1)*CC(NP(4,2),NP(1,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 53 N1 = 1,NOQR1                                                        
      DO 53 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(1,1).AND.IOQB(1,N2).EQ.NP(3,1)) GOTO 54              
   53 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 53'                                    
      STOP                                                                      
   54 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 55 N1 = 1,NOQR2                                                        
      DO 55 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(2,2).AND.IOQB(2,N2).EQ.NP(4,2)) GOTO 56              
   55 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 55'                                    
      STOP                                                                      
   56 CONTINUE                                                                  
      SUM=SUM+CC(NP(4,1),NP(2,1),1)*CC(NP(3,2),NP(1,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 63 N1 = 1,NOQR1                                                        
      DO 63 N2 = 1,NOQR1                                                        
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(1,1).AND.IOQB(1,N2).EQ.NP(3,1)) GOTO 64              
   63 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 63'                                    
      STOP                                                                      
   64 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 67 N1 = 1,NOQR2                                                        
      DO 67 N2 = 1,NOQR2                                                        
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(1,2).AND.IOQB(2,N2).EQ.NP(4,2)) GOTO 68              
   67 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 67'                                    
      STOP                                                                      
   68 CONTINUE                                                                  
      SUM=SUM-CC(NP(4,1),NP(2,1),1)*CC(NP(3,2),NP(2,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 133 N1 = 1,NOQR1                                                       
      DO 133 N2 = 1,NOQR1                                                       
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(1,1).AND.IOQB(1,N2).EQ.NP(3,1)) GOTO 134             
  133 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 133'                                   
      STOP                                                                      
  134 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 135 N1 = 1,NOQR2                                                       
      DO 135 N2 = 1,NOQR2                                                       
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(1,2).AND.IOQB(2,N2).EQ.NP(3,2)) GOTO 136             
  135 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 135'                                   
      STOP                                                                      
  136 CONTINUE                                                                  
      SUM=SUM+CC(NP(4,1),NP(2,1),1)*CC(NP(4,2),NP(2,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 143 N1 = 1,NOQR1                                                       
      DO 143 N2 = 1,NOQR1                                                       
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(1,1).AND.IOQB(1,N2).EQ.NP(3,1)) GOTO 144             
  143 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 143'                                   
      STOP                                                                      
  144 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 145 N1 = 1,NOQR2                                                       
      DO 145 N2 = 1,NOQR2                                                       
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(2,2).AND.IOQB(2,N2).EQ.NP(3,2)) GOTO 146             
  145 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 145'                                   
      STOP                                                                      
  146 CONTINUE                                                                  
      SUM=SUM-CC(NP(4,1),NP(2,1),1)*CC(NP(4,2),NP(1,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 117 N1 = 1,NOQR1                                                       
      DO 117 N2 = 1,NOQR1                                                       
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(2,1).AND.IOQB(1,N2).EQ.NP(3,1)) GOTO 118             
  117 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 117'                                   
      STOP                                                                      
  118 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 119 N1 = 1,NOQR2                                                       
      DO 119 N2 = 1,NOQR2                                                       
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(2,2).AND.IOQB(2,N2).EQ.NP(4,2)) GOTO 120             
  119 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 119'                                   
      STOP                                                                      
  120 CONTINUE                                                                  
      SUM=SUM-CC(NP(4,1),NP(1,1),1)*CC(NP(3,2),NP(1,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 127 N1 = 1,NOQR1                                                       
      DO 127 N2 = 1,NOQR1                                                       
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(2,1).AND.IOQB(1,N2).EQ.NP(3,1)) GOTO 128             
  127 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 127'                                   
      STOP                                                                      
  128 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 129 N1 = 1,NOQR2                                                       
      DO 129 N2 = 1,NOQR2                                                       
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(1,2).AND.IOQB(2,N2).EQ.NP(4,2)) GOTO 130             
  129 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 129'                                   
      STOP                                                                      
  130 CONTINUE                                                                  
      SUM=SUM+CC(NP(4,1),NP(1,1),1)*CC(NP(3,2),NP(2,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 137 N1 = 1,NOQR1                                                       
      DO 137 N2 = 1,NOQR1                                                       
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(2,1).AND.IOQB(1,N2).EQ.NP(3,1)) GOTO 140             
  137 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 137'                                   
      STOP                                                                      
  140 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 138 N1 = 1,NOQR2                                                       
      DO 138 N2 = 1,NOQR2                                                       
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(1,2).AND.IOQB(2,N2).EQ.NP(3,2)) GOTO 139             
  138 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 138'                                   
      STOP                                                                      
  139 CONTINUE                                                                  
      SUM=SUM-CC(NP(4,1),NP(1,1),1)*CC(NP(4,2),NP(2,2),2)*HM6(N12,N34)          
C                                                                               
      N12 = 0                                                                   
      DO 147 N1 = 1,NOQR1                                                       
      DO 147 N2 = 1,NOQR1                                                       
      N12 = N12 + 1                                                             
      IF (IOQB(1,N1).EQ.NP(2,1).AND.IOQB(1,N2).EQ.NP(3,1)) GOTO 148             
  147 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 147'                                   
      STOP                                                                      
  148 CONTINUE                                                                  
      N34 = 0                                                                   
      DO 149 N1 = 1,NOQR2                                                       
      DO 149 N2 = 1,NOQR2                                                       
      N34 = N34 + 1                                                             
      IF (IOQB(2,N1).EQ.NP(2,2).AND.IOQB(2,N2).EQ.NP(3,2)) GOTO 150             
  149 CONTINUE                                                                  
      WRITE(1,*) ' STOP: FEHLER IN H2N2P 149'                                   
      STOP                                                                      
  150 CONTINUE                                                                  
      SUM=SUM+CC(NP(4,1),NP(1,1),1)*CC(NP(4,2),NP(1,2),2)*HM6(N12,N34)          
C                                                                               
      H2N2P=SUM                                                                 
      END                                                                       
C                                                                               
C                                                                               
      DOUBLE PRECISION FUNCTION CLEBG(AJ,AM,BJ,BM,CJ,CM)                        
C                                                                               
      PARAMETER (NIE=127)                                                       
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      INTEGER I(11)                                                             
C                                                                               
      COMMON /GFVLOG/ G(NIE),LOGTAB                                             
C                                                                               
C     COMPUTE LOGARITHMS OF FACTORIALS                                          
c      G(1)=0.0D00                                                              
c      DO 20 J=1,LOGTAB                                                         
c      X=J                                                                      
c      G(J+1)=G(J)+DLOG(X)                                                      
c   20 CONTINUE                                                                 
C                                                                               
   30 CLEB=0                                                                    
C     CONVERT THE ARGUMENT TO INTEGER                                           
      J1=NINT(2.*AJ)                                                            
      J2=NINT(2.*BJ)                                                            
      J3=NINT(2.*CJ)                                                            
      M1=NINT(2.*AM)                                                            
      M2=NINT(2.*BM)                                                            
      M3=NINT(2.*CM)                                                            
C     TEST M1+M2=M3                                                             
      IF (J1.LT.IABS(M1)) GOTO 500                                              
      IF (J2.LT.IABS(M2)) GOTO 500                                              
      IF (J3.LT.IABS(M3)) GOTO 500                                              
      IF(M1+M2-M3) 500,40,500                                                   
C     TEST TABLE SIZE                                                           
   40 I(10)=(J1+J2+J3)/2+2                                                      
      I(11)=J3+2                                                                
      N=I(10)                                                                   
      IF(N-LOGTAB) 70,70,50                                                     
   50 PRINT 60, LOGTAB,N,AJ,BJ,CJ,AM,BM,CM                                      
   60 FORMAT(' *** TABLE SIZE (',I3,') PROBLEM ***',I5,6F5.1)                   
      GO TO 500                                                                 
   70 I(1)=J1+J2-J3                                                             
      I(2)=J2+J3-J1                                                             
      I(3)=J3+J1-J2                                                             
      I(4)=J1-M1                                                                
      I(5)=J1+M1                                                                
      I(6)=J2-M2                                                                
      I(7)=J2+M2                                                                
      I(8)=J3-M3                                                                
      I(9)=J3+M3                                                                
C     CHECK I(J)=EVEN,TRIANGULAR INEQUALITY, M LESS J, FIND NO. OF TERMS        
      DO 110 J=1,9                                                              
      K=I(J)/2                                                                  
      IF(I(J)-2*K) 500,80,500                                                   
   80 IF(K) 500,90,90                                                           
   90 IF(K-N) 100,110,110                                                       
  100 N=K                                                                       
  110 I(J)=K+1                                                                  
      IF(M3) 115,400,115                                                        
  115 IL=0                                                                      
      LA=I(1)-I(5)                                                              
      LB=I(1)-I(6)                                                              
      IF(IL-LA) 120,130,130                                                     
  120 IL=LA                                                                     
  130 IF(IL-LB) 140,145,145                                                     
  140 IL=LB                                                                     
C     FORM COEFFICIENT OF SUM                                                   
  145 C=(G(I(11))-G(I(11)-1)+G(I(1))+G(I(2))+G(I(3))-G(I(10))                   
     1 +G(I(4))+G(I(5))+G(I(6))+G(I(7))+G(I(8))+G(I(9)))/2.D00                  
      J1=I(1)-IL                                                                
      J2=I(4)-IL                                                                
      J3=I(7)-IL                                                                
      M1=IL+1                                                                   
      M2=IL-LA+1                                                                
      M3=IL-LB+1                                                                
      C=C-G(J1)-G(J2)-G(J3)-G(M1)-G(M2)-G(M3)                                   
      C=DEXP(C)                                                                 
      IF(MOD(IL,2).NE.0) C=-C                                                   
      IF(N) 500,150,160                                                         
  150 CLEB=C                                                                    
      GO TO 500                                                                 
C     FORM  SUM                                                                 
  160 A=J1-1                                                                    
      B=J2-1                                                                    
      H=J3-1                                                                    
      D=M1                                                                      
      E=M2                                                                      
      F=M3                                                                      
      S=1                                                                       
      Q=N-1                                                                     
      DO 170 J=1,N                                                              
      T=(A-Q)/(D+Q)*(B-Q)/(E+Q)*(H-Q)/(F+Q)                                     
      S=1.D00-S*T                                                               
      Q=Q-1.D00                                                                 
  170 CONTINUE                                                                  
      CLEB=C*S                                                                  
      GOTO 500                                                                  
C                                                                               
C     SPECIAL FORMULA FOR M3=0 AND M1=0 OR 1/2                                  
  400 K=I(10)/2                                                                 
      IF(I(10)-2*K) 410,420,410                                                 
  410 K=1                                                                       
      GO TO 430                                                                 
  420 K=0                                                                       
  430 IF(M1) 115,440,460                                                        
  440 L=0                                                                       
      IF(K) 500,480,500                                                         
  460 IF(M1-1) 115,470,115                                                      
  470 L=1                                                                       
  480 X=L                                                                       
      M=I(3)+(I(1)+K+1)/2-L                                                     
      M1=I(10)/2+K                                                              
      M2=I(4)+I(5)                                                              
      M3=I(6)+I(7)                                                              
      J1=(I(1)+1-K)/2                                                           
      J2=(I(2)+1+K-L)/2                                                         
      J3=(I(3)+1+K-L)/2                                                         
      CLEB=DEXP((G(I(11))-G(I(11)-1)+G(I(1))+G(I(2))+G(I(3))-G(I(10)))          
     1 /2.D00+G(M1)-G(J1)-G(J2)-G(J3)                                           
     2 +X*(G(3)-(G(M2)-G(M2-1)+G(M3)-G(M3-1))/2.D00))                           
      IF(2*(M/2).NE.M) CLEB=-CLEB                                               
  500 CLEBG=CLEB                                                                
      END                                                                       
C                                                                               
C                                                                               
      DOUBLE PRECISION FUNCTION SDFUNC(AJ,AM,AK,COSBET)                         
C                                                                               
C CALLED BY SUBS. INTEG; SDFSET.                                                
C                                                                               
      PARAMETER (NIE=127)                                                       
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      COMMON /GFVLOG/ G(NIE),LOGTAB                                             
C                                                                               
      SMD=0.D0                                                                  
      JA=NINT(2.*AJ)                                                            
      MA=NINT(2.*AM)                                                            
      KA=NINT(2.*AK)                                                            
      JPM=(JA+MA)/2+1                                                           
      JMM=(JA-MA)/2+1                                                           
      JPK=(JA+KA)/2+1                                                           
      JMK=(JA-KA)/2+1                                                           
      IF (MIN(JPM,JMM).LT.1) GOTO 4                                             
      IF (MIN(JPK,JMK).LT.1) GOTO 4                                             
      MMK=(MA-KA)/2                                                             
      Y=1.D0-COSBET                                                             
      IF (Y.EQ.0.D0) GOTO 5                                                     
      MPK=(MA+KA)/2                                                             
      X=1.D0+COSBET                                                             
      IF (X.EQ.0.D0) GOTO 6                                                     
      MM=MAX(JPM,JMM)                                                           
      KK=MAX(JPK,JMK)                                                           
      M=MAX(MM,KK)                                                              
      IF (M.LE.LOGTAB) GOTO 2                                                   
      WRITE(*,1) M                                                              
    1 FORMAT(/,' TABLE SIZE IN SUBR GFV TOO SMALL. CHANGE TO ',I3)              
      STOP                                                                      
    2 KKM=JPM+JMK-2                                                             
      NMIN=MAX(0,MMK)                                                           
      NMAX=MIN(JPM,JMK)-1                                                       
      IF (NMIN.GT.NMAX) GOTO 4                                                  
      F=.5D0*(G(JPM)+G(JMM)+G(JPK)+G(JMK))                                      
      X=DSQRT(.5D0*X)                                                           
      Y=DSQRT(.5D0*Y)                                                           
      Z=1.D0                                                                    
      IF (MOD(NMIN,2).EQ.0) Z=-Z                                                
      DO 3 N=NMIN,NMAX                                                          
      Z=-Z                                                                      
      M=N+N                                                                     
      KK=KKM-M                                                                  
      MM=M-MMK                                                                  
      M=N+1                                                                     
      SMD=SMD+Z*X**KK*Y**MM*DEXP(F-G(JPM-N)-G(JMK-N)-G(M)-G(M-MMK))             
    3 CONTINUE                                                                  
    4 SDFUNC=SMD                                                                
      RETURN                                                                    
C     SPECIAL VALUES                                                            
    5 IF (MMK.NE.0) GOTO 4                                                      
      SDFUNC=1.                                                                 
      RETURN                                                                    
    6 IF (MPK.NE.0) GOTO 4                                                      
      SDFUNC=1.                                                                 
      IF (MOD(JPM,2).EQ.0) SDFUNC=-1.                                           
      END                                                                       
C                                                                               
C                                                                               
