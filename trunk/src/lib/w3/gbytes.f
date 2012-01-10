      SUBROUTINE GBYTE (IN,IOUT,ISKIP,NBYTE)
	integer in(*), iout(*)
      CALL GBYTES (IN,IOUT,ISKIP,NBYTE,0,1)
      RETURN
      END

      SUBROUTINE GBYTES (IN,IOUT,ISKIP,NBYTE,NSKIP,N)
C          Get bytes - unpack bits:  Extract arbitrary size values from a
C          packed bit string, right justifying each value in the unpacked
C          array.
      DIMENSION IN(*), IOUT(*)
C            IN    = packed array input
C            IO    = unpacked array output
C            ISKIP = initial number of bits to skip
C            NBYTE = number of bits to take
C            NSKIP = additional number of bits to skip on each iteration
C            N     = number of iterations
C************************************** MACHINE SPECIFIC CHANGES START HERE
C          Machine dependent information required:
C            LMWD   = Number of bits in a word on this machine
C            MASKS  = Set of word masks where the first element has only the
C                     right most bit set to 1, the second has the two, ...
C            LEFTSH = Shift left bits in word M to the by N bits
C            RGHTSH = Shift right
C            OR     = Logical OR (add) on this machine.
C            AND    = Logical AND (multiply) on this machine
C          This is for Sun UNIX Fortran, DEC Alpha, and RS6000
      PARAMETER (LMWD=32)
      DIMENSION MASKS(LMWD)
      SAVE      MASKS
      DATA      MASKS /Z'1',Z'3',Z'7',Z'F', Z'1F',Z'3F',Z'7F',Z'FF',
     +Z'1FF',Z'3FF',Z'7FF',Z'FFF', Z'1FFF',Z'3FFF',Z'7FFF',Z'FFFF',
     +Z'1FFFF',       Z'3FFFF',       Z'7FFFF',       Z'FFFFF',
     +Z'1FFFFF',      Z'3FFFFF',      Z'7FFFFF',      Z'FFFFFF',
     +Z'1FFFFFF',     Z'3FFFFFF',     Z'7FFFFFF',     Z'FFFFFFF',
     +Z'1FFFFFFF',    Z'3FFFFFFF',    Z'7FFFFFFF',    Z'FFFFFFFF'/
C    +Z'1FFFFFFFF',   Z'3FFFFFFFF',   Z'7FFFFFFFF',   Z'FFFFFFFFF',
C    +Z'1FFFFFFFFF',  Z'3FFFFFFFFF',  Z'7FFFFFFFFF',  Z'FFFFFFFFFF',
C    +Z'1FFFFFFFFFF', Z'3FFFFFFFFFF', Z'7FFFFFFFFFF', Z'FFFFFFFFFFF',
C    +Z'1FFFFFFFFFFF',Z'3FFFFFFFFFFF',Z'7FFFFFFFFFFF','FFFFFFFFFFFF',
C    +Z'1FFFFFFFFFFFF',   Z'3FFFFFFFFFFFF',   Z'7FFFFFFFFFFFF',
C    +                                        Z'FFFFFFFFFFFFF',
C    +Z'1FFFFFFFFFFFFF',  Z'3FFFFFFFFFFFFF',  Z'7FFFFFFFFFFFFF',
C                                             Z'FFFFFFFFFFFFFF',
C    +Z'1FFFFFFFFFFFFFF', Z'3FFFFFFFFFFFFFF', Z'7FFFFFFFFFFFFFF',
C                                             Z'FFFFFFFFFFFFFFF',
C    +Z'1FFFFFFFFFFFFFFF',Z'3FFFFFFFFFFFFFFF',Z'7FFFFFFFFFFFFFFF',
C                                             'FFFFFFFFFFFFFFFF'/
C          IBM PC using Microsoft Fortran uses different syntax:
C     DATA MASKS/16#1,16#3,16#7,16#F,16#1F,16#3F,16#7F,16#FF,
C    + 16#1FF,16#3FF,16#7FF,16#FFF,16#1FFF,16#3FFF,16#7FFF,16#FFFF,
C    + 16#1FFFF,16#3FFFF,16#7FFFF,16#FFFFF,16#1FFFFF,16#3FFFFF,
C    + 16#7FFFFF,16#FFFFFF,16#1FFFFFF,16#3FFFFFF,16#7FFFFFF,16#FFFFFFF,
C    + 16#1FFFFFFF,16#3FFFFFFF,16#7FFFFFFF,16#FFFFFFFF/
      INTEGER RGHTSH, OR, AND
      LEFTSH(M,N) = ISHFT(M,N)
      RGHTSH(M,N) = ISHFT(M,-N)
C     OR(M,N)  = M.OR.N
C     AND(M,N) = M.AND.N
C     OR(M,N)  = IOR(M,N)
C     AND(M,N) = IAND(M,N)
C************************************** MACHINE SPECIFIC CHANGES END HERE
C          History:  written by Robert C. Gammill, jul 1972.


C          NBYTE must be less than or equal to LMWD
      ICON = LMWD-NBYTE
      IF (ICON.LT.0) RETURN
      MASK = MASKS (NBYTE)
C          INDEX  = number of words into IN before the next "byte" appears
C          II     = number of bits the "byte" is from the left side of the word
C          ISTEP  = number of bits from the start of one "byte" to the next
C          IWORDS = number of words to skip from one "byte" to the next
C          IBITS  = number of bits to skip after skipping IWORDS
C          MOVER  = number of bits to the right, a byte must be moved to be
C                   right adjusted
      INDEX = ISKIP/LMWD
      II    = MOD (ISKIP,LMWD)
      ISTEP = NBYTE+NSKIP
      IWORDS= ISTEP/LMWD
      IBITS = MOD (ISTEP,LMWD)

      DO 6 I=1,N                                                                
      MOVER = ICON-II
      IF (MOVER) 2,3,4

C          The "byte" is split across a word break.
    2 MOVEL = -MOVER
      MOVER = LMWD-MOVEL
      NP1 = LEFTSH (IN(INDEX+1),MOVEL)
      NP2 = RGHTSH (IN(INDEX+2),MOVER)
      IOUT(I) = IAND (IOR (NP1,NP2) , MASK)
      GO TO 5                                                                   

C          The "byte" is already right adjusted.
    3 IOUT(I) = IAND (IN (INDEX+1) , MASK)
      GO TO 5                                                                   

C          Right adjust the "byte".
    4 IOUT(I) = IAND (RGHTSH (IN (INDEX+1),MOVER) , MASK)

    5 II = II+IBITS
      INDEX = INDEX+IWORDS
      IF (II .LT. LMWD) GO TO 6
      II = II-LMWD
      INDEX = INDEX+1
    6 CONTINUE                                                                  

      RETURN                                                                    
      END                                                                       

