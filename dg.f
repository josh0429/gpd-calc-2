      DOUBLE PRECISION FUNCTION DGAUS1(NPS,A,B,PE)                      
CCCCCC                                                                  
CCCCCC                  INTEGRALE UNIDIMENSIONALE                       
CCCCCC                                                                  
CCCCCC              NPS - NUMERO DI PUNTI DI GAUSS                      
CCCCCC              PE - VETTORE CONTENENTE LA FUNZIONE INTEGRANDA      
CCCCCC                                                                  
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION PE(1),W(48),WE(48)                                      
C     WRITE(8,662)                                                      
      CALL GWT(W,NPS)                                                   
C     WRITE(8,663)                                                      
      AB=0.5D0*(B-A)                                                    
      DGAUS1 = 0.D0                                                     
      DO 30 I=1,NPS                                                     
      WE(I) = W(I)*AB                                                   
   30 DGAUS1 = DGAUS1+WE(I)*PE(I)                                       
      RETURN                                                            
C661   FORMAT(2(D10.5,2X))                                              
C662   FORMAT('SUB.DGAUS1')                                             
C663   FORMAT('ESCE DA GWT')                                            
      END                                                               
CCCCCC                                                                  
      SUBROUTINE DINTER(E0,A,B,NPS)                                     
      IMPLICIT REAL*8 (A-H,O-Z)                                         
CCCCCC                                                                  
CCCCCC               IMPOSTA I PUNTI DI INTEGRAZIONE                    
CCCCCC                                                                  
      DIMENSION E0(1),POINT(48)                                         
C     CALL DCHK(A,B,NPS)                                                
C     WRITE(6,610)                                                      
C     WRITE(8,610)                                                      
      CALL GPT(POINT,NPS)                                               
      DF=0.5D0*(B-A)                                                    
      DS=0.5D0*(B+A)                                                    
      DO 20 L=1,NPS                                                     
      E0(L)=DF*POINT(L)+DS                                              
C     WRITE(6,600)E0(L)                                                 
C     WRITE(8,600)E0(L)                                                 
20    CONTINUE                                                          
C     WRITE(6,1000)                                                     
C1000 FORMAT(//,10X,'WARNING : VERIFY INTEGRATION POINTS')              
C600   FORMAT(D10.5)                                                    
C610   FORMAT('ENTRA SUB. DINTER ')                                     
      RETURN                                                            
      END                                                               
      SUBROUTINE GPT(P,N)                                               
CCCCCCC                                                                 
CCCCCCC        IMPOSTA I PUNTI DI GAUSS-LEGENDRE NEL VETTORE P          
CCCCCCC                    (N=8,12,16,24,36,48)                         
CCCCCCC                                                                 
      DOUBLE PRECISION P(1)                                             
      DATA JP/-1/                                                       
      L=N/4                                                             
      GO TO (1,2,3,4,1,6,1,1,9,1,1,12),L                                
    1 CONTINUE                                                          
      RETURN                                                            
    2 CONTINUE                                                          
      P(1)=0.9602898564975362D0                                         
      P(2)=0.7966664774136267D0                                         
      P(3)=0.5255324099163290D0                                         
      P(4)=0.1834346424956498D0                                         
      GO TO 40                                                          
    3 CONTINUE                                                          
      P(1)=0.9815606342467193D0                                         
      P(2)=0.9041172563704749D0                                         
      P(3)=0.7699026741943047D0                                         
      P(4)=0.5873179542866175D0                                         
      P(5)=0.3678314989981802D0                                         
      P(6)=0.1252334085114689D0                                         
      GO TO 40                                                          
    4 CONTINUE                                                          
      P(1)=0.9894009349916499D0                                         
      P(2)=0.9445750230732326D0                                         
      P(3)=0.8656312023878317D0                                         
      P(4)=0.7554044083550030D0                                         
      P(5)=0.6178762444026438D0                                         
      P(6)=0.4580167776572274D0                                         
      P(7)=0.2816035507792589D0                                         
      P(8)=0.9501250983763744D-01                                       
      GO TO 40                                                          
    6 CONTINUE                                                          
      P(1)=0.9951872199970214D0                                         
      P(2)=0.9747285559713095D0                                         
      P(3)=0.9382745520027328D0                                         
      P(4)=0.8864155270044010D0                                         
      P(5)=0.8200019859739029D0                                         
      P(6)=0.7401241915785544D0                                         
      P(7)=0.6480936519369756D0                                         
      P(8)=0.5454214713888395D0                                         
      P(9)=0.4337935076260451D0                                         
      P(10)=0.3150426796961634D0                                        
      P(11)=0.1911188674736163D0                                        
      P(12)=0.6405689286260563D-01                                      
      GO TO 40                                                          
    9 CONTINUE                                                          
      P(1)=0.9978304624840858D0                                         
      P(2)=0.9885864789022122D0                                         
      P(3)=0.9720276910496980D0                                         
      P(4)=0.9482729843995076D0                                         
      P(5)=0.9174977745156591D0                                         
      P(6)=0.8799298008903971D0                                         
      P(7)=0.8358471669924753D0                                         
      P(8)=0.7855762301322065D0                                         
      P(9)=0.7294891715935566D0                                         
      P(10)=0.6680012365855211D0                                        
      P(11)=0.6015676581359805D0                                        
      P(12)=0.5306802859262452D0                                        
      P(13)=0.4558639444334203D0                                        
      P(14)=0.3776725471196892D0                                        
      P(15)=0.2966849953440283D0                                        
      P(16)=0.2135008923168656D0                                        
      P(17)=0.1287361038093848D0                                        
      P(18)=0.4301819847370861D-01                                      
      GO TO 40                                                          
   12 CONTINUE                                                          
      P(1)=0.9987710072524261D+00                                       
      P(2)=0.9935301722663508D+00                                       
      P(3)=0.9841245837228269D+00                                       
      P(4)=0.9705915925462473D+00                                       
      P(5)=0.9529877031604309D+00                                       
      P(6)=0.9313866907065543D+00                                       
      P(7)=0.9058791367155697D+00                                       
      P(8)=0.8765720202742479D+00                                       
      P(9)=0.8435882616243935D+00                                       
      P(10)=0.8070662040294426D+00                                      
      P(11)=0.7671590325157403D+00                                      
      P(12)=0.7240341309238147D+00                                      
      P(13)=0.6778723796326639D+00                                      
      P(14)=0.6288673967765136D+00                                      
      P(15)=0.5772247260839727D+00                                      
      P(16)=0.5231609747222330D+00                                      
      P(17)=0.4669029047509584D+00                                      
      P(18)=0.4086864819907167D+00                                      
      P(19)=0.3487558862921607D+00                                      
      P(20)=0.2873624873554556D+00                                      
      P(21)=0.2247637903946891D+00                                      
      P(22)=0.1612223560688917D+00                                      
      P(23)=0.9700469920946270D-01                                      
      P(24)=0.3238017096286936D-01                                      
   40 M=N/2+1                                                           
      J=N+JP*M+2                                                        
      DO 50 I=M,N                                                       
      J=J+JP                                                            
      P(I)=P(J)                                                         
   50 P(J)=DFLOAT(JP)*P(J)                                              
      RETURN                                                            
      END                                                               
      SUBROUTINE DCHK(A,B,NPASSI)                                       
CCCCCC                                                                  
CCCCCC                  CHECK SUL VALORE DI NPASSI                      
CCCCCC                  (NPASSI=8,12,16,24,36,48)                       
CCCCCC                                                                  
      N=NPASSI/4                                                        
      GO TO (1,2,2,2,1,2,1,1,2,1,1,2), N                                
    1 WRITE(6,1000)                                                     
 1000 FORMAT(//,30X,'NUMERO DI PASSI NON PREVISTO',//)                  
      CALL EXIT                                                         
    2 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE GWT(W,N)                                               
CCCCCCC                                                                 
CCCCCCC            IMPOSTA I PESI DI GAUSS-LEGENDRE NEL VETTORE W       
CCCCCCC                                                                 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION W(1)                                                    
      DATA JP/-1/                                                       
      L=N/4                                                             
      GO TO (1,2,3,4,1,6,1,1,9,1,1,12),L                                
    1 CONTINUE                                                          
      RETURN                                                            
    2 CONTINUE                                                          
      W(1)=0.1012285362903763D0                                         
      W(2)=0.2223810344533745D0                                         
      W(3)=0.3137066458778873D0                                         
      W(4)=0.3626837833783620D0                                         
      GO TO 40                                                          
    3 CONTINUE                                                          
      W(1)=0.4717533638651183D-01                                       
      W(2)=0.1069393259953184D0                                         
      W(3)=0.1600783285433462D0                                         
      W(4)=0.2031674267230659D0                                         
      W(5)=0.2334925365383548D0                                         
      W(6)=0.2491470458134028D0                                         
      GO TO 40                                                          
    4 CONTINUE                                                          
      W(1)=0.2715245941175410D-01                                       
      W(2)=0.6225352393864789D-01                                       
      W(3)=0.9515851168249279D-01                                       
      W(4)=0.1246289712555339D0                                         
      W(5)=0.1495959888165767D0                                         
      W(6)=0.1691565193950025D0                                         
      W(7)=0.1826034150449236D0                                         
      W(8)=0.1894506104550685D0                                         
      GO TO 40                                                          
    6 CONTINUE                                                          
      W(1)=0.1234122979998720D-01                                       
      W(2)=0.2853138862893367D-01                                       
      W(3)=0.4427743881741981D-01                                       
      W(4)=0.5929858491543678D-01                                       
      W(5)=0.7334648141108031D-01                                       
      W(6)=0.8619016153195328D-01                                       
      W(7)=0.9761865210411389D-01                                       
      W(8)=0.1074442701159656D0                                         
      W(9)=0.1155056680537256D0                                         
      W(10)=0.1216704729278034D0                                        
      W(11)=0.1258374563468283D0                                        
      W(12)=0.1279381953467522D0                                        
      GO TO 40                                                          
    9 CONTINUE                                                          
      W(1)=0.5565719664245045D-02                                       
      W(2)=0.1291594728406557D-01                                       
      W(3)=0.2018151529773547D-01                                       
      W(4)=0.2729862149856878D-01                                       
      W(5)=0.3421381077030723D-01                                       
      W(6)=0.4087575092364490D-01                                       
      W(7)=0.4723508349026598D-01                                       
      W(8)=0.5324471397775992D-01                                       
      W(9)=0.5886014424532482D-01                                       
      W(10)=0.6403979735501549D-01                                      
      W(11)=0.6874532383573644D-01                                      
      W(12)=0.7294188500565306D-01                                      
      W(13)=0.7659841064587068D-01                                      
      W(14)=0.7968782891207160D-01                                      
      W(15)=0.8218726670433971D-01                                      
      W(16)=0.8407821897966194D-01                                      
      W(17)=0.8534668573933863D-01                                      
      W(18)=0.8598327567039475D-01                                      
      GO TO 40                                                          
   12 CONTINUE                                                          
      W(1)=0.3153346052305839D-02                                       
      W(2)=0.7327553901276262D-02                                       
      W(3)=0.1147723457923454D-01                                       
      W(4)=0.1557931572294385D-01                                       
      W(5)=0.1961616045735553D-01                                       
      W(6)=0.2357076083932438D-01                                       
      W(7)=0.2742650970835695D-01                                       
      W(8)=0.3116722783279809D-01                                       
      W(9)=0.3477722256477044D-01                                       
      W(10)=0.3824135106583071D-01                                      
      W(11)=0.4154508294346475D-01                                      
      W(12)=0.4467456085669428D-01                                      
      W(13)=0.4761665849249047D-01                                      
      W(14)=0.5035903555385447D-01                                      
      W(15)=0.5289018948519367D-01                                      
      W(16)=0.5519950369998416D-01                                      
      W(17)=0.5727729210040322D-01                                      
      W(18)=0.5911483969839564D-01                                      
      W(19)=0.6070443916589388D-01                                      
      W(20)=0.6203942315989266D-01                                      
      W(21)=0.6311419228625403D-01                                      
      W(22)=0.6392423858464819D-01                                      
      W(23)=0.6446616443595008D-01                                      
      W(24)=0.6473769681268392D-01                                      
   40 M=N/2+1                                                           
      J=N+JP*M+2                                                        
      DO 50 I=M,N                                                       
      J=J+JP                                                            
   50 W(I)=W(J)                                                         
      RETURN                                                            
      END                                                               