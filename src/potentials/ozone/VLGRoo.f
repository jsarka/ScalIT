c       program VLGRoo
c       use ssplin
c       implicit real*8(a-h,o-z)
* calcul de l'energie  Electrostatique +  Dispersion pour l'etat fondamental
* gam en degree, r et R en bohr, energies en cm-1
* pour 1.8 < r < 3.3 bohr
c       data tocm/219474.63d0/
c       common/drob/q20O2,rooi
c       pi=dacos(-1.d0)
c       torad=pi/180.d0
c       rooi=0.d0
* Lecture des data
c       call inpdat
*
c       write(6,*)'# r=dist(O-O)   R=dist(O-O2), gam=angle(OG) '
c       write(6,*)'# roo(a0)  R(a0)  Gam(deg)  En(X)   En(X+SO) (cm-1)'
c       do roo = 1.8,3.3,0.2
c       do gam=0.d0,90.d0,30.d0
c       do rr=8.d0,20.d0,0.2d0
c       call vpot(roo,rr,gam,vx,vxso)
c      write(7,*)
c       write(7,*)'Apres Diagonalisation'
c       write(6,'(5g16.8)')roo,rr,gam,vx,vxso
c       enddo
c       write(6,*)
c       enddo
c       enddo
c       end
       subroutine buildmatrix(roo,rr,gam,vmf,vmfso)
       USe ssplin
       implicit real*8(a-h,o-z)
* 7 constantes: q_1^0(OH), q_2^0(OH),q_2^2(OH),Q_0(O),DelE(3P_0-^3P_1)
*  Del1O=DE(^3P_1-^3P_0), Del2O=DE(^3P_2-^3P_0),DelOH=DE(^2Pi_1/2-^2Pi_3/2)
       parameter(Nd=9,Ld=9,nfj=9,lmax=2)
       dimension VM(Nd,Nd),VMF(Nd,Nd),VMFSO(Nd,Nd)
       dimension E(9,9),F(9,9),G(9,9)
       dimension Q20(9,9),Q21(9,9),Q22(9,9),Q2m1(9,9),Q2m2(9,9)
       dimension A(9,9),B(9,9),C(9,9)
       dimension EFS(Nd,Nd),E1(Nd,Nd),E2(Nd,Nd)
       dimension xi(nr), yi(nr), bq(nr),cq(nr),dq(nr)
       common/cne/cnee1(nfj,nfj,0:lmax,0:2*lmax+1)
       common/c6j/c6jj1(nfj,nfj,0:lmax,0:2*lmax+1)
       common/drob/q20O2,rooi
* Moments en a.u. ; q20_OH=Q20(OH) 
       data Q0O/-1.90d0/
c      data q10O2,q20O2/0.0d0,-0.4546d0/
       data q10O2/0.0d0/
* Energies en cm-1
       data del1O,del2O/158.265d0,226.977d0/
       data autocm/219474.625d0/
* alpha en a0^3, energie d'ionisation in eV
       autoev=27.2113957d0
       pi=dacos(-1.d0)
       torad=pi/180.d0
       gm=gam*torad
* Initialisation des Cn et q20 en fct de roo
       if(dabs(roo-rooi).gt.1.e-06)then
       call init(roo)
* evaluation de q20O2  pour roo
       call initspline2(xi,yi,bq,cq,dq)
       q20O2=ispline(roo, xi, yi, bq, cq, dq, nr)
c       write(11,*)roo,q10O2,q20O2
c       write(11,*)
       endif
       rooi=roo
* construction de la matrice 
* element <+1| |+1>
       ipp=1
       mmb=0
c      print *,'ipp=',ipp
       do i=1,nfj
       do j=1,nfj
       xe=0.d0
       xe2=0.d0
       do il=iabs(mmb),lmax
       do mm=-il,il
       imm=mm+il
       imm2=-mm+il
       imb=mmb+il
       call matd(il,imm2,imb,gm,dd)
       xcn=cnee1(i,j,il,imm)
       xe=xe+xcn*dd/rr**(il+3)
       xc6=c6jj1(i,j,il,imm)
       xe2=xe2+xc6*dd/rr**6
*      if(dabs(xcn).gt.1.e-5)then
*      print '(4i3)',i,j,il,mm
*      print '(4f16.8)',gm,xcn,dd,xe
*      endif
       enddo
       enddo
       E1(i,j)=xe
       E2(i,j)=xe2
       enddo
       enddo
* Construction des matrices A,B,C,E,F,G,EFS
       s=dsin(gm)
       cc=dcos(gm)
       s2=dsin(gm/2.d0)
       c2=dcos(gm/2.d0)
       do i=1,9
       do j=1,9
       E(i,j)=0.d0
       F(i,j)=0.d0
       G(i,j)=0.d0
       EFS(i,j)=0.d0
       enddo
       enddo
       do i=1,Nd
       do j=1,Nd
       vm(i,j)=0.d0
       enddo
       enddo
* Matrix EFS (Structure Fine)
       EFS(1,1)=0.d0
       EFS(2,2)=0.d0
       EFS(3,3)=0.d0
       EFS(4,4)=0.d0
       EFS(5,5)=0.d0
       EFS(6,6)=del1O
       EFS(7,7)=del1O
       EFS(8,8)=del1O
       EFS(9,9)=del2O
* Matrix E
       E(1,1)=-1.d0
       E(2,2)=0.5d0
       E(3,3)=1.d0
       E(4,4)=0.5d0
       E(5,5)=-1.d0
       E(6,6)=0.5d0
       E(7,7)=-1.d0
       E(8,8)=0.5d0
       E(2,6)=-1.5d0
       E(3,9)=-dsqrt(2.d0)
       E(4,8)=1.5d0
       E(6,2)=-1.5d0
       E(8,4)=1.5d0
       E(9,3)=-dsqrt(2.d0)
* Matrix F
      F(1,2)=1.d0/dsqrt(2.d0)
      F(2,3)=1.d0/dsqrt(12.d0)
      F(3,4)=-1.d0/dsqrt(12.d0)
      F(4,5)=-1.d0/dsqrt(2.d0)
      F(6,7)=-0.5d0
      F(7,8)=0.5d0
      F(1,6)=-1.d0/dsqrt(2.d0)
      F(2,7)=0.5d0
      F(3,8)=dsqrt(3.d0)/2.d0
      F(2,9)=-dsqrt(2.d0/3.d0)
      F(6,3)=dsqrt(3.d0)/2.d0
      F(7,4)=0.5d0
      F(8,5)=-1.d0/dsqrt(2.d0)
      F(9,4)=dsqrt(2.d0/3.d0)
* Matrix G
      G(1,3)=-1.d0/dsqrt(6.d0)
      G(2,4)=-0.5d0
      G(3,5)=-1.d0/dsqrt(6.d0)
      G(1,7)=1.d0/dsqrt(2.d0)
      G(2,8)=0.5d0
      G(1,9)=-1.d0/dsqrt(3.d0)
      G(6,4)=-0.5d0
      G(7,5)=-1.d0/dsqrt(2.d0)
      G(6,8)=0.5d0
      G(9,5)=-1.d0/dsqrt(3.d0)
* Matrix Q20
       do i=1,9
       do j=1,9
       Q20(i,j)=Q0O*E(i,j)/4.d0
       enddo
       enddo 
* Matrix Q21
       do i=1,9
       do j=1,9
       Q21(i,j)=3.d0*Q0O*F(i,j)/(2.d0*dsqrt(2.d0))
       enddo
       enddo
* Matrix Q22
       do i=1,9
       do j=1,9
       Q22(i,j)=3.d0*Q0O*G(i,j)
       enddo
       enddo
* Matrix Q2m1
       do i=1,9
       do j=1,9
       Q2m1(i,j)=(-1.d0)*Q21(j,i)
       enddo
       enddo
* Matrix Q2m2
       do i=1,9
       do j=1,9
       Q2m2(i,j)=Q22(j,i)
       enddo
       enddo
* Matrix A
       do i=1,9
       do j=1,9
       a1=q10O2/(2.d0*rr**4)*(6.d0*cc*Q20(i,j)+s*(Q21(i,j)
     & -Q2m1(i,j)))
       a2=2.d0*q20O2/(48.d0*rr**5)*(72.d0*(3.d0*cc*cc-1.d0)*Q20(i,j)
     & +48.d0*s*cc*(Q21(i,j)-Q2m1(i,j))+3.d0*s*s*(Q22(i,j)
     & +Q2m2(i,j)))
*      print '(2i4,3g16.8)',i,j,a1,a2,a1+a2
       A(i,j)=a1+a2
       enddo
       enddo
      do i=1,9
      do j=1,9
      vm(i,j)=A(i,j)
      enddo
      enddo
      do i=1,Nd
      do j=1,Nd
      if(dabs(vm(i,j)-E1(i,j)).gt.1.e-8)then
c      write(6,*)'erreur de calcul de E1'
c      write(6,*)i,j,vm(i,j),E1(i,j)
      endif
      vmfso(i,j)=EFS(i,j)+(vm(i,j)+E2(i,j))*autocm
      vmf(i,j)=(vm(i,j)+E2(i,j))*autocm
      enddo
      enddo
       return
       end
      SUBROUTINE HDIAG(H,N,IVEC,U)                                      
      IMPLICIT REAL*8(A-H,O-Z)                                          
      parameter (Ld=9)
      DIMENSION H(Ld,Ld),U(Ld,Ld),DIA(Ld)                               
*     CALL ERRSET(208,300,-1)                                           
C     DIAGONALISATION EN DOUBLE PRECISION                               
      NT=0                                                              
      NM=N-1                                                            
      U(N,N)=1.0                                                        
      IF(NM)99,99,12                                                    
   12 IF(IVEC)161,10,161                                                
   10 DO 100 I=1,NM                                                     
      U(I,I)=1.0                                                        
      IPU=I+1                                                           
      DO 100 J=IPU,N                                                    
      U(I,J)=0.                                                         
      U(J,I)=0.                                                         
  100 CONTINUE                                                          
 161  DO 1022 IX=1,N                                                    
 1022  DIA(IX)=H(IX,IX)                                                 
      DO 1020 IX=1,NM                                                   
      IXP=IX+1                                                          
      DO 102 JX=IXP,N                                                   
      ABH= DABS(H(IX,JX))                                               
  104 H12=H(IX,JX)                                                      
      H(IX,JX)=0.                                                       
      H1=H(IX,IX)                                                       
      H2=H(JX,JX)                                                       
      DIFR=H1-H2                                                        
      IF(DIFR)108,107,108                                               
  107 H(IX,IX)=H1+H12                                                   
      H(JX,JX)=H1-H12                                                   
      COSN=0.707106781186547                                            
      SINUS=COSN                                                        
      GO TO 109                                                         
  108 DH12=H12+H12                                                      
      TANG= DABS(DIFR)+ DSQRT(DIFR*DIFR+DH12**2)                        
      TANG=DH12/ DSIGN(TANG,DIFR)                                       
      COSV=1.+TANG**2                                                   
      H(IX,IX)=(H1+(DH12+TANG*H2)*TANG)/COSV                            
      H(JX,JX)=(H2-(DH12-TANG*H1)*TANG)/COSV                            
      COSN=1./ DSQRT(COSV)                                              
      SINUS=TANG*COSN                                                   
  109 CONTINUE                                                          
      DO 1021 IC=1,N                                                    
      IXV=IX                                                            
      ICV=IC                                                            
      JXV=JX                                                            
      ICW=IC                                                            
      IF(IC-IX)204,103,205                                              
  204 IXV=IC                                                            
      ICV=IX                                                            
  205 HGA=H(IXV,ICV)                                                    
      IF(IC-JX)202,103,106                                              
  202 ICW=JX                                                            
      JXV=IC                                                            
  106 HGAJ=H(JXV,ICW)                                                   
      H(IXV,ICV)=HGA*COSN+HGAJ*SINUS                                    
      H(JXV,ICW)=HGAJ*COSN-HGA*SINUS                                    
  103 CONTINUE                                                          
      IF(IVEC)1021,210,1021                                             
  210 UGA=U(IC,IX)                                                      
      U(IC,IX)=SINUS*U(IC,JX)+COSN*UGA                                  
      U(IC,JX)=-SINUS*UGA+COSN*U(IC,JX)                                 
 1021 CONTINUE                                                          
  102 CONTINUE                                                          
 1020 CONTINUE                                                          
      DO 13 IX=1,N                                                      
      IF( DABS(H(IX,IX)-DIA(IX))-1. E-10)13,13,170                      
   13 CONTINUE                                                          
      GO TO191                                                          
  170 NT=NT+1                                                           
  169 IF(NT-25)161,161,19                                               
   19 CONTINUE                                                          
 191  DO 209 I=1,NM                                                     
      IP=I+1                                                            
      DO 209 J=IP,N                                                     
      IF(H(I,I) -H(J,J))212,209,209                                     
  212 TH=H(I,I)                                                         
      H(I,I)=H(J,J)                                                     
      H(J,J)=TH                                                         
      IF(IVEC)209,218,209                                               
  218 DO 208 K=1,N                                                      
      UT=U(K,J)                                                         
      U(K,J)=U(K,I)                                                     
      U(K,I)=UT                                                         
  208 CONTINUE                                                          
  209 CONTINUE                                                          
      DO 14 IX=1,N                                                      
      HBS= DABS(H(IX,IX))                                               
      NBS=HBS+6.E-10                                                    
      ABE=NBS                                                           
      IF( DABS(ABE-HBS)-6.E-10)15,15,16                                 
   15 H(IX,IX)= DSIGN(ABE,H(IX,IX))                                     
   16 DO 17 JX=1,N                                                      
      HBS= DABS(U(IX,JX))                                               
      NBS=HBS+4.D-10                                                    
      ABE=NBS                                                           
      IF( DABS(ABE-HBS)-5.E-10)18,18,17                                 
   18 U(IX,JX)= DSIGN(ABE,U(IX,JX))                                     
   17 CONTINUE                                                          
   14 CONTINUE                                                          
   99 RETURN                                                            
      END                                                               
      subroutine matd(il,imm,imb,gm,dd)
      implicit real*8(a-h,o-z)
      dcs=dcos(gm)
      dsn=dsin(gm)
      mb=imb-il
      mm=imm-il
      cc=1.d0
      if(mb.gt.mm)then
      cc=(-1)**(mb-mm)
      ii=mm
      mm=mb
      mb=ii
      endif
      if((iabs(mb).gt.iabs(mm)).or.(mb.lt.0.and.mm.lt.0))then
      ii=mm
      mm=-mb
      mb=-ii 
      endif
      if(il.eq.0)then
      dd=1.d0
      else if(il.eq.1)then
      	if(mm.eq.0.and.mb.eq.0)then
      	dd=dcos(gm)
      	else if(mm.eq.1.and.mb.eq.-1)then
      	dd=(1.d0-dcos(gm))/2.d0
      	else if(mm.eq.1.and.mb.eq.0)then
      	dd=(-1.d0)*dsin(gm)/dsqrt(2.d0)
      	else if(mm.eq.1.and.mb.eq.1)then
      	dd=(1.d0+dcos(gm))/2.d0
      	else
      	print *,'erreur de mm,mb',il,mm,mb
      	stop
      	endif
      else if(il.eq.2)then
      	if(mm.eq.0.and.mb.eq.0)then
      	dd=1.5d0*dcos(gm)*dcos(gm)-0.5d0
      	else if(mm.eq.1.and.(mb+1).eq.0)then
      	dd=(1.d0-dcos(gm))*(2.d0*dcos(gm)+1.d0)/2.d0
      	else if(mm.eq.1.and.mb.eq.0)then
      	dd=(-1.d0)*dsqrt(1.5d0)*dsin(gm)*dcos(gm)
      	else if(mm.eq.1.and.mb.eq.1)then
      	dd=(1.d0+dcos(gm))*(2.d0*dcos(gm)-1.d0)/2.d0
      	else if(mm.eq.2.and.(mb+2).eq.0)then
      	dd=(1.d0-dcos(gm))**2/4.d0
      	else if(mm.eq.2.and.(mb+1).eq.0)then
      	dd=(-1.d0)*(1.d0-dcos(gm))*dsin(gm)/2.d0
      	else if(mm.eq.2.and.mb.eq.0)then
      	dd=dsqrt(6.d0)*dsin(gm)*dsin(gm)/4.d0
      	else if(mm.eq.2.and.mb.eq.1)then
      	dd=(-1.d0)*(1.d0+dcos(gm))*dsin(gm)/2.d0
      	else if(mm.eq.2.and.mb.eq.2)then
      	dd=(1.d0+dcos(gm))**2/4.d0
      	else
      	print *,'erreur de mm,mb',il,mm,mb
      	stop
      	endif
      else
      print *,'erreur de il'
      stop
      endif
      dd=dd*cc
      return
      end
      subroutine init(roo)
       implicit real*8(a-h,o-z)
       parameter(nfj=9,lmax=2)
       common/cne/cnee1(nfj,nfj,0:lmax,0:2*lmax+1)
       common/c6j/c6jj1(nfj,nfj,0:lmax,0:2*lmax+1)
       common/dcc/ic
       common/droo/rooc
       ic=0
       rooc=roo
* initialisation des matrices Cn en couplage jj
       ipp=1
       call cnjel(ipp,cnee1) 
* initialisation des matrices C6 en couplage jj
       ipp=1
       call c6jmat(ipp,c6jj1) 
      return
      end
      subroutine vpot(roo,rr,gam,vx,vxso)
       implicit real*8(a-h,o-z)
       parameter(Nd=9,Ld=9)
       dimension VM(Ld,Ld),VL(Ld,Ld),VMSO(Ld,Ld)
       ivec=1
       call buildmatrix(roo,rr,gam,VM,VMSO)
       call hdiag(VM,Nd,ivec,Vl)
       vx=vm(9,9)
c       write(16,'(21g16.8)')roo,rr,gam,(vm(i,i),i=1,Nd)
       call hdiag(VMSO,Nd,ivec,Vl)
       vxso=vmso(9,9)
c       write(17,'(21g16.8)')roo, rr,gam,(vmso(i,i),i=1,Nd)
      return
      end
