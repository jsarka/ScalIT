!c*********************************************************************
!c
        subroutine potff3(r,phi,vpot,inp,mod)
!c
!c---------------------------------------------------------------------
!c
!c   Compute the potential form bent XY_2 molecules using the force
!c   field constant model by Esa Kauppi and Lauri Halonen, see their
!c   papaer in J Phys Chem. 94, 5779(1990).
!c   The potential energy function is expanded in terms of the
!c   curviliner internal displacement coordinate theta=phi-phi_e for
!c   the bond and in terms of the Morse variables y_i=1-exp(-a*r_i)
!c   (i=1 or 2) for the stretches. r_i=R_i-R_ei. R_i and phi are the
!c   instantaneous bond length and valence angle, and R_ie and phi_e
!c   are the corresponding equilibrium values.
!c
!c   The potential energy parameters fo SO_2 are tacken from the paper
!c   by the same author above in J Chem Phys 96, 2933(1992).
!c
!c   Since the force field constants are in the units of aJ=10^(-18),
!c   we convert the final result to eV by dividing it by 0.1602.
!c
!c   Input:  R(i) (i=1 and 2, in angstrom), phi (in rad); 
!c   Output: vpot (in eV)
!c*********************************************************************
!c
        implicit real*8 (a-h,o-z)
        dimension r(2),y(2),y2(2)
        parameter (aJtoeV=1.60217733d-1)
        common /cost/pi,pi2
        character inp*3, mod*2

!c*********************************************************************
!c   parameters of potential energy for SO_2 with model AB
!c   taken from the paper by E. Kauppi and L. Halonen in
!c   J. Chem. Phys. 96(4), 2933(1992)
!c*********************************************************************

!c       data re,phie/1.4308d0,119.33d0/              ! phie in degree
!c       data re,phie/1.4308d0,2.08270139640484d0/    ! phie in rad
!c       data De,a/0.97665d0,2.2706d0/
!c       data T1,T2/0.d0,0.06206d0/
!c       data frr,frrr,frrrr/10.070d0,-68.59d0,402.99d0/
!c       data fss,fsss,fssss/1.69355d0,-2.4965d0,6.690d0/
!c       data frrp,frrrp/-0.2017d0,1.668d0/
!c       data frs,frss,frrs/-0.1669d0,-0.9742d0,2.872d0/
!c       data frrps,frrss/0.d0,-12.502d0/
!c       mod='AB'
!c       write(10,*)'AB potential model'

!c*********************************************************************
!c   parameters of potential energy for SO_2 with model FV
!c   taken from the paper by E. Kauppi and L. Halonen in
!c   J. Chem. Phys. 96(4), 2933(1992)
!c*********************************************************************

!c       data re,phie/1.4308d0,119.33d0/
        data re,phie/1.4308d0,2.08270139640484d0/    ! phie in rad
        data De,a/1.0120d0,2.2545d0/
        data T1,T2/0.d0,0.0285d0/
        data frr,frrr,frrrr/10.288d0,-69.58d0,383.72d0/
        data fss,fsss,fssss/1.65551d0,-2.6111d0,6.34d0/
        data frrp,frrrp/-0.0202d0,0.d0/
        data frs,frss,frrs/0.2511d0,-1.089d0,0.d0/
        data frrps,frrss/0.d0,-8.74d0/
!c       write(10,*)'FV potential model'
        mod='FV'

        a2=a*a
        a3=a2*a
        a4=a2*a2

        if (inp .eq. 'DTT') then
           goto 10
        elseif (inp .eq. 'frr') then
           De=frr/a2/2.d0
           T1=frrr/a3/6.d0+De
           T2=(frrrr/a4-14.d0*De+3.d0*T1)/24.d0
        else
           write(*,100)inp
100        format(' The input parameter inp should be "DTT" or "frr"',
     &      /' for inputing De,T1,T2 or frr,frrr,frrrr, respectively',
     &      /' but now it is ',a3)
           stop
        endif

 10     continue

        do i=1,2
           rmre=r(i)-re
           y(i)=1.d0-dexp(-a*rmre)
           y2(i)=y(i)*y(i)
!c          write(*,*)'r,y=',r(i),rmre,y(i)
        enddo

        sita=phi-phie
        sita2=sita*sita
!c       write(*,*)' sita=',sita,phie

        v0 = De*(y2(1)+y2(2)) + frrp/a2*y(1)*y(2)
     &     + T1*(y(1)**3+y(2)**3) + T2*(y(1)**4+y(2)**4)
     &     + 0.5d0*(frrrp/a3+frrp/a2)*(y2(1)*y(2)+y(1)*y2(2))
        v1 = fss*sita2/2.d0 + fsss*sita**3/6.d0 + fssss*sita**4/24.d0
        v2 = frs/a*(y(1)+y(2))*sita
        v3 = 0.5d0/a*frss*(y(1)+y(2))*sita2
        v4 = 0.25d0*(frrss/a2+frss/a)*(y2(1)+y2(2))*sita2
        v5 = 0.5d0*(frrs/a2+frs/a)*(y2(1)+y2(2))*sita
     &     + frrps/a2*y(1)*y(2)*sita

        vpot=v0+v1+v2+v3+v4+v5
!c       write(*,*)'v=',v0,v1,v2,v3,v4,v5,vpot
        vpot=vpot/aJtoeV

        return
        end

!c***********************************************************************
