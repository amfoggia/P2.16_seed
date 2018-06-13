c ================================
	program D1Q3
c       advection-diffusion-reaction  
c       drho/dt + U drho/dx = D d^2 rho/dx^2 + ar*rho-br*rho*rho
c       three discrete speeds 0,+1,-1  = f0,f1,f2 
c ================================
	implicit double precision(a-h,o-z)
	include 'd1q3.par'
c -----------------------------------------------------------
c  static frames: fort.51,52 ...
c  movie: movie.out
        iout  = 50
        iseed = 7953167
        call input	

        call init(iseed)   

        call hydrovar
        call equili
    	do istep = 1,nsteps
          call pbc
          call move
          call hydrovar
          call equili
          call colli
          call force
          if (mod(istep,nsteps/100).eq.1) then
            write(6,*) 'calling diagno at time ',istep
            call diagno(istep)
            call movie
          endif
c 3 profiles only
	  if (mod(istep,nsteps/2).eq.1) then
            call profil(iout)
	  endif
        enddo

        stop
	end
c ==================================
	subroutine input
c ==================================
	implicit double precision(a-h,o-z)
	include 'd1q3.par'
c---------------------------------------------------
        iunit=7
        open(unit=iunit,file='d1q3.inp') 
        open(99,file='movie.out')

        read(iunit,*) dif                  ! diffusivity 
        read(iunit,*) ar,dr,br,cr          ! reaction coefficients
        read(iunit,*) sigma                ! width of the initial density distribution
        read(iunit,*) rho0,rhoinl,rhoout   ! initial density and inlet outlet densities for open flow sim's
        read(iunit,*) u0                   ! initial velocity
        read(iunit,*) efield               ! external electric field (gravity)
        read(iunit,*) nsteps,ndiag       

        omega = 1./(0.5+3.0*dif)           ! LB relaxation parameter
        rho1  = 2*rho0
        write(6,*) 'dif and omega',dif,omega
        write(6,*) 'sigma and amu',sigma,amu
  	
	return
	end
c ===============================
	subroutine init(iseed)
c initial conditions for the populations
c ===============================
	implicit double precision(a-h,o-z)
	include 'd1q3.par'
c--------------------------------------------------
        iseed = 9353179
        pi   = 4.0*atan(1.0)
        w(0) = 2./3.
        w(1) = 1./6.
        w(2) = w(1)

        amp= 1./sqrt(2.*pi)/sigma
    	do i=1,nx
          x = float(i-nx/2)/sigma
          form   = exp(-0.5*x*x)
	  if (i.lt.nx/2) then
	   rho(i) = rho0+0.5d0	!amp*form
	  else
	   rho(i) = rho0-0.5d0
	  end if
	  u(i)   = u0
          uinit(i)=u(i)
c we use the local equilibrium with the input density and flow fields
          f(0,i)=rho(i)*w(0)
          f(1,i)=rho(i)*w(1)*(1+3.0*u(i))
          f(2,i)=rho(i)*w(2)*(1-3.0*u(i))
c local chemical rate
          rr = 2.0*ranpang(iseed)-1.
          arl(i) = ar !+ dr*rr
          brl(i) = br
	  crl(i) = cr
          write (6,*) 'ar,br,cr',arl(i),brl(i),crl(i)
        enddo
        
	return	
	end
c ==================================
	subroutine pbc
c ==================================
	implicit double precision(a-h,o-z)
	include 'd1q3.par'
c --------------------------------------
c periodic BC
        f(1,0)    = f(1,nx)
        f(2,nx+1) = f(2,1)
        return
        end
c ==================================
	subroutine move
c ==================================
! Stream	
	implicit double precision(a-h,o-z)
	include 'd1q3.par'
c---------------------------------------------
        do i = nx,1,-1      ! note that the loop runs oposite to the particle motion
	  f(1,i) = f(1,i-1)
        end do
        do i = 1,nx
	  f(2,i) = f(2,i+1)
        end do

	return	
	end
c =====================================
	subroutine hydrovar
c =====================================
	implicit double precision(a-h,o-z)
	include 'd1q3.par'
c form the macroscopic fields (hydro-variables)
c----------------------------------------
	do i = 1, nx
          rho(i)=f(0,i)+f(1,i)+f(2,i) 
          u(i)=(f(1,i)-f(2,i))/rho(i) 
	enddo

	return
	end
c ========================================
	subroutine equili
c ========================================
	implicit double precision(a-h,o-z)
	include 'd1q3.par'
c form local equilibria
c-------------------------------------------------
        cs2  = 1./3.
        do i = 1, nx
c since we deal with advectio in an imposed field ui is given as an input
         ui       = uinit(i)/cs2
         feq(0,i) = w(0)*rho(i)
         feq(1,i) = w(1)*rho(i)*(1.0d0 + ui)
         feq(2,i) = w(2)*rho(i)*(1.0d0 - ui)
        enddo

	return
	end
c =======================================================
	subroutine colli
c =======================================================
	implicit double precision(a-h,o-z)
	include 'd1q3.par'
c----------------------------------------------------------
        do i = 1, nx
	 f(0,i) = (1.0-omega)*f(0,i)+omega*feq(0,i)
	 f(1,i) = (1.0-omega)*f(1,i)+omega*feq(1,i)
	 f(2,i) = (1.0-omega)*f(2,i)+omega*feq(2,i)
	enddo

	return 
	end  
c =======================================================
	subroutine force
c =======================================================
	implicit double precision(a-h,o-z)
	include 'd1q3.par'
c----------------------------------------------------------
        frce = efield
c external force and chemical source
        do i=1,nx
         chr    = arl(i)-brl(i)*rho(i)+crl(i)*rho(i)*rho(i)
         srce   = chr*rho(i)
         f(0,i) = f(0,i)+w(0)*srce
         f(1,i) = f(1,i)+w(1)*(srce+frce)
         f(2,i) = f(2,i)+w(2)*(srce-frce)
        end do

	return 
	end
c ===================================
	subroutine diagno(istep)
c ===================================
	implicit double precision(a-h,o-z)
        include 'd1q3.par'
c diagnostic output
c----------------------------------------------------------
        t = float(istep)
        rhomass=0.
        umass  = 0.
        do i = 1,nx
         rhomass=rhomass+(rho(i)-rho0)   ! subtract the background
         umass  = umass+u(i)
        end do
        write(6,*) 'mass ',rhomass,umass
c moments
        ave  = 0.
        ptot = 0.
        do i = 1,nx
         x = float(i)
         pi = (rho(i)-rho0)/rhomass
         ptot=ptot+pi
         ave = ave + pi*x
        end do
        var=0.
        do i = 1,nx
         x = float(i)-ave
         pi = (rho(i)-rho0)/rhomass
         var = var + pi*x*x
        end do
c analytical
        varan = sigma*sigma+dif*t
        write(6,*) 'ptot, ave, var', ptot,ave,var,varan

	return
	end
c ========================
	subroutine profil(iout)
c ========================
	implicit double precision(a-h,o-z)
        include 'd1q3.par'
c density profiles
c----------------------------------------------------------
        do i = 1,nx
         write(iout,*) i,rho(i)-rho0,u(i)
        end do
        write(iout,'(bn)')

        return
        end
c ========================
	subroutine movie
c ========================
	implicit double precision(a-h,o-z)
        include 'd1q3.par'
c----------------------------------------------------------
        do i = 1,nx
         write(99,*) i,rho(i)-rho0,u(i)
        end do
        write(99,'(bn)')
        write(99,'(bn)')

        return
        end
c ============================          
         function ranpang(iseed)
c random number generator, from T. Pang, p.47         
	 implicit double precision(a-h,o-z)
c ============================          
         i2e30 = 2**30
         ia=16807
         ic=i2e30-1    ! ic=2**31-1, but 2**31 is a overflow
         ic=ic+i2e30
          
         iq=ic/ia
         ir=mod(ic,ia)

         ih=iseed/iq
         il=mod(iseed,iq)

         it=ia*il-ir*ih
         if(it.gt.0) then
           iseed=it
         else
           iseed=it+ic
         endif

         ranpang=iseed/float(ic)

         return
         end
