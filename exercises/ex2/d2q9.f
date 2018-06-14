css     Two-dimensional Lattice BGK 
css     sauro succi, adapted for ac274, fall 2007
c ***********************************************
c D2Q9 lattice, BGK version
c 0=rest particles, 1-4, nearest-neigh (nn), 5-8 (nnn)
c ===================
       program d2q9
c ===================
       implicit double precision(a-h,o-z)
       include'd2q9.par'
c ---------------------------------------------------------------
        open(5,file='d2q9.inp')
        open(35,file='porous')
c ------------------------------------------------------
	call input         ! input parameters
	call inithydro     ! initialize hydro fields
	call equil         ! local equils
        call initpop       ! initial Lb populations

c ------- time loop
        iconf=0
	do 10 istep = 1,nsteps
c periodic boundary conditions
css     call pbc
c mixed boundary conditions (Poiseuille flow)
         call mbc
         call move         ! streaming
         call hydro        ! hydrofields
         call equil        ! local equila
         call colli
         if (iforce)  then
          call force(istep,frce)        ! volume force
         endif
c obstacle ?
	 if (iobst.eq.1) call obst      ! thin plate
	 if (iobst.eq.2) call porous    ! porous media
c 0-dim diagnostic
	 if (mod(istep,ndiag).eq.0) then
	   call diag0D(istep)           ! volume-average diagnostics
	 endif
c 1d movie for gnuplot
         if(mod(istep,500).eq.0) then
          call movie                    ! fluid configuration
         endif
c 1d profiles
	 if (mod(istep,nout).eq.0) then
	    call profil(istep,frce)    ! selected flow profiles
	 endif
c 2d configs
	 if (mod(istep,nout).eq.0) then
	    call velvor(istep,iconf)
            call string(istep,iconf)
            iconf=iconf+1
	 endif
c-------- end of time loop
10	continue
	end
c =========================
	subroutine input
c =========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c---------------------------------------------------
	print*,' Number of steps'
        read(5,*)nsteps
	print*,' Number of steps between printing profile'
        read(5,*)nout
	print*,' Number of steps between performing diagnostics'
	read(5,*)ndiag
        print*,' Relaxation frequency omega'
        read(5,*)omega
	print*,'Applied force  (.TRUE. or .FALSE.) ?'
	read(5,*)iforce 
	print*,' Initial density and velocity for the Poiseuille force'
	read(5,*)rho0,u0,v0
	print*,' Final velocity for the Poise force'
	read(5,*) uf
	print*,' Linear obstacle (T or F?'
	read(5,*) iobst
	print*,' Obstacle heigth ?'
	read(5,*) nobst
        print*,' Obstacle id',iobst
        print*,' Length of the obstacle (multiple of 2)',nobst
        print*,' File for output: 5 chars'
        read(5,'(A)')fileout

        open(10,file=fileout//'.uy')
        open(11,file=fileout//'.vy')
        open(12,file=fileout//'.rhoy')
        open(14,file=fileout//'.uvx')
        open(16,file=fileout//'.pop')
        open(50,file=fileout//'.probe')

	print*,'*****************************************'
	print*,' Lattice BGK model, 2D with 9 velocities'
	print*,'*****************************************'
	print*,'Number of cells :',nx,'*',ny
	print*,'Nsteps :',nsteps
	print*,'Relaxation frequency :',omega
	print*,'Initial velocity for this Poiseuille force :',u0
	print*,'Initial density :',rho0
	print*,'Applied force :',iforce
	if (iobst.eq.1) then
	    print*,' Linear Obstacle with length :',nobst
	endif
	write(6,'(A)')'Output file :',fileout
c lattice weights
        w0  = 4.0d0/9.0d0
	w1  = 1.0/9.0d0
	w2  = 1.0/36.0d0
c sound-speed and related constants
	cs2  = 1.0d0 / 3.0d0
	cs22 = 2.0d0 * cs2
	cssq = 2.0d0 / 9.0d0
c viscosity and nominal Reynolds
	visc = (1.0d0 / omega - 0.5d0) * cs2
        rey  = u0*ny/visc
	print*,' Viscosity and nominal Reynolds:',visc,rey

        if(visc.lt.0) stop "OMEGA OUT of (0,2) interval!!'"

c Applied force (based on Stokes problem)

	fpois = 8.0d0 * visc * uf / dfloat(ny) / dfloat(ny)
        fpois = rho0*fpois/6.  ! # of biased populations
	print*,' Intensity of the applied force ',fpois
	
	return
	end
c ==============================
	subroutine inithydro
c ==============================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c---------------------------------------------------
        write(6,*) 'u0',u0
	do j = 0, ny+1
	  do i = 0, nx+1
           rho(i,j)  = rho0
	   u(i,j)    = u0
	   v(i,j)    = v0
	  enddo
        enddo
        if(iobst.eq.2) then
         write(6,*) 'iobst from inithydro',iobst
         call iniporous
        endif
	return	
	end
c ============================
	subroutine initpop
c ============================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c---------------------------------------------------
        iseed = 15391
        ramp  = 0.0      ! random amplitude
	do j = 0, ny+1
	 do i = 0, nx+1
          do ip=0,npop-1
           f(ip,i,j)=feq(ip,i,j)
           rr=2.*rand(iseed)-1.
           f(ip,i,j)=(1.+ramp*rr)*rho0/float(npop)
           iseed = iseed+1
          end do
         end do
        end do

        return
        end
c ========================
	subroutine move
c ========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c---------------------------------------------
	do j = ny,1,-1
	 do i = 1, nx
          f(2,i,j) = f(2,i,j-1)
          f(6,i,j) = f(6,i+1,j-1)
	 enddo
        enddo
	do j = ny,1,-1
	 do i = nx,1,-1
	  f(1,i,j) = f(1,i-1,j)
	  f(5,i,j) = f(5,i-1,j-1)
	 enddo
	enddo
	do j = 1,ny
	 do i = nx,1,-1
          f(4,i,j) = f(4,i,j+1)
	  f(8,i,j) = f(8,i-1,j+1)
	 enddo
	enddo
	do j = 1,ny
	 do i = 1,nx
	  f(3,i,j) = f(3,i+1,j)
	  f(7,i,j) = f(7,i+1,j+1)
         enddo
	enddo

	return
	end
c =========================
	subroutine hydro
c =========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c--------------------------------------------
c hydro variables
	do j = 1,ny
	 do i = 1,nx
	  rho(i,j)=f(1,i,j)+f(2,i,j)+f(3,i,j) 
     .            +f(4,i,j)+f(5,i,j)+f(6,i,j) 
     .            +f(7,i,j)+f(8,i,j)+f(0,i,j) 
          rhoi=1./rho(i,j)

	  u(i,j)=(f(1,i,j)-f(3,i,j)+f(5,i,j)- 
     .            f(6,i,j)-f(7,i,j)+f(8,i,j))*rhoi 

	  v(i,j)=(f(5,i,j)+f(2,i,j)+f(6,i,j)
     .          - f(7,i,j)-f(4,i,j)-f(8,i,j))*rhoi
	  enddo
	enddo
		
	return
	end
c ========================
	subroutine equil
c ========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c-------------------------------------------------
c equils are written explicitly to avoid multiplications by zero
	do j = 0,ny+1
	 do i = 0,nx+1
      rl = rho(i,j)
	  ul = u(i,j)/cs2
	  vl = v(i,j)/cs2
	  uv = ul*vl
	  usq = u(i,j)*u(i,j) 
	  vsq = v(i,j)*v(i,j)
	  sumsq  = (usq+vsq)/cs22
	  sumsq2 = sumsq*(1.0d0-cs2)/cs2
	  u2 = usq/cssq 
          v2 = vsq/cssq

	  feq(0,i,j) = w0*(1.0d0 - sumsq)

	  feq(1,i,j) = w1*(1.0d0 - sumsq  + u2 + ul)
	  feq(2,i,j) = w1*(1.0d0 - sumsq  + v2 + vl)
	  feq(3,i,j) = w1*(1.0d0 - sumsq  + u2 - ul)
	  feq(4,i,j) = w1*(1.0d0 - sumsq  + v2 - vl)

	  feq(5,i,j) = w2*(1.0d0 + sumsq2 + ul + vl + uv)
	  feq(6,i,j) = w2*(1.0d0 + sumsq2 - ul + vl - uv)
	  feq(7,i,j) = w2*(1.0d0 + sumsq2 - ul - vl + uv)
	  feq(8,i,j) = w2*(1.0d0 + sumsq2 + ul - vl - uv)
	 enddo
	enddo
		
	return
	end
c =========================
	subroutine colli
c =========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c----------------------------------------------------------
	do k = 0,npop-1
	 do j = 1,ny
	  do i = 1,nx
	   f(k,i,j)=f(k,i,j)*(1.0d0-omega)+omega*feq(k,i,j)
	  enddo
	 enddo
	enddo
	    
	return 
	end  
c ===================================
	subroutine force(it,frce)
c ===================================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c--------------------------------------------------------
        frce = fpois
	do j = 1,ny
	 do i = 1,nx
	  f(1,i,j) = f(1,i,j) + frce
	  f(5,i,j) = f(5,i,j) + frce
	  f(8,i,j) = f(8,i,j) + frce
	  f(3,i,j) = f(3,i,j) - frce
          f(6,i,j) = f(6,i,j) - frce
          f(7,i,j) = f(7,i,j) - frce
	 enddo
	enddo

	return
	end
c =========================
	subroutine pbc
c =========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c-----------------------------------------------------------
c EAST
	do j = 1,ny
	 f(1,0,j) = f(1,nx,j)
	 f(5,0,j) = f(5,nx,j)
	 f(8,0,j) = f(8,nx,j)
	enddo
c WEST 
	do j = 1,ny
	 f(3,nx+1,j) = f(3,1,j)
	 f(6,nx+1,j) = f(6,1,j)
	 f(7,nx+1,j) = f(7,1,j)
	enddo
c NORTH 
	do i = 1,nx
	 f(2,i,0) = f(2,i,ny)
	 f(5,i,0) = f(5,i,ny)
         f(6,i,0) = f(6,i,ny)
	enddo
c SOUTH 
	do i = 1,nx
	 f(4,i,ny+1) = f(4,i,1)
         f(7,i,ny+1) = f(7,i,1)
         f(8,i,ny+1) = f(8,i,1)
	enddo

	return
	end
c ========================
	subroutine mbc
c ========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c mixed boundary: periodic inlet/outlet, solid walls north/south
c-------------------------------------------------------------
c WEST inlet
	do j = 1,ny
         f(1,0,j) = f(1,nx,j)
         f(5,0,j) = f(5,nx,j)
         f(8,0,j) = f(8,nx,j)
        enddo
c EAST outlet
	do j = 1,ny
         f(3,nx+1,j) = f(3,1,j)
         f(6,nx+1,j) = f(6,1,j)
         f(7,nx+1,j) = f(7,1,j)
        enddo
c NORTH solid
	do i = 1,nx
	 f(4,i,ny+1) = f(2,i,ny)
	 f(8,i,ny+1) = f(6,i+1,ny)
	 f(7,i,ny+1) = f(5,i-1,ny)
	enddo
c SOUTH solid
	do i = 1,nx
         f(2,i,0) = f(4,i,1)
         f(6,i,0) = f(8,i-1,1)
         f(5,i,0) = f(7,i+1,1)
        enddo
c corners bounce-back
        f(8,0,ny+1)    = f(6,1,ny)
        f(5,0,0)       = f(7,1,1)
        f(7,nx+1,ny+1) = f(5,nx,ny)
        f(6,nx+1,0)    = f(8,nx,1)        

	return
	end
c ==========================
	subroutine obst
c ==========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c--------------------------------------------------------
	i    = nx / 4	
        jbot = ny/2-nobst/2
        jtop = ny/2+nobst/2+1
        do j = ny/2-nobst/2,ny/2+nobst/2+1
         f(1,i,j) = f(3,i+1,j)
         f(5,i,j) = f(7,i+1,j+1)
         f(8,i,j) = f(6,i+1,j-1)
         f(3,i,j) = f(1,i-1,j)
         f(7,i,j) = f(5,i-1,j)
         f(6,i,j) = f(8,i-1,j)
        enddo
c top
        f(2,i,jtop)=f(4,i,jtop+1)
        f(6,i,jtop)=f(8,i-1,jtop+1)
c bot
        f(4,i,jbot)=f(2,i,jbot-1)
        f(7,i,jbot)=f(5,i-1,jbot-1)

	return
	end
c ==========================
	subroutine iniporous
c ==========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c -----------------------------------------
        iseed = 9753571
        por = 0.1
        npore = (nx-1)*(ny-1)*por
        do j=1,ny
         do i=1,nx
           iflag(i,j) = 0
         end do
         end do

         write(6,*) 'npore from iniporous', npore
         write(35,*) npore
         do np=1,npore
          ip = nx*rand(iseed)
          iseed = iseed-1
          jp = ny*rand(iseed)
          iseed = iseed-1
          jp = ny*rand(iseed)
          iflag(ip,jp)=1
          write(35,*) np,ip,jp,iflag(ip,jp)
         end do

        return
        end
c ==========================
	subroutine porous
c ==========================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c -----------------------------------------
        do j=2,ny-1
        do i=1,nx-1
         ifl = iflag(i,j)
         if(ifl.eq.1.and.iflag(i+1,j).eq.0) then
          f(3,i+1,j)=f(1,i,j) 
         endif
         if(ifl.eq.1.and.iflag(i-1,j).eq.0) then
          f(1,i-1,j)=f(3,i,j) 
         endif
         if(ifl.eq.1.and.iflag(i,j+1).eq.0) then
          f(4,i,j+1)=f(2,i,j) 
         endif
         if(ifl.eq.1.and.iflag(i,j-1).eq.0) then
          f(2,i,j+1)=f(4,i,j) 
         endif

         if(ifl.eq.1.and.iflag(i+1,j+1).eq.0) then
          f(7,i+1,j+1)=f(5,i,j) 
         endif
         if(ifl.eq.1.and.iflag(i-1,j+1).eq.0) then
          f(8,i-1,j+1)=f(6,i,j) 
         endif
         if(ifl.eq.1.and.iflag(i-1,j-1).eq.0) then
          f(7,i-1,j-1)=f(5,i,j) 
         endif
         if(ifl.eq.1.and.iflag(i+1,j-1).eq.0) then
          f(6,i-1,j-1)=f(8,i,j) 
         endif

        end do
        end do

        return 
        end
c ===============================
	subroutine movie
c ===============================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c----------------------------------------------------------
        open(99,file='movie.out')
	do j = 1,ny
          write(99,*) j,u(nx/4,j),u(nx/2,j),u(3*nx/4,j)
	enddo
	write(99,'(bn)')
	write(99,'(bn)')

        return
        end
c ===================================
	subroutine profil(it,frce)
c ===================================
	implicit double precision(a-h,o-z)
	include'd2q9.par'
c----------------------------------------------------------
        write(6,*) 'ucenter,force',u(nx/2,ny/2),frce
	do j = 1,ny
	  write(10,*) j,u(nx/4,j),u(nx/2,j),u(3*nx/4,j)
	  write(11,*) j,v(nx/4,j),v(nx/2,j),v(3*nx/4,j)
	  write(12,*) j,rho(nx/2,j)
          write(99,*) j,u(nx/2,j)
	enddo
	write(10,'(bn)')
	write(11,'(bn)')
	write(12,'(bn)')

	do i = 1,nx
	  write(14,*) i,u(i,ny/2),v(i,ny/2)
	  write(16,*) i,f(1,i,ny/2),f(3,i,ny/2)
	enddo
	write(14,'(bn)')

        write(50,*) it,u(nx/2,ny/2)
   
	return
	end
c ================================
	subroutine diag0D(istep)
c ================================
	implicit double precision(a-h,o-z)
        include'd2q9.par'
c----------------------------------------------------------
	densit = 0.0d0
	do k = 0,npop-1
	 do j= 1,ny
	   do i = 1,nx
	     densit = densit + f(k,i,j)
	   enddo
	 enddo
	enddo
	densit = densit / dfloat(nx*ny) /dfloat(npop)
	umoy = 0.0d0
	vmoy = 0.0d0
	do j = 1,ny
	 do i = 1,nx
	    umoy = umoy + u(i,j)
	    vmoy = vmoy + v(i,j)
	 enddo
	enddo
	umoy = umoy / dfloat(nx*ny)
	vmoy = vmoy / dfloat(nx*ny)

	print*,'diagnostic 0D : istep density umoy and vmoy',
     .          istep,densit,umoy,vmoy

	return
	end
c ======================================
	subroutine velvor(istep,iconf) 
c ======================================
	implicit double precision(a-h,o-z)
        include'd2q9.par'
c -------------------------------------------
        iout=60+iconf
	do j = 2,ny-1
	  do i = 2,nx-1
            vor = u(i,j+1)-u(i,j-1)-v(i+1,j)+v(i-1,j)
            write(iout,*) i,j,u(i,j),v(i,j),vor
	   enddo
           write(iout,'(bn)')
	enddo
        write(6,*) 'configuration at time and file>>>>>> ',istep,iout
	
	return
	end
c ======================================
	subroutine string(istep,iconf) 
c ======================================
	implicit double precision(a-h,o-z)
        include'd2q9.par'
        dimension fx(0:npop),fy(0:npop),cx(0:npop)
c renumber counterclockwise to plot the contour
c  0  1  2  3  4  5  6  7  8
c  0  1  5  2  6  3  7  4  8
c -------------------------------------------
            iout=70+iconf
            i=nx/2
            j=ny-3
            fx(0)=0.
            fy(0)=0.
            fx(1)=f(1,i,j)
            fy(1)=0.0
            fx(2)=f(5,i,j)
            fy(2)=f(5,i,j)
            fx(3)=0.0
            fy(3)=f(2,i,j)
            fx(4)=-f(6,i,j)
            fy(4)=+f(6,i,j)

            fx(5)=-f(3,i,j)
            fy(5)=0.0
            fx(6)=-f(7,i,j)
            fy(6)=-f(7,i,j)
            fx(7)=0.0
            fy(7)=-f(4,i,j)
            fx(8)=+f(8,i,j)
            fy(8)=-f(8,i,j)
c close the contour
            fx(9) = fx(1)
            fy(9) = fy(1)
            do k=1,npop
             write(iout,*) i,j,fx(k),fy(k)
            end do
            write(iout,'(bn)')
            write(6,*) 'string at time and file>>>>>> ',istep,iout
	
	return
	end
