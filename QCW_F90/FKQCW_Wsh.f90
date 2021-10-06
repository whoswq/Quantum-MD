program main
   use mpi
   use mt19937
   implicit none
   integer, parameter :: P = 192 ! number of beads
   real*8, parameter :: hbar = 1d0, pi = acos(-1d0)
   real*8, parameter :: temp = 14d0, kb = 3.16683d-6, beta = 1d0/temp/kb
   integer , parameter :: natom = 125
   integer , parameter :: n = natom*3
   real*8, parameter :: m(n) = 1837.0d0 * 2.0d0

   !  -------------------------------------------
   ! SG-H2 potential parameter (SI vlaues)
   ! use atomic unit
   real*8,parameter::ah = 1.713d0
   real*8,parameter::bh = 1.5671d0
   real*8,parameter::rh = 0.00993d0
   real*8,parameter::rc = 8.32d0
   real*8,parameter::c6 = 12.14d0
   real*8,parameter::c8 = 215.2d0
   real*8,parameter::c9 = 143.1d0
   real*8,parameter::c10 = 4813.9d0

   real*8,parameter:: dens0 = 1.0d0/25.60d0    ! (mol/cm3)
   real*8,parameter:: clength = 5.291772d-11 ! Length (m/au)
   real*8,parameter:: avog = 6.02214d23   ! Avogadro constant 
   real*8,parameter:: dens = dens0 * 1.0d6 * clength**3 * avog              

   real*8,parameter:: box = (dble(natom)/dens) ** (1.0d0/3.0d0)  ! a.u.
   real*8,parameter:: inv_box = 1.0d0 / box

   real*8,parameter:: rcut = (box-1.0d-15*box)*0.50d0 
   real*8,parameter:: rcutsq = rcut ** 2
   real*8, parameter :: angperau = 0.52917721092d0
   !  -----------------------------------------------

   integer :: count
   real*8, parameter :: betap = beta / P
   real*8, parameter :: omegaj(0:P-1) = [0d0,(1d0/betap/hbar*2*sin(count*pi/P),count=1,P-1)] !freqencies of normal modes
   real*8, parameter :: gamad = 1d0 !adiabatic parameter
   real*8, parameter :: lgam(0:P-1) = [0.017d0,(omegaj(count),count=1,P-1)] !gamma for Langevin thermostat
   real*8, parameter :: lambda = 0.5d0 !???
   real*8, parameter :: lgamad(0:P-1) = [0d0,(omegaj(count)/sqrt(gamad),count=1,P-1)]*lambda*2 !ad gamma for Langevin thermostat
   real*8, parameter :: dt = 40d0 !timestep
   real*8, parameter :: dtr = dt*sqrt(gamad)
   real*8 :: mqj(n,0:P-1) = reshape([(m,count=0,P-1)],[n,P]) !mass of normal modes
   real*8 :: mqjad(n,0:P-1) = reshape([m,(m*gamad,count=1,P-1)],[n,P]) !ad mass of normal modes
   real*8 :: mxj(n,P) = reshape([(m,count=0,P-1)],[n,P]) !???
   real*8 :: c1(0:P-1), c2(0:P-1) !constant in Langevin thermostat
   real*8 :: c1ad(0:P-1), c2ad(0:P-1) !adiabatic version of those
   integer, parameter :: eqstep = 1d4 ! steps for equilibrium
   integer, parameter :: nout = dt / dtr !???
   integer, parameter :: ncons = 1000

   integer, parameter :: navg = 0
   integer, parameter :: nprint = 100
   integer, parameter :: ninit = 19
   integer, parameter :: samples = 1!50
   logical, parameter :: fcent = .false.
   logical, parameter :: ftrj = .true.

   integer, parameter :: nfft = 18
   integer, parameter :: corstep = 2**(nfft-1) - 1
   integer, parameter :: trj_step = 2**nfft
   integer, parameter :: nstep = trj_step + 2*navg
   integer*8, parameter :: sampstep = nstep * nout

   real*8 :: cmat(P,0:P-1)
   integer :: idof = 0, jdof = 0, ierr
   integer :: j = 0, k = 0, l = 0, sc1
   integer :: nprocs, idproc, icons
   integer*8 :: i = 0
   real*8 :: qj(n,0:P-1) = 0d0 ! normal-mode coordinates
   real*8 :: pqj(n,0:P-1) = 0d0 ! normal-mode momenta
   real*8 :: qj0(n,0:P-1) = 0d0
   real*8 :: pqj0(n,0:P-1) = 0d0
   real*8 :: fqj0(n,0:P-1) = 0d0
   real*8 :: qj1(n,0:P-1) = 0d0
   real*8 :: pqj1(n,0:P-1) = 0d0
   real*8 :: pxj(n,P) = 0d0
   real*8 :: fqj(n,0:P-1) = 0d0 ! normal-mode forces
   real*8 :: xj(n,P) = 0d0 ! Cartesian coordinates
   real*8 :: fxj(n,P) = 0d0 ! Cartesian forces
   real*8 :: xc(n) ! centroid coordinates
   real*8 :: pc(n) ! centroid momenta
   real*8 :: fxc(n) ! forces on the centroid
   real*8 :: vec1(n), vec2(n), sumfxc(n), sumhes(n,n), avghes(n,n)
   real*8 :: xeq(n)
   real*8 :: eigvec(n,n), eigval(n), eigvec_old(n,n)
   real*8 :: work(3*n-1)
   real*8 :: ct1, eta
   character*30 :: cc

   real*8 :: fxc_rec(n,-navg:navg), hes_rec(n,n,-navg:navg)
   real*8 :: omec(n), qcor(n), bhw(n), epsil = 1d-2
   real*8 :: xtil(n), ptil(n)
   real*8 :: xtmp(n), ptmp(n)
   real*8 :: xphy(n), pphy(n), vtherm(n)
   real*8 :: omec_map(n)
   real*8 :: tmp1, vc(3), lc(3), summ

   real*8 :: u0(n,P) = 0d0, Ep(P) = 0d0, mtherm(n,n), mtherm_diag(3,3,natom)
   real*8 :: tmpmat(n,n)

   logical :: zmaster
   integer :: ii
   real*8 :: estfxc(n), esthes(n,n)
   real*8, allocatable :: pqj_mpi(:,:,:), qj_mpi(:,:,:), fqj_mpi(:,:,:)
   real*8, allocatable :: hes_mpi(:,:,:), fxc_mpi(:,:)

   CALL MPI_INIT(ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,idproc,ierr)
   CALL SYSTEM_CLOCK( sc1 )
   CALL CPU_TIME( ct1 )
   CALL init_genrand_real(mod(sc1,100000)+idproc*139)
   eta = genrand_real()
   WRITE(*,*) 'The 1st random number of', idproc, eta
   eta = genrand_real()
   WRITE(*,*) 'The 2nd random number of', idproc, eta
   zmaster = (idproc == 0)
   allocate(pqj_mpi(n,0:P-1,nprocs))
   allocate(qj_mpi(n,0:P-1,nprocs))
   allocate(fqj_mpi(n,0:P-1,nprocs))
   allocate(fxc_mpi(n,nprocs))
   allocate(hes_mpi(n,n,nprocs))

   call init_cmat(cmat) !generate a transformming matrix of normal mode

   do k=1, samples
      if(zmaster) then
         write(*,'(i,a,i)') idproc, ' sample=', k

         !-----------------initial condition-------------------
         write(cc,*) ninit + k
         open(13,file='Init_'//trim(adjustl(cc)))
         read(13,*) qj
         close(13)

         call dgemm('N','T',n,P,P,1d0,qj,n,cmat,P,0d0,xj,n)
         do j = 1, P
            call getforce(xj(:,j), fxj(:,j), Ep(j) )
         end do
         fqj = matmul(fxj,cmat)

         do j = 0, P - 1
            do idof = 1, n
               call gasdev(eta)
               pqj(idof,j) = sqrt(mqj(idof,j)/betap) * eta
            end do
         end do
         !----------------------------------------------------

         do i = 1, eqstep
            call NMPIMD(pqj,qj,xj,fqj,fxj)
         end do

         open(100,file='Final_'//trim(adjustl(cc)))
         write(100,*) qj
         close(100)

         pqj(:,1:P-1) = pqj(:,1:P-1) * sqrt(gamad)

         write(cc,*) ninit + k
         if(fcent) then
            open(11,file='cent_'//trim(adjustl(cc)))
         end if
         if(ftrj) then
            open(21,file='trj_qcw_'//trim(adjustl(cc)),form='unformatted')
            open(22,file='trj_map_'//trim(adjustl(cc)),form='unformatted')
         end if
         !open(123,file='freq_'//trim(adjustl(cc)))

         sumfxc = 0d0
         sumhes = 0d0
      end if

      do i = 1, sampstep
         if(mod(i,nprocs) .eq. 1) then
            if(zmaster) then
               do ii = 1, nprocs
                  call PACMD(pqj,qj,xj,fqj,fxj)
                  pqj_mpi(:,:,ii) = pqj
                  qj_mpi(:,:,ii) = qj
                  fqj_mpi(:,:,ii) = fqj
               end do
            end if
            call MPI_SCATTER(pqj_mpi(:,:,:),n*P,MPI_DOUBLE_PRECISION,pqj(:,:),n*P,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_SCATTER(qj_mpi(:,:,:),n*P,MPI_DOUBLE_PRECISION,qj(:,:),n*P,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_SCATTER(fqj_mpi(:,:,:),n*P,MPI_DOUBLE_PRECISION,fqj(:,:),n*P,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

            esthes = 0d0
            estfxc = 0d0
            do icons = 1, ncons
               call CONS_TRPMDCOR(dt)
               do j = 1, P - 1
                  vec1 = -fqj(:,j) / sqrt(m)
                  vec2 = qj(:,j) * sqrt(m)
                  call dsyr('L',n,betap/(P-1),vec1,1,esthes,n)
                  call dsyr2('L',n,-omegaj(j)**2/2*betap/(P-1),vec1,1,vec2,1,esthes,n)
               end do
               fxc = fqj(:,0) / sqrt(dble(P))
               estfxc = estfxc + fxc
            end do

            call MPI_GATHER(estfxc(:),n,MPI_DOUBLE_PRECISION,fxc_mpi(:,:),n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(esthes(:,:),n*n,MPI_DOUBLE_PRECISION,hes_mpi(:,:,:),n*n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         end if

         if(zmaster) then
            pqj = pqj_mpi(:,:,mod(i-1,nprocs)+1)
            qj = qj_mpi(:,:,mod(i-1,nprocs)+1)
            fqj = fqj_mpi(:,:,mod(i-1,nprocs)+1)
            sumhes = sumhes + hes_mpi(:,:,mod(i-1,nprocs)+1)
            sumfxc = sumfxc + fxc_mpi(:,mod(i-1,nprocs)+1)

            if(mod(i,nout) .eq. 0) then
               hes_rec(:,:,-navg:navg-1) = hes_rec(:,:,-navg+1:navg)
               hes_rec(:,:,navg) = sumhes / nout / ncons
               sumhes = 0d0
               fxc_rec(:,-navg:navg-1) = fxc_rec(:,-navg+1:navg)
               fxc_rec(:,navg) = sumfxc / nout / ncons
               sumfxc = 0d0

               write(cc,*) n*3 + n**2*1
               if(fcent) write(11,'('//trim(adjustl(cc))//'F)') xc, pc, fxc, hes_rec(:,:,navg)
            end if

            if(i .eq. nout*(2*navg+1)) then
               !xphy = xj(:,1)
               xc = qj(:,0)/sqrt(dble(P))
               pc = pqj(:,0)/sqrt(dble(P))

               avghes = sum(hes_rec,3)/(navg*2+1)
               eigvec = avghes
               CALL dsyev('V','L',n,eigvec,n,eigval,work,3*n-1,ierr)
               if(ierr/=0) then
                  write(*,*) 'Diagonalization fails!'
                  write(*,*) 'Error information:', ierr
                  stop
               end if
               eigvec_old = eigvec
               do idof = 1, n
                  if(eigval(idof) > 0d0) then
                     omec(idof) = sqrt(eigval(idof))
                     bhw(idof) = beta*hbar*omec(idof)
                     qcor(idof) = beta*hbar*omec(idof)/2/tanh(beta*hbar*omec(idof)/2)
                  else
                     omec(idof) = 0d0
                     bhw(idof) = 0d0
                     qcor(idof) = 1d0
                  end if
               end do

               !call dgemv('T',n,n,1d0,eigvec,n,(xphy-xc)*sqrt(m),1,0d0,xtmp,1)
               !do idof = 1, n
               !   if(bhw(idof) > epsil) then
               !      xtil(idof) = xtmp(idof) / sqrt(2*(qcor(idof)-1d0)/beta/omec(idof)**2)
               !   else
               !      xtil(idof) = xtmp(idof) / &
               !         sqrt(2*hbar**2*beta*(1d0/12d0-bhw(idof)**2/720d0+bhw(idof)**4/30240d0-bhw(idof)**6/1209600d0))
               !   end if
               !end do

               do idof = 1, n
                  call gasdev(eta)
                  ptil(idof) = eta/sqrt(2d0)
                  call gasdev(eta)
                  xtil(idof) = eta/sqrt(2d0)
                  !ptmp(idof) = sqrt(2*(qcor(idof)-1d0)/beta)*ptil(idof)
               END do
               !call dgemv('N',n,n,1d0,eigvec,n,ptmp,1,0d0,pphy,1)
               !pphy = pc + pphy * sqrt(m)
               do idof = 1, natom*3
                  if(bhw(idof) > epsil) then
                     xtmp(idof) = xtil(idof) * sqrt(2*(qcor(idof)-1d0)/beta/omec(idof)**2)
                  else
                     xtmp(idof) = xtil(idof) * &
                        sqrt(2*hbar**2*beta*(1d0/12d0-bhw(idof)**2/720d0+bhw(idof)**4/30240d0-bhw(idof)**6/1209600d0))!a polynominial approximation
                  end if
                  ptmp(idof) = sqrt(2*(qcor(idof)-1d0)/beta)*ptil(idof)
               end do
               call dgemv('N',n,n,1d0,eigvec,n,xtmp,1,0d0,xphy,1)
               xphy = xc + xphy / sqrt(m)
               call dgemv('N',n,n,1d0,eigvec,n,ptmp,1,0d0,pphy,1)
               pphy = pc + pphy * sqrt(m)

               do idof = 1, n
                  tmpmat(:,idof) = eigvec(:,idof) / qcor(idof)
               end do
               call dgemm('N','T',n,n,n,1d0,tmpmat,n,eigvec,n,0d0,mtherm,n)
               do idof = 1, n
                  do jdof = 1, n
                     mtherm(idof,jdof) = mtherm(idof,jdof) / sqrt(m(idof)*m(jdof))
                  end do
               end do
               do idof = 1, natom
                  mtherm_diag(:,:,idof) = mtherm(3*idof-2:3*idof,3*idof-2:3*idof)
               end do
               call dsymv('U',n,1d0,mtherm,n,pphy,1,0d0,vtherm,1)

               if(ftrj) then
                  write(cc,*) n*8
                  write(21) xphy, pphy, xc, pc,&
                     vtherm, mtherm_diag
               end if
            end if

            if(mod(i,nout).eq.0 .and. i.gt.nout*(2*navg+1)) then

               do idof = 1, n
                  tmp1 = xtil(idof)*cos(omec(idof)*dt) + ptil(idof)*sin(omec(idof)*dt) !integrate dimensionless vector
                  ptil(idof) = -xtil(idof)*sin(omec(idof)*dt) + ptil(idof)*cos(omec(idof)*dt)
                  xtil(idof) = tmp1
               end do

               xc = qj(:,0) / sqrt(dble(P))
               pc = pqj(:,0) / sqrt(dble(P))
               avghes = sum(hes_rec,3)/(2*navg+1)
               eigvec = avghes
               call dsyev('V','L',n,eigvec,n,eigval,work,3*n-1,ierr)
               if(ierr/=0) then
                  write(*,*) 'Diagonalization fails!'
                  write(*,*) 'Error information:', ierr
                  stop
               end if
               do idof = 1, n
                  if(eigval(idof) > 0d0) then
                     omec(idof) = sqrt(eigval(idof))
                     bhw(idof) = hbar*omec(idof)*beta
                     qcor(idof) = hbar*omec(idof)*beta/2/tanh(hbar*omec(idof)*beta/2)
                  else
                     omec(idof) = 0d0
                     bhw(idof) = 0d0
                     qcor(idof) = 1d0
                  end if
               end do
               call dgemv('N',n,n,1d0,eigvec_old,n,xtil,1,0d0,xtmp,1)
               call dgemv('T',n,n,1d0,eigvec,n,xtmp,1,0d0,xtil,1)
               call dgemv('N',n,n,1d0,eigvec_old,n,ptil,1,0d0,ptmp,1)
               call dgemv('T',n,n,1d0,eigvec,n,ptmp,1,0d0,ptil,1)!rescaling dimensionless vector
               eigvec_old = eigvec

               do idof = 1, natom*3
                  if(bhw(idof) > epsil) then
                     xtmp(idof) = xtil(idof) * sqrt(2*(qcor(idof)-1d0)/beta/omec(idof)**2)
                  else
                     xtmp(idof) = xtil(idof) * &
                        sqrt(2*hbar**2*beta*(1d0/12d0-bhw(idof)**2/720d0+bhw(idof)**4/30240d0-bhw(idof)**6/1209600d0))
                  end if
                  ptmp(idof) = sqrt(2*(qcor(idof)-1d0)/beta)*ptil(idof)!ptmp of FK-QCW(1)
               end do
               call dgemv('N',n,n,1d0,eigvec,n,xtmp,1,0d0,xphy,1)
               xphy = xc + xphy / sqrt(m)
               call dgemv('N',n,n,1d0,eigvec,n,ptmp,1,0d0,pphy,1)
               pphy = pc + pphy * sqrt(m)!and the pphy of FK-QCW(1)

               do idof = 1, n
                  tmpmat(:,idof) = eigvec(:,idof) / qcor(idof)
               end do
               call dgemm('N','T',n,n,n,1d0,tmpmat,n,eigvec,n,0d0,mtherm,n)
               do idof = 1, n
                  do jdof = 1, n
                     mtherm(idof,jdof) = mtherm(idof,jdof) / sqrt(m(idof)*m(jdof))!inverse of thermal mass Mtherm=M^-1/2*U*Q(u)^-1*U^T*M^-1/2
                  end do
               end do
               do idof = 1, natom
                  mtherm_diag(:,:,idof) = mtherm(3*idof-2:3*idof,3*idof-2:3*idof)
               end do
               call dsymv('U',n,1d0,mtherm,n,pphy,1,0d0,vtherm,1)!vtherm=fbeta=Mtherm^-1 p

               if(ftrj) then
                  write(cc,*) n*8
                  write(21) xphy, pphy, xc, pc,&
                     vtherm, mtherm_diag
               end if

            end if

            if(mod(i,nout).eq.0 .and. i.ge.nout*(2*navg+1)) then

               fxc = sum(fxc_rec,2)/(2*navg+1)
               do idof = 1, n
                  if(eigval(idof) > 0d0) then
                     omec_map(idof) = sqrt(eigval(idof))
                     bhw(idof) = beta*hbar*omec_map(idof)
                     qcor(idof) = bhw(idof)/2/tanh(bhw(idof)/2)
                  else
                     omec_map(idof) = -sqrt(-eigval(idof))
                     bhw(idof) = beta*hbar*omec_map(idof)
                     qcor(idof) = tanh(bhw(idof)/2)/(bhw(idof)/2)
                  end if
               end do
               !write(123,*) omec_map

               call dgemv('T',n,n,1d0,eigvec,n,fxc/sqrt(m),1,0d0,xtmp,1)
               do idof = 1, n
                  if(bhw(idof) > epsil) then
                     xtmp(idof) = xtmp(idof) * (sqrt(qcor(idof))-1d0)/omec_map(idof)**2
                  else if(bhw(idof) > 0d0) then
                     xtmp(idof) = xtmp(idof) * (1d0/24d0-bhw(idof)**2/640d0+79d0*bhw(idof)**4/967680d0)*hbar**2*beta**2
                  else if(bhw(idof) > -epsil) then
                     xtmp(idof) = xtmp(idof) * (1d0/24d0-19d0*bhw(idof)**2/5760d0+55d0*bhw(idof)**4/193536d0)*hbar**2*beta**2
                  else
                     xtmp(idof) = xtmp(idof) * (sqrt(qcor(idof))-1d0)/-omec_map(idof)**2
                  end if
               end do
               call dgemv('N',n,n,1d0,eigvec,n,xtmp,1,0d0,xphy,1)
               xphy = xphy / sqrt(m) + xc
               call dgemv('T',n,n,1d0,eigvec,n,pc/sqrt(m),1,0d0,ptmp,1)
               call dgemv('N',n,n,1d0,eigvec,n,ptmp*sqrt(qcor),1,0d0,pphy,1)
               pphy = pphy * sqrt(m)

               do idof = 1, n
                  tmpmat(:,idof) = eigvec(:,idof) / qcor(idof)
               end do
               call dgemm('N','T',n,n,n,1d0,tmpmat,n,eigvec,n,0d0,mtherm,n)
               do idof = 1, n
                  do jdof = 1, n
                     mtherm(idof,jdof) = mtherm(idof,jdof) / sqrt(m(idof)*m(jdof))
                  end do
               end do
               do idof = 1, natom
                  mtherm_diag(:,:,idof) = mtherm(3*idof-2:3*idof,3*idof-2:3*idof)
               end do
               call dsymv('U',n,1d0,mtherm,n,pphy,1,0d0,vtherm,1)

               if(ftrj) then
                  write(cc,*) n*8
                  write(22) xphy, pphy, xc, pc,&
                     vtherm, mtherm_diag
               end if

            end if

            if (mod(i,nprint) .eq. 0) then
               write(*,'(i,a,i10,a,i10)') idproc, ' Step:', i, '/', sampstep
               !write(cc,*) n
               !cc = '(i,a,'//trim(adjustl(cc))//'F)'
               !write(*,cc) idproc, ' xc', xc
               !write(*,cc) idproc, ' pc', pc
               !write(*,cc) idproc, ' fxc', fxc
               !vc = 0d0
               !summ = 0d0
               !do idof = 1, n / 3
               !   vc = vc + pqj(idof*3-2:idof*3,0)
               !   summ = summ + m(idof*3)
               !end do
               !vc = vc / summ
               !lc = 0d0
               !do idof = 1, n / 3
               !   lc = lc + cross_product(qj(idof*3-2:idof*3,0),pqj(idof*3-2:idof*3,0))
               !end do
               !write(*,'(i,a,3f)') idproc, ' vel COM', vc
               !write(*,'(i,a,3f)') idproc, ' ang mom COM:', lc
               CALL GET_TIME(dble(sampstep+eqstep)*samples,dble(sampstep+eqstep)*(k-1)+dble(i+eqstep))
            end if
         end if
      end do
      if(zmaster) then
         if(fcent) close(11)
         if(ftrj) close(21)
         if(ftrj) close(22)
         !close(123)
      end if

   end do

   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
   CALL MPI_FINALIZE(ierr)


contains


   ! =====================================================================================

   function utwo(x)
      !    ---------------------------------------------
      !    pair-wise potential form in term of distances
      !    ---------------------------------------------
      implicit none
      real*8::x, utwo

      if (x .le. rc)then
         utwo = exp(ah-bh*x-rh*x**2)-  &
            (c6/x**6+c8/x**8-c9/x**9+c10/x**10)*exp(-(rc/x-1.0d0)**2)
      else
         utwo = exp(ah-bh*x-rh*x**2)-(c6/x**6+c8/x**8-c9/x**9+c10/x**10)
      end if

      return
   End function utwo

   ! =========================================================

   function dutwo(x)
      !    ------------------------------------------------
      !    potential first derivative in term of distances
      !    ------------------------------------------------
      implicit none
      real*8::x, dutwo

      if (x .le. rc)then
         dutwo = exp(ah-bh*x-rh*x**2)*(-bh-2.0d0*rh*x) +      &
            (6.0d0*c6/x**7+8.0d0*c8/x**9-9.0d0*c9/x**10  &
            +10.0d0*c10/x**11)*exp(-(rc/x-1.0d0)**2) -      &
            2.0d0*exp(-(rc/x-1.0d0)**2)*rc*(rc/x-1.0d0)/x**2*   &
            (c6/x**6+c8/x**8-c9/x**9+c10/x**10)
      else
         dutwo = exp(ah-bh*x-rh*x**2)*(-bh-2.0d0*rh*x) +      &
            (6.0d0*c6/x**7+8.0d0*c8/x**9-9.0d0*c9/x**10  &
            +10.0d0*c10/x**11)
      end if

      return
   End function dutwo

   ! ==========================================================
   ! ------------------------------------------------------------------
   !
   subroutine inverse_stage(ubd, xbd)
      !     get the matrix to transform ubd to xbd
      implicit none
      real*8,intent(in):: ubd(n, P)
      real*8,intent(out):: xbd(n, P)
      integer:: i, j
      !
      xbd(:, 1) = ubd(:, 1)
      xbd(:, P) = ubd(:, P) + ubd(:, 1)
      do i = P-1, 2, -1
         xbd(:, i) = ubd(:, i) + dble(i-1)/dble(i)*xbd(:, i+1) + ubd(:, 1)/dble(i)
      end do
      !
      return
   end subroutine inverse_stage
   !
   ! -------------------------------------------------------------------------------------------------------------------
   ! ==============================================================================
   !
   subroutine getforce(xps, force_xps, Etot_xps)
      !    -----------------------------------------------
      !    pair-wise force on each particle
      !    -----------------------------------------------
      !
      implicit none
      real*8,intent(in):: xps(n)
      real*8,intent(out):: force_xps(n), Etot_xps
      integer:: i,j
      real*8:: dx,dy,dz
      real*8:: fx,fy,fz
      real*8:: rsq, dut, inv_rsq
      !
      Etot_xps = 0.0d0 
      force_xps = 0.0d0
      !
      do i = 1,natom
         do j = i+1,natom
            dx = xps(3*i-2) - xps(3*j-2)
            dy = xps(3*i-1) - xps(3*j-1)
            dz = xps(3*i)   - xps(3*j)
            dx = dx - box*anint(dx*inv_box)
            dy = dy - box*anint(dy*inv_box)
            dz = dz - box*anint(dz*inv_box)

            rsq = sqrt(dx*dx+dy*dy+dz*dz)
            if(rsq .le. rcut)then
               dut = dutwo(rsq)
               Etot_xps = Etot_xps + utwo(rsq)
            else
               dut = 0.0d0
            endif
            !
            fx = dx/rsq*dut
            fy = dy/rsq*dut
            fz = dz/rsq*dut
            force_xps(3*i-2) = force_xps(3*i-2) + fx
            force_xps(3*i-1) = force_xps(3*i-1) + fy
            force_xps(3*i)   = force_xps(3*i) + fz
            force_xps(3*j-2) = force_xps(3*j-2) - fx
            force_xps(3*j-1) = force_xps(3*j-1) - fy
            force_xps(3*j)   = force_xps(3*j) - fz
         end do
      end do
      !
      return
   end subroutine getforce
   !
   ! ===========================================================================================


   subroutine init_cmat(cmat)
      real*8, intent(inout) :: cmat(P,0:P-1)

      if(mod(P,2)==0)then
         cmat(:,0)=sqrt(1d0/P)
         do k=1, P/2-1
            do j=1, P
               cmat(j,k)=sqrt(2d0/P)*cos(2*pi*j*k/P)
            end do
         end do
         do j=1, P
            cmat(j,P/2)=sqrt(1d0/P)*(-1)**j
         end do
         do k=P/2+1, P-1
            do j=1, P
               cmat(j,k)=sqrt(2d0/P)*sin(2*pi*j*k/P)
            end do
         end do
      else
         cmat(:,0)=sqrt(1d0/P)
         do k=1, (P-1)/2
            do j=1, P
               cmat(j,k)=sqrt(2d0/P)*cos(2*pi*j*k/P)
            end do
         end do
         do k=(P+1)/2, P-1
            do j=1, P
               cmat(j,k)=sqrt(2d0/P)*sin(2*pi*j*k/P)
            end do
         end do
      end if

   end subroutine init_cmat


   subroutine CONS_TRPMDCOR(dtime)
      implicit none
      real*8, intent(in) :: dtime
      do j = 1, P-1
         pqj(:,j) = pqj(:,j) - (fqj(:,j) + mqj(:,j)*omegaj(j)**2*qj(:,j)) * dtime/2
      end do

      !!pxj = matmul(pqj/mqj,transpose(cmat)) * mxj
      !call dgemm('N','T',n,P,P,1d0,pqj/mqj,n,cmat,P,0d0,pxj,n)
      !pxj = pxj * mxj
      !do j = 1, P
      !   call rot_trans_corr(mxj(:,j), xj(:,j), pxj(:,j))
      !end do
      !!pqj = matmul(pxj/mxj,cmat) * mqj
      !call dgemm('N','N',n,P,P,1d0,pxj/mxj,n,cmat,P,0d0,pqj,n)
      !pqj = pqj * mqj

      qj(:,1:P-1) = qj(:,1:P-1) + pqj(:,1:P-1) / mqj(:,1:P-1) * dtime/2

      c1 = exp(-lgam*dtime)
      c2 = sqrt(1d0-c1**2)
      do j = 1, P - 1
         do idof = 1, n
            call gasdev(eta)
            pqj(idof,j) = c1(j)*pqj(idof,j) + c2(j)*sqrt(mqj(idof,j)/betap)*eta
         end do
      end do

      !!xj = matmul(qj,transpose(cmat))
      !call dgemm('N','T',n,P,P,1d0,qj,n,cmat,P,0d0,xj,n)
      !!pxj = matmul(pqj/mqj,transpose(cmat)) * mxj
      !call dgemm('N','T',n,P,P,1d0,pqj/mqj,n,cmat,P,0d0,pxj,n)
      !pxj = pxj * mxj
      !do j = 1, P
      !   call rot_trans_corr(mxj(:,j), xj(:,j), pxj(:,j))
      !end do
      !!pqj = matmul(pxj/mxj,cmat) * mqj
      !call dgemm('N','N',n,P,P,1d0,pxj/mxj,n,cmat,P,0d0,pqj,n)
      !pqj = pqj * mqj

      qj(:,1:P-1) = qj(:,1:P-1) + pqj(:,1:P-1) / mqj(:,1:P-1) * dtime/2

      !xj = matmul(qj,transpose(cmat))
      call dgemm('N','T',n,P,P,1d0,qj,n,cmat,P,0d0,xj,n)
      do j = 1, P
         call getforce(xj(:,j), fxj(:,j), Ep(j) )
      end do
      !fqj = matmul(fxj,cmat)
      call dgemm('N','N',n,P,P,1d0,fxj,n,cmat,P,0d0,fqj,n)

      do j = 1, P-1
         pqj(:,j) = pqj(:,j) - (fqj(:,j) + mqj(:,j)*omegaj(j)**2*qj(:,j)) * dtime/2
      end do
   end subroutine CONS_TRPMDCOR


   subroutine NMPIMD(pqj,qj,xj,fqj,fxj) ! NMPIMD
      implicit none
      real*8, intent(inout) :: pqj(n,0:P-1), qj(n,0:P-1)
      real*8, intent(inout) :: xj(n,P), fqj(n,0:P-1), fxj(n,P)
      pqj(:,0) = pqj(:,0) - fqj(:,0) * dt/2
      do j = 1, P-1
         pqj(:,j) = pqj(:,j) - (fqj(:,j) + mqj(:,j)*omegaj(j)**2*qj(:,j)) * dt/2
      end do

      !!pxj = matmul(pqj/mqj,transpose(cmat)) * mxj
      !call dgemm('N','T',n,P,P,1d0,pqj/mqj,n,cmat,P,0d0,pxj,n)
      !pxj = pxj * mxj
      !do j = 1, P
      !   call rot_trans_corr(mxj(:,j), xj(:,j), pxj(:,j))
      !end do
      !!pqj = matmul(pxj/mxj,cmat) * mqj
      !call dgemm('N','N',n,P,P,1d0,pxj/mxj,n,cmat,P,0d0,pqj,n)
      !pqj = pqj * mqj

      qj = qj + pqj / mqj * dt/2

      c1 = exp(-lgam*dt)
      c2 = sqrt(1d0-c1**2)
      do j = 0, P - 1
         do idof = 1, n
            call gasdev(eta)
            pqj(idof,j) = c1(j)*pqj(idof,j) + c2(j)*sqrt(mqj(idof,j)/betap)*eta
         end do
      end do

      !!xj = matmul(qj,transpose(cmat))
      !call dgemm('N','T',n,P,P,1d0,qj,n,cmat,P,0d0,xj,n)
      !!pxj = matmul(pqj/mqj,transpose(cmat)) * mxj
      !call dgemm('N','T',n,P,P,1d0,pqj/mqj,n,cmat,P,0d0,pxj,n)
      !pxj = pxj * mxj
      !do j = 1, P
      !   call rot_trans_corr(mxj(:,j), xj(:,j), pxj(:,j))
      !end do
      !!pqj = matmul(pxj/mxj,cmat) * mqj
      !call dgemm('N','N',n,P,P,1d0,pxj/mxj,n,cmat,P,0d0,pqj,n)
      !pqj = pqj * mqj

      qj = qj + pqj / mqj * dt/2

      !xj = matmul(qj,transpose(cmat))
      call dgemm('N','T',n,P,P,1d0,qj,n,cmat,P,0d0,xj,n)
      do j = 1, P
         call getforce(xj(:,j), fxj(:,j), Ep(j) )
      end do
      !fqj = matmul(fxj,cmat)
      call dgemm('N','N',n,P,P,1d0,fxj,n,cmat,P,0d0,fqj,n)

      pqj(:,0) = pqj(:,0) - fqj(:,0) * dt/2
      do j = 1, P-1
         pqj(:,j) = pqj(:,j) - (fqj(:,j) + mqj(:,j)*omegaj(j)**2*qj(:,j)) * dt/2
      end do
   end subroutine NMPIMD


   subroutine PACMD(pqj,qj,xj,fqj,fxj) ! PACMD
      implicit none
      real*8, intent(inout) :: pqj(n,0:P-1), qj(n,0:P-1)
      real*8, intent(inout) :: xj(n,P), fqj(n,0:P-1), fxj(n,P)
      pqj(:,0) = pqj(:,0) - fqj(:,0) * dtr/2
      do j = 1, P-1
         pqj(:,j) = pqj(:,j) - (fqj(:,j) + mqj(:,j)*omegaj(j)**2*qj(:,j)) * dtr/2
      end do

      !!pxj = matmul(pqj/mqjad,transpose(cmat)) * mxj
      !call dgemm('N','T',n,P,P,1d0,pqj/mqjad,n,cmat,P,0d0,pxj,n)
      !pxj = pxj * mxj
      !do j = 1, P
      !   call rot_trans_corr(mxj(:,j), xj(:,j), pxj(:,j))
      !end do
      !!pqj = matmul(pxj/mxj,cmat) * mqjad
      !call dgemm('N','N',n,P,P,1d0,pxj/mxj,n,cmat,P,0d0,pqj,n)
      !pqj = pqj * mqjad

      qj(:,0) = qj(:,0) + pqj(:,0) / mqj(:,0) * dtr/2
      qj(:,1:P-1) = qj(:,1:P-1) + pqj(:,1:P-1) / (gamad*mqj(:,1:P-1)) * dtr/2

      c1ad = exp(-lgamad*dtr)
      c2ad = sqrt(1d0-c1ad**2)
      do j = 1, P - 1
         do idof = 1, n
            call gasdev(eta)
            pqj(idof,j) = c1ad(j)*pqj(idof,j) + c2ad(j)*sqrt(mqj(idof,j)*gamad/betap)*eta
         end do
      end do

      !!xj = matmul(qj,transpose(cmat))
      !call dgemm('N','T',n,P,P,1d0,qj,n,cmat,P,0d0,xj,n)
      !!pxj = matmul(pqj/mqjad,transpose(cmat)) * mxj
      !call dgemm('N','T',n,P,P,1d0,pqj/mqjad,n,cmat,P,0d0,pxj,n)
      !pxj = pxj * mxj
      !do j = 1, P
      !   call rot_trans_corr(mxj(:,j), xj(:,j), pxj(:,j))
      !end do
      !!pqj = matmul(pxj/mxj,cmat) * mqjad
      !call dgemm('N','N',n,P,P,1d0,pxj/mxj,n,cmat,P,0d0,pqj,n)
      !pqj = pqj * mqjad

      qj(:,0) = qj(:,0) + pqj(:,0) / mqj(:,0) * dtr/2
      qj(:,1:P-1) = qj(:,1:P-1) + pqj(:,1:P-1) / (gamad*mqj(:,1:P-1)) * dtr/2

      !xj = matmul(qj,transpose(cmat))
      call dgemm('N','T',n,P,P,1d0,qj,n,cmat,P,0d0,xj,n)
      do j = 1, P
         call getforce(xj(:,j), fxj(:,j), Ep(j) )
      end do
      !fqj = matmul(fxj,cmat)
      call dgemm('N','N',n,P,P,1d0,fxj,n,cmat,P,0d0,fqj,n)
      call rot_trans_corr(mqj(:,0), qj(:,0), pqj(:,0), fqj(:,0), 1)
      ! test angular momentum
      !write(*,*) 'angular momentum:', cross_product(qj(1:3,0),pqj(1:3,0))+cross_product(qj(4:6,0),pqj(4:6,0)) 

      pqj(:,0) = pqj(:,0) - fqj(:,0) * dtr/2
      do j = 1, P-1
         pqj(:,j) = pqj(:,j) - (fqj(:,j) + mqj(:,j)*omegaj(j)**2*qj(:,j)) * dtr/2
      end do
   end subroutine PACMD


   subroutine rot_trans_corr(m_,x_,p,dv,flag)!???
      real*8, intent(in) :: m_(:), x_(:)
      real*8, intent(inout) :: p(:)
      real*8, intent(inout), optional :: dv(:)
      integer, intent(in), optional :: flag
      real*8, allocatable :: m(:), x(:,:), xm(:,:)
      real*8, allocatable :: v(:,:), vc2(:,:)
      real*8, allocatable :: f(:,:), fc2(:,:)
      real*8 :: cen(3), vc1(3), fc1(3), cmi(3,3), cam(3)
      real*8 :: cav(3), ct(3), caa(3), cmi_inv(3,3)
      integer :: iatom, natom, i1, i2, ipiv(3), info
      logical :: fv, ff

      if(present(dv)) then
         ff = .true.
      else
         ff = .false.
      end if
      if(present(flag).and.flag==0) then
         fv = .false.
      else
         fv = .true.         
      end if

      natom = size(x_) / 3
      allocate(m(natom))
      allocate(x(3,natom))
      allocate(xm(3,natom))
      if(fv) then
         allocate(v(3,natom))
         allocate(vc2(3,natom))
      end if
      if(ff) then
         allocate(f(3,natom))
         allocate(fc2(3,natom))
      end if

      do iatom = 1, natom
         m(iatom) = m_(iatom*3)
         x(:,iatom) = x_(iatom*3-2:iatom*3)
         if(fv) v(:,iatom) = p(iatom*3-2:iatom*3) / m_(iatom*3-2:iatom*3)
         if(ff) f(:,iatom) = -dv(iatom*3-2:iatom*3)
      end do

      cen = 0d0
      do iatom = 1, natom
         cen = cen + m(iatom)*x(:,iatom)
      end do
      cen = cen / sum(m) !calculate centroid of atoms(not beads!)
      do iatom = 1, natom
         xm(:,iatom) = x(:,iatom) - cen !reduced position
      end do

      cmi = 0d0 !??? inverse of center of mass?
      do iatom = 1, natom
         do i1 = 1, 3
            cmi(i1,i1) = cmi(i1,i1) + m(iatom)*sum(xm(:,iatom)**2)
         end do
         do i1 = 1, 3
            do i2 = 1, 3
               cmi(i1,i2) = cmi(i1,i2) - m(iatom)*xm(i1,iatom)*xm(i2,iatom)
            end do
         end do
      end do
      call pseudo_inv(cmi,cmi_inv,1d-10)

      if(fv) then
         vc1 = 0d0 !velocity of centroid(of atoms)
         do iatom = 1, natom
            vc1 = vc1 + m(iatom)*v(:,iatom)
         end do
         vc1 = vc1 / sum(m)

         cam = 0d0
         do iatom = 1, natom
            cam = cam + cross_product(xm(:,iatom),m(iatom)*(v(:,iatom)-vc1))
         end do

         cav = matmul(cmi_inv,cam)

         do iatom = 1, natom
            vc2(:,iatom) = cross_product(cav,xm(:,iatom))
            v(:,iatom) = v(:,iatom) - vc1 - vc2(:,iatom)
         end do

      end if

      if(ff) then
         fc1 = sum(f,2) / sum(m)

         ct = 0d0
         do iatom = 1, natom
            ct = ct + cross_product(xm(:,iatom),f(:,iatom)-fc1)
         end do

         caa = matmul(cmi_inv,ct)

         do iatom = 1, natom
            fc2(:,iatom) = m(iatom)*cross_product(caa,xm(:,iatom))
            f(:,iatom) = f(:,iatom) - fc1 - fc2(:,iatom)
         end do
      end if

      do iatom = 1, natom
         if(fv) p(iatom*3-2:iatom*3) = v(:,iatom) * m(iatom)
         if(ff) dv(iatom*3-2:iatom*3) = -f(:,iatom)
      end do

      deallocate(m)
      deallocate(x)
      deallocate(xm)
      if(fv) then
         deallocate(v)
         deallocate(vc2)
      end if
      if(ff) then
         deallocate(f)
         deallocate(fc2)
      end if

   end subroutine rot_trans_corr


   function cross_product(vec1,vec2)
      real*8 :: vec1(3),vec2(3),cross_product(3)
      cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
      cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
      cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
   end function cross_product


   subroutine pseudo_inv(a,b,eps)
      real*8, intent(in) :: a(:,:), eps
      real*8, intent(out) :: b(:,:)
      real*8, allocatable :: acp(:,:), s(:), u(:,:), vt(:,:), work(:)
      integer :: m, n, mn, mx, lwork, info, i

      m = size(a,1)
      n = size(a,2)
      mn = minval([m, n])
      mx = maxval([m, n])
      lwork = 5*mn + mx

      allocate(acp(m,n))
      allocate(s(mn))
      allocate(u(m,m))
      allocate(vt(n,n))
      allocate(work(lwork))

      acp = a
      call dgesvd('A','A',m,n,acp,m,s,u,m,vt,n,work,lwork,info)
      if(info /= 0) then
         write(*,*) 'SVD decomposition fails.'
         write(*,*) 'Info:', info
      end if

      b = 0d0
      do i = 1, mn
         if(s(i) .gt. eps) then
            b(i,i) = 1d0 / s(i)
         else
            b(i,i) = 0d0
         end if
      end do
      b = matmul(matmul(transpose(vt),b),transpose(u))

      deallocate(acp)
      deallocate(s)
      deallocate(u)
      deallocate(vt)
      deallocate(work)

   end subroutine pseudo_inv


   subroutine FFT(ain,nsiz,aout)
      implicit none
      complex*16, intent(in) :: ain(:)
      integer, intent(in) :: nsiz
      complex*16, intent(out) :: aout(:)
      real*8, parameter :: pi = acos(-1d0)
      complex*16, allocatable :: vtmp1(:), vtmp2(:)
      complex*16 :: cj=(0.0d0,1.0d0)
      integer :: i1, i2, i3, nsamp

      nsamp = 2**nsiz
      allocate(vtmp1(0:nsamp-1))
      allocate(vtmp2(0:nsamp-1))
      vtmp1 = ain(1:nsamp)

      do i1 = 0, nsiz-1 ! fast Fourier transform
         do i2 = 0, 2**i1 - 1
            do i3 = 0, nsamp/(2**i1)/2 - 1
               vtmp2(i2+2*i3*2**i1) = vtmp1(i2+i3*2**i1) + (cos(2*pi*i2/(2**i1*2)) - cj*sin(2*pi*i2/(2**i1*2))) * vtmp1(i2+i3*2**i1+nsamp/2)
               vtmp2(i2+(2*i3+1)*2**i1) = vtmp1(i2+i3*2**i1) - (cos(2*pi*i2/(2**i1*2)) - cj*sin(2*pi*i2/(2**i1*2))) * vtmp1(i2+i3*2**i1+nsamp/2)
            end do
         end do
         vtmp1 = vtmp2
      end do

      aout(1:nsamp) = vtmp2(:)
      deallocate(vtmp1)
      deallocate(vtmp2)

   end subroutine FFT


   SUBROUTINE GET_TIME(ttot,tnow)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: ttot
      REAL*8, INTENT(IN) :: tnow
      REAL*8 :: ct, ct2
      INTEGER :: time_d, time_h, time_m, time_s, time_ms

      CALL CPU_TIME( ct2 )

      ct=ct2-ct1 !CPU TIME (s)

      time_s=dint(ct)
      time_m=time_s/60
      time_h=time_m/60
      time_d=time_h/24

      time_ms=dint( (ct-time_s)*1d3 )
      time_s=mod(time_s,60)
      time_m=mod(time_m,60)
      time_h=mod(time_h,24)

      WRITE(*,*) "Simulation Time" 
      WRITE(*,"(a,i5,a,i5,a,i5,a,i5,a,i5,a)") &
         & " Time",time_d," d",time_h," h",time_m," m",time_s," s",time_ms," ms" 

      ct = ct2 - ct1 !CPU TIME (s)  
      ct = ct * ( ttot - tnow ) / tnow
      time_s = dint(ct)
      time_m = time_s/60
      time_h = time_m/60
      time_d = time_h/24

      time_ms=dint( (ct-time_s)*1d3 )
      time_s=mod(time_s,60)
      time_m=mod(time_m,60)
      time_h=mod(time_h,24)

      WRITE(*,*) "Remaining Time"
      WRITE(*,"(a,i5,a,i5,a,i5,a,i5,a,i5,a)") &
         & " Time",time_d," d",time_h," h",time_m," m",time_s," s",time_ms," ms"  

      ct = ct2 - ct1 !CPU TIME (s) 
      ct = ct * ttot / tnow
      time_s = dint(ct)
      time_m = time_s/60
      time_h = time_m/60
      time_d = time_h/24

      time_ms=dint( (ct-time_s)*1d3 )
      time_s=mod(time_s,60)
      time_m=mod(time_m,60)
      time_h=mod(time_h,24)

      WRITE(*,*) "Total Time"
      WRITE(*,"(a,i5,a,i5,a,i5,a,i5,a,i5,a)") &
         & " Time",time_d," d",time_h," h",time_m," m",time_s," s",time_ms," ms"  


      WRITE(*,*)

   END SUBROUTINE GET_TIME


end program main
