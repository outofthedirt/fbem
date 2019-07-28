module multi_green
    
    implicit none
    
    !private
    !public :: init, normalize
    
    ! Number of layers
    integer :: Nl
    ! Arrays containing lamé constants for each layer and density + poisson ratio
    ! All values are not normalized
    complex*16 :: mu(100), lb(100), rho(100), nu(100)
    real*8 :: zi(100), hi(100)
    ! Normalization refs
    real*8 :: rho0, mu0, k0, w
    ! Normalized arrays
    complex*16 :: mj(100), lj(100), rhj(100), nuj(100), aj(100), bj(100)
    real*8 :: zj(100), hj(100)
    
 contains
    
    subroutine init(N,Vs,Vp,xis,xip,dens,hl,Nref)
        
        implicit none
        ! Number of layers
        integer :: N
        ! Arrays containing lamé constants for each layer and density + poisson ratio
        ! All values are not normalized
        real*8 :: Vs(:), Vp(:), xis(:), xip(:), dens(:), hl(:)
        ! Reference layer for Normalization
        integer :: Nref
        ! Loop integer
        integer :: i
        ! Temps vars
        real*8 :: G, M, lev
        complex*16 :: mjj, ljj, nujj
        ! Define imag number
        complex*16 :: j
        j = (0.0D+00,1.0D+00)
        lev = 0.0D+00
        
        ! Fill all values
        do i=1,N
            ! Compute modulus from waves
            G = dens(i)*Vs(i)*Vs(i)
            M = dens(i)*Vp(i)*Vp(i)
            ! Compute complex lamé values
            mjj = G*(1 + 2*j*xis(i))
            ljj = (M - 2*G)*(1 + 2*j*xip(i)) + 4*j*G*(xip(i) - xis(i))
            ! Compute depth
            lev = lev + hl(i)
            ! Fill arrays
            mu(i) = mjj
            lb(i) = ljj
            rho(i) = dens(i)
            nu(i) = 0.5D+00*ljj/(ljj + mjj)
            zi(i) = lev
            hi(i) = hl(i)
        end do
        
        ! Reference values
        Nl = N
        rho0 = dens(Nref)
        mu0 = rho0*Vs(Nref)*Vs(Nref)
        
        return
        
    end subroutine
    
    subroutine normalize(freq)
        
        implicit none
        ! excitation frequency
        real*8 :: freq
        
        ! excitation pulse
        w = 8*atan(1.0D+00)*freq
        
        ! compute 
        k0 = w*sqrt(rho0/mu0)
        
        ! normalize
        zj = k0*zi
        hj = k0*hi
        mj = mu/mu0
        lj = lb/mu0
        nuj = nu
        rhj = rho/rho0
        
        return
        
    end subroutine
    
    function Fij(i,j,zt)
        
        implicit none
        ! layers no
        integer :: i,j
        ! transform param
        complex*16 :: zt
        ! output value
        complex*16 :: Fij
        
        ! compute
        Fij = 2*zt*(4*zt*zt*(mj(j)-mj(i))*(mj(j)-mj(i)) * (zt*zt - aj(i)*bj(i)) &
            + 2*(mj(j)-mj(i))*(2*rhj(i)*zt*zt - rhj(j)*(zt*zt - aj(i)*bj(i))) &
            - rhj(i)*(rhj(j)-rhj(i)))
        
        return
    end function Fij
    
    function Gij(i,j,zt,k1)
        
        implicit none
        ! layers no
        integer :: i,j
        ! transfer param
        complex*16 :: zt
        ! param value
        real*8 :: k1
        ! output value
        complex*16 :: Gij
        
        ! compute
        Gij = 2*zt*zt*(mj(j)-mj(i))*((1+k1)*(bj(j)-bj(i))+(1-k1)*(aj(j)-aj(i))) &
            + (1+k1)*(rhj(i)*bj(j)+rhj(j)*bj(i)) + (1-k1)*(rhj(i)*aj(j)+rhj(j)*aj(i))
            
        return
    end function Gij
    
    function Hij(i,j,zt)
        
        implicit none
        ! layers no
        integer :: i,j
        ! transform param
        complex*16 :: zt
        ! output value
        complex*16 :: Hij
        
        ! compute
        Hij = 2*zt*(2*(mj(j)-mj(i))*(zt*zt* aj(j)*bj(i))-(rhj(j)-rhj(i)))
        
        return
    end function Hij
    
    function Sij(i,j,zt,k1,k2,k3)
        
        implicit none
        ! layers no
        integer :: i,j
        ! transfer param
        complex*16 :: zt
        ! param value
        real*8 :: k1,k2,k3
        ! output value
        complex*16 :: Sij
        
        ! compute
        Sij = 4*zt*zt*(mj(j)-mj(i))*(mj(j)-mj(i))*(zt*zt-k1*aj(i)*bj(i))*(zt*zt-k2*aj(j)*bj(j)) &
            + 4*zt*zt*(mj(j)-mj(i))*(rhj(i)*(zt*zt-k2*aj(j)*bj(j))-rhj(j)*(zt*zt-k1*aj(i)*bj(i))) &
            + zt*zt*(rhj(j)-rhj(i))*(rhj(j)-rhj(i))-k2*(rhj(i)*bj(j)+k2*k3*rhj(j)*bj(i)) &
            * (rhj(i)*aj(j)+k1*k3*rhj(j)*aj(i))
        
        return
    end function Sij
    
    function Hstep(x)
        
        implicit none
        ! input value
        real*8 :: x
        ! small tol value
        real*8 :: tol=1e-9
        ! output
        real*8 :: Hstep
        
        if (abs(x) < tol) then
            Hstep = 0.5D+00
        else if (x > 0.0D+00) then
            Hstep = 1.0D+00
        else
            Hstep = 0.0D+00
        end if
        
        return
    end function Hstep
    
    pure function matinv2(A) result(B)
        !! Performs a direct calculation of the inverse of a 2×2 matrix.
        complex*16, intent(in) :: A(2,2)   !! Matrix
        complex*16 :: B(2,2)   !! Inverse matrix
        complex*16 :: detinv
        
        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
        
        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * A(2,2)
        B(2,1) = -detinv * A(2,1)
        B(1,2) = -detinv * A(1,2)
        B(2,2) = +detinv * A(1,1)
        
        return
    end function matinv2
    
    subroutine wjdmwjum(zt,ss,zz,rr,th0,wj)
        
        implicit none
        ! transfom parameter
        complex*16,intent(in) :: zt
        ! position of load
        real*8,intent(in) :: ss
        ! position of observation point
        real*8,intent(in) :: zz,rr,th0
        ! output array
        complex*16,intent(out) :: wj(6)
        ! loop integer
        integer :: i
        ! intermediate variables
        complex*16 :: ks1
        real*8 :: s, z, r
        integer :: Ns, Nz
        ! intermediate matrixes
        complex*16 :: Ehj(2,2),Ru0(2,2),Rdn(2,2),Ruj(2,2),Tuj(2,2),Rdj(2,2),Tdj(2,2),Id(2,2)
        complex*16 :: Eh(Nl,2,2),Rj(Nl,2,2),Tj(Nl,2,2)
        
        ! for testing only
        Rj = (0.0D+00,0.0D+00)
        Tj = (0.0D+00,0.0D+00)
        
        ! Identity matrix
        Id(1,1) = (1.0D+00,0.0D+00)
        Id(2,2) = Id(1,1)
        Id(1,2) = (0.0D+00,0.0D+00)
        Id(2,1) = Id(1,2)
        
        ! Normalize
        s = k0*ss
        r = k0*rr
        z = k0*zz
        
        ! Start with surface layer
        i = 1
        ! wave numbers
        aj(i) = sqrt(zt*zt-rhj(i)/(lj(i)+2*mj(i)))
        bj(i) = sqrt(zt*zt-rhj(i)/mj(i))
        aj(i+1) = sqrt(zt*zt-rhj(i+1)/(lj(i+1)+2*mj(i+1)))
        bj(i+1) = sqrt(zt*zt-rhj(i+1)/mj(i+1))
        ! kpj = sqrt(rhj(i)/(lj(i)+2*mj(i)))
        ks1 = rhj(i)/mj(i)
        ! surface reflection matrix
        Ru0(1,1) = (bj(i)*bj(i)+zt*zt)**2 + 4*zt*zt*aj(i)*bj(i)
        Ru0(1,2) = 4*zt*bj(i)*(bj(i)*bj(i)+zt*zt)
        Ru0(2,1) = 4*zt*aj(i)*(bj(i)*bj(i)+zt*zt)
        Ru0(2,2) = Ru0(1,1)
        Ru0 = Ru0/((2*zt*zt-ks1)**2 - 4*zt*zt*aj(i)*bj(i))
        ! exponation matrix
        Ehj(1,1) = exp(-aj(i)*zj(i))
        Ehj(2,1) = (0.0D+00,0.0D+00)
        Ehj(1,2) = Ehj(2,1)
        Ehj(2,2) = exp(-bj(i)*zj(i))
        Eh(i,:,:) = Ehj
        ! Local reflection matrix for upward waves
        Ruj(1,1) = Sij(i,i+1,zt,1.0D+00,-1.0D+00,-1.0D+00)
        Ruj(1,2) = bj(i+1)*Fij(i,i+1,zt)
        Ruj(2,1) = aj(i+1)*Fij(i,i+1,zt)
        Ruj(2,2) = Sij(i,i+1,zt,1.0D+00,-1.0D+00,1.0D+00)
        Ruj = -Ruj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)      
        ! Local transmission matrix for upward waves
        Tuj(1,1) = aj(i+1)*Gij(i,i+1,zt,1.0D+00)
        Tuj(1,2) = bj(i+1)*Hij(i,i+1,zt)
        Tuj(2,1) = -aj(i+1)*Hij(i+1,i,zt)
        Tuj(2,2) = bj(i+1)*Gij(i,i+1,zt,-1.0D+00)
        Tuj = -rhj(i+1)*Tuj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)
        ! Local reflection matrix for downward waves
        Rdj(1,1) = Sij(i,i+1,zt,-1.0D+00,1.0D+00,1.0D+00)
        Rdj(1,2) = -bj(i)*Fij(i+1,i,zt)
        Rdj(2,1) = -aj(i)*Fij(i+1,i,zt)
        Rdj(2,2) = Sij(i,i+1,zt,-1.0D+00,1.0D+00,-1.0D+00)
        Rdj = -Rdj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)      
        ! Local transmission matrix for downward waves
        Tdj(1,1) = aj(i)*Gij(i,i+1,zt,1.0D+00)
        Tdj(1,2) = -bj(i)*Hij(i+1,i,zt)
        Tdj(2,1) = aj(i)*Hij(i,i+1,zt)
        Tdj(2,2) = bj(i)*Gij(i,i+1,zt,-1.0D+00)
        Tdj = -rhj(i)*Tdj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)
        ! Generalized reflection matrix for upward waves
        Tj(i,:,:) = matmul(matinv2(Id - matmul(Rdj,matmul(Ehj,matmul(Ru0,Ehj)))),Tuj)
        Rj(i,:,:) = Ruj + matmul(Tdj,matmul(Ehj,matmul(Ru0,matmul(Ehj,Tj(i,:,:)))))!                 
        ! Check if obervation layer
        if (zj(i) > z) then
            Nz = i
        end if
        ! Increase i
        i = i + 1
        
        ! loop while not in loaded layer
        do while (zj(i) < s)
            
            ! wave numbers
            aj(i+1) = sqrt(zt*zt-rhj(i+1)/(lj(i+1)+2*mj(i+1)))
            bj(i+1) = sqrt(zt*zt-rhj(i+1)/mj(i+1))
            ! exponation matrix
            Ehj(1,1) = exp(-aj(i)*zj(i))
            Ehj(2,1) = (0.0D+00,0.0D+00)
            Ehj(1,2) = Ehj(2,1)
            Ehj(2,2) = exp(-bj(i)*zj(i))
            Eh(i,:,:) = Ehj
            ! Local reflection matrix for upward waves
            Ruj(1,1) = Sij(i,i+1,zt,1.0D+00,-1.0D+00,-1.0D+00)
            Ruj(1,2) = bj(i+1)*Fij(i,i+1,zt)
            Ruj(2,1) = aj(i+1)*Fij(i,i+1,zt)
            Ruj(2,2) = Sij(i,i+1,zt,1.0D+00,-1.0D+00,1.0D+00)
            Ruj = -Ruj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)      
            ! Local transmission matrix for upward waves
            Tuj(1,1) = aj(i+1)*Gij(i,i+1,zt,1.0D+00)
            Tuj(1,2) = bj(i+1)*Hij(i,i+1,zt)
            Tuj(2,1) = -aj(i+1)*Hij(i+1,i,zt)
            Tuj(2,2) = bj(i+1)*Gij(i,i+1,zt,-1.0D+00)
            Tuj = -rhj(i+1)*Tuj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)
            ! Local reflection matrix for downward waves
            Rdj(1,1) = Sij(i,i+1,zt,-1.0D+00,1.0D+00,1.0D+00)
            Rdj(1,2) = -bj(i)*Fij(i+1,i,zt)
            Rdj(2,1) = -aj(i)*Fij(i+1,i,zt)
            Rdj(2,2) = Sij(i,i+1,zt,-1.0D+00,1.0D+00,-1.0D+00)
            Rdj = -Rdj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)      
            ! Local transmission matrix for downward waves
            Tdj(1,1) = aj(i)*Gij(i,i+1,zt,1.0D+00)
            Tdj(1,2) = -bj(i)*Hij(i+1,i,zt)
            Tdj(2,1) = aj(i)*Hij(i,i+1,zt)
            Tdj(2,2) = bj(i)*Gij(i,i+1,zt,-1.0D+00)
            Tdj = -rhj(i)*Tdj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)
            ! Generalized reflection matrix for upward waves
            Tj(i,:,:) = matmul(matinv2(Id - matmul(Rdj,matmul(Ehj,matmul(Rj(i-1,:,:),Ehj)))),Tuj)
            Rj(i,:,:) = Ruj + matmul(Tdj,matmul(Ehj,matmul(Rj(i-1,:,:),matmul(Ehj,Tj(i,:,:)))))
            ! Check if obervation layer
            if (zj(i) > z) then
                Nz = i
            end if
            ! increase layer
            i = i + 1
            
        end do
        
        ! Save loaded layer
        Ns = i
        
        ! Continue with bottom layer
        i = Nl - 1
        ! wave numbers
        aj(i) = sqrt(zt*zt-rhj(i)/(lj(i)+2*mj(i)))
        bj(i) = sqrt(zt*zt-rhj(i)/mj(i))
        aj(i+1) = sqrt(zt*zt-rhj(i+1)/(lj(i+1)+2*mj(i+1)))
        bj(i+1) = sqrt(zt*zt-rhj(i+1)/mj(i+1))
        ! halfspace reflection matrix
        Rdn = (0.0D+00,0.0D+00)
        ! exponation matrix
        Ehj(1,1) = exp(-aj(i+1)*zj(i+1))
        Ehj(2,1) = (0.0D+00,0.0D+00)
        Ehj(1,2) = Ehj(2,1)
        Ehj(2,2) = exp(-bj(i+1)*zj(i+1))
        Eh(i+1,:,:) = Ehj
        ! Local reflection matrix for upward waves
        Ruj(1,1) = Sij(i,i+1,zt,1.0D+00,-1.0D+00,-1.0D+00)
        Ruj(1,2) = bj(i+1)*Fij(i,i+1,zt)
        Ruj(2,1) = aj(i+1)*Fij(i,i+1,zt)
        Ruj(2,2) = Sij(i,i+1,zt,1.0D+00,-1.0D+00,1.0D+00)
        Ruj = -Ruj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)      
        ! Local transmission matrix for upward waves
        Tuj(1,1) = aj(i+1)*Gij(i,i+1,zt,1.0D+00)
        Tuj(1,2) = bj(i+1)*Hij(i,i+1,zt)
        Tuj(2,1) = -aj(i+1)*Hij(i+1,i,zt)
        Tuj(2,2) = bj(i+1)*Gij(i,i+1,zt,-1.0D+00)
        Tuj = -rhj(i+1)*Tuj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)
        ! Local reflection matrix for downward waves
        Rdj(1,1) = Sij(i,i+1,zt,-1.0D+00,1.0D+00,1.0D+00)
        Rdj(1,2) = -bj(i)*Fij(i+1,i,zt)
        Rdj(2,1) = -aj(i)*Fij(i+1,i,zt)
        Rdj(2,2) = Sij(i,i+1,zt,-1.0D+00,1.0D+00,-1.0D+00)
        Rdj = -Rdj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)      
        ! Local transmission matrix for downward waves
        Tdj(1,1) = aj(i)*Gij(i,i+1,zt,1.0D+00)
        Tdj(1,2) = -bj(i)*Hij(i+1,i,zt)
        Tdj(2,1) = aj(i)*Hij(i,i+1,zt)
        Tdj(2,2) = bj(i)*Gij(i,i+1,zt,-1.0D+00)
        Tdj = -rhj(i)*Tdj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)
        ! Generalized reflection matrix for downward waves
        Tj(i,:,:) = matmul(matinv2(Id - matmul(Ruj,matmul(Ehj,matmul(Rdn,Ehj)))),Tdj)
        Rj(i,:,:) = Rdj + matmul(Tuj,matmul(Ehj,matmul(Rdn,matmul(Ehj,Tj(i,:,:)))))
        ! Check if obervation layer
        if (zj(i) <= z) then
            Nz = i+1
        end if
        ! Decrease i
        i = i - 1
        
        ! loop while not in loaded layer
        do while (zj(i) >= s)

            ! wave numbers
            aj(i) = sqrt(zt*zt-rhj(i)/(lj(i)+2*mj(i)))
            bj(i) = sqrt(zt*zt-rhj(i)/mj(i))
            ! exponation matrix
            Ehj(1,1) = exp(-aj(i+1)*zj(i+1))
            Ehj(2,1) = (0.0D+00,0.0D+00)
            Ehj(1,2) = Ehj(2,1)
            Ehj(2,2) = exp(-bj(i+1)*zj(i+1))
            Eh(i+1,:,:) = Ehj
            ! Local reflection matrix for upward waves
            Ruj(1,1) = Sij(i,i+1,zt,1.0D+00,-1.0D+00,-1.0D+00)
            Ruj(1,2) = bj(i+1)*Fij(i,i+1,zt)
            Ruj(2,1) = aj(i+1)*Fij(i,i+1,zt)
            Ruj(2,2) = Sij(i,i+1,zt,1.0D+00,-1.0D+00,1.0D+00)
            Ruj = -Ruj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)      
            ! Local transmission matrix for upward waves
            Tuj(1,1) = aj(i+1)*Gij(i,i+1,zt,1.0D+00)
            Tuj(1,2) = bj(i+1)*Hij(i,i+1,zt)
            Tuj(2,1) = -aj(i+1)*Hij(i+1,i,zt)
            Tuj(2,2) = bj(i+1)*Gij(i,i+1,zt,-1.0D+00)
            Tuj = -rhj(i+1)*Tuj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)
            ! Local reflection matrix for downward waves
            Rdj(1,1) = Sij(i,i+1,zt,-1.0D+00,1.0D+00,1.0D+00)
            Rdj(1,2) = -bj(i)*Fij(i+1,i,zt)
            Rdj(2,1) = -aj(i)*Fij(i+1,i,zt)
            Rdj(2,2) = Sij(i,i+1,zt,-1.0D+00,1.0D+00,-1.0D+00)
            Rdj = -Rdj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)      
            ! Local transmission matrix for downward waves
            Tdj(1,1) = aj(i)*Gij(i,i+1,zt,1.0D+00)
            Tdj(1,2) = -bj(i)*Hij(i+1,i,zt)
            Tdj(2,1) = aj(i)*Hij(i,i+1,zt)
            Tdj(2,2) = bj(i)*Gij(i,i+1,zt,-1.0D+00)
            Tdj = -rhj(i)*Tdj/Sij(i,i+1,zt,1.0D+00,1.0D+00,1.0D+00)
            ! Generalized reflection matrix for downward waves
            Tj(i,:,:) = matmul(matinv2(Id - matmul(Ruj,matmul(Ehj,matmul(Rj(i+1,:,:),Ehj)))),Tdj)
            Rj(i,:,:) = Rdj + matmul(Tuj,matmul(Ehj,matmul(Rj(i+1,:,:),matmul(Ehj,Tj(i,:,:)))))
            ! Decrease i
            i = i - 1
        
        end do
        
        !!!!! NEXT TO DO: Eh layer for loaded medium
        !                 Check in which layer observation point is
        !                 First try of vdla 
        !                 Check if Ns works for first layer (not the case)
        !                 
        
        do i=1,Nl
            print*,'   PRINTING FOR LAYER ',i
            print*,Rj(i,:,:)
            print*,Tj(i,:,:)
            print*,Eh(i,:,:)
        end do
        
        print*,'Loaded layer is ',Ns
        print*,'Obervation layer is ',Nz
        
    end subroutine
    
    
end module multi_green
    
    
