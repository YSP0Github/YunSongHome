
program test_Radon
    implicit none
    
    integer*4            :: nt, nx, nq, flag, k, i, j
    real*8               :: noise, dt
    real*8               :: dq, dx, RR(2,10)
    real*8               :: A(3,3),B(2,3),C(3,3),C1(2,2)
    real*8,allocatable   :: d(:,:), x(:,:), q(:,:),m(:,:)
    complex*16            :: R(10,10), toep(10,10), ufpp(10,40),C2(2,2),R1(80,69)
    
    
    nt=1001; nx=69; nq=70; noise=0.1; flag=1
    allocate(d(nt,nx), x(1,nx), q(1,nq),m(nt,nq))
    dt=0.01; dq=0.00001; dx=30
    d=0;q=0;x=0;m=0
    do i=1,nx
        x(1,i)=(i-1)*dx+30 !偏移距或地震道的位置
    end do
    do i=1,nq
        q(1,i)=(i-1)*dq ! q参数,q=射线参数(flag=1 )或抛物线曲率(flag=2)
    end do
    !
    !M(40, nq/3)=1; M(40, nq/3-1)=0.5; M(40, nq/3+1)=0.5
    !M(44,2*nq/3)=1; M(44,2*nq/3-1)=0.5; M(44,2*nq/3+1)=0.5
    !call Inv_Radon(nt, nx, nq, m, dt, x, q, flag, d)
    !open(10,file="t-x.txt")
    !do i=1,nt
    !    write(10,"(30f9.5)")   (d(i,j),j=1,nx)
    !end do
    !close(10)

    open(10,file="test3.dat")
    do j=1,nx
        do i=1,nt
            read(10,*) d(i,j)
        end do
    end do
    close(10)
    open(10,file="t-x.txt")
    do i=1,nt
        write(10,"(100f9.5)")   (d(i,j),j=1,nx)
    end do
    close(10)
    !Radon变换 t-x --> t-p
    call For_Radon(nt, nx, nq, d, dt, x, q, flag, noise, m)
    open(11,file="t-p.txt")
    do i=1,nt
        write(11,"(100f)")   (m(i,j),j=1,nq)
    end do    
    close(11)
    !Radon反变换 t-p --> t-x
    call Inv_Radon(nt, nx, nq, m, dt, x, q, flag, d)
    open(10,file="111t-x.txt")
    do i=1,nt
        write(10,"(100f9.5)")   (d(i,j),j=1,nx)
    end do
    close(10)
    

    deallocate(d,x,q,m)

 end program test_Radon
    