subroutine For_Radon(nt, nx, nq, d, dt, x, q, flag, noise, m)  ! t-x --> t-p
!% FOR－RADON 正 Radon 变换。 输入地震数据(剥去文件头、剥去道头,按道序排列的数据矩阵,每一列为一道的数据) 。
!% 此函数返回在 Radon 域中的数据。
!% m＝for_radon(d,dt,x,q,flag,noise)
!% d：输入地震数据(d[nt,nx ],nt 是时间采样点数,nx 是道数)
!% dt：时间采样间隔,以 s 为单位
!% x:偏移距或地震道的位置,以m为单位(x(nx))
!% q:射线参数(flag=1 )或抛物线曲率(flag=2)
!% flag:1 tau-p变换; 2 tau-q变换
!% noise:白噪系数,请勿设置为0
!% 输出m:变换到Radon域中的数据
    use operator_i
    implicit none 

    integer*4         :: nt, nx, nq, flag, k, i, j
    real*8            :: d(nt,nx), x(1,nx), q(1,nq),m(nt,nq)   
    real*8            :: noise, omega, dt
    real*8,allocatable :: d_real(:),d_image(:)
    complex*16,allocatable :: R(:,:), toep(:,:), ufpp(:,:), Matrix(:,:), g(:,:), uf(:,:), temp(:,:),RR(:,:)

    
    
    allocate(d_real(nt),d_image(nt))
    allocate(R(nq,nx), toep(nq,nx), ufpp(nt,nq), Matrix(nq,nq), g(nq,1), uf(nt,nx), temp(nq,nq),RR(nx,nq))
    do j=1,nx
        d_real=d(:,j);d_image=0.0
        CALL FFT(d_real,d_image,nt,2)
        uf(:,j)=cmplx(d_real(:),d_image(:))
    end do
    do k=1,nt
        if(k<=nt/2+1)then
            omega=6.28318530717959*(k-1)/(nt*dt)
            R=exp(matmul(cmplx(0,1)*omega*transpose(q),x**flag))
            RR=conjg(transpose(R))
            Matrix=matmul(R,RR)
            temp=Matrix+noise*eye(nq)
            toep=matmul(bcinv(temp),R)
            g=matmul(toep,reshape(uf(k:k,1:nx),(/nx,1/)))
            ufpp(k:k,:)=transpose(g)
        else
            ufpp(k,:)=conjg(ufpp(nt+2-k,:))
        end if
    end do
!    RR=conjg(transpose(R))
open(10,file="Matrix.txt")
    do i=1,nq
        write(10,*)   (Matrix(i,j),j=1,nq)
    end do
    close(10)
    
    
    
    do j=1,nq
        d_real(1:nt)=real(ufpp(:,j));d_image=imag(ufpp(:,j))
        CALL FFT(d_real,d_image,nt,1)
        m(:,j)=d_real(1:nt)
    end do
    deallocate(d_real,d_image)
    deallocate(R, toep, ufpp, Matrix, g, uf, temp,RR)
end subroutine For_Radon
    
    
    
subroutine Inv_Radon(nt, nx, nq, m, dt, x, q, flag, d) !t-= --> t-x
!% b=inv_radon(m, dt, x, q, flag)
!%m:输入Radon域中的数据(m[nt, nq], nt是时间采样点数, nq是q道的数目)
!% dt:时间采样间隔，以s为单位
!% x:偏移距地震道的位置,以m为单位(x(nx))
!% q:射线参数(flag=1 )或抛物线曲率(flag=2)
!% flag:1 tau-p变换; 2 tau-q变换
!%输出b:返回到x-t域中的地震数据(b[nt, nx])
    use operator_i
    implicit none 

    integer*4         :: nt, nx, nq, flag, k, i, j
    real*8            :: m(nt,nq), x(1,nx), q(1,nq),d(nt,nx)
    real*8            :: omega, dt
     real*8,allocatable :: m_real(:),m_image(:)
    complex*16,allocatable :: R(:,:), temp(:,:), fp(:,:), ufxx(:,:)

    
    allocate(m_real(nt),m_image(nt),R(nq,nx), temp(nx,1), fp(nt,nq), ufxx(nt,nx))
    temp=0;ufxx=0
    nt=size(m,1); nx=size(x,2); nq=size(m,2)
    do j=1,nq
        m_real=m(:,j);m_image=0.0
        CALL FFT(m_real,m_image,nt,2)
        fp(:,j)=cmplx(m_real(:),m_image(:))
    end do
    
    do k=1,nt
        if(k<nt/2+1)then
            omega=6.28318530717959*(k-1)/(nt*dt)
            R=exp(matmul(cmplx(0,1)*omega*transpose(q),x**flag))
            temp=matmul(conjg(transpose(R)),reshape(fp(k:k,1:nq),(/nq,1/)))
            ufxx(k:k,:)=transpose(temp)
        else
            ufxx(k,:)=conjg(ufxx(nt+2-k,:))
        end if
    end do

    do j=1,nx
        m_real(1:nt)=real(ufxx(:,j));m_image(1:nt)=imag(ufxx(:,j))
        CALL FFT(m_real,m_image,nt,1)
        d(:,j)=m_real(1:nt)
    end do

   deallocate(m_real,m_image,R, temp, fp, ufxx)
end subroutine Inv_Radon    
    
    
SUBROUTINE FFT(Xreal_,Ximage_,MM,nf)
    
!******************************************************************c
!                                                                  c
!  功能：复数组1-D快速Fourier变换
!
!  输入参数说明：
!      Freal(n): 输入数据实部          
!     Fimage(n): 输入数据虚部     
!             N: 点数（M 必须为 2 的幂次方）
!            NF: 正、反变换标志量(1:反变换;2:正变换)
!
!  输出参数说明：
!      Freal(n): 变换后频谱实部          
!     Fimage(n): 变换后频谱虚部     
!
!  对应频率分布说明：
!      数据Freal(n)和Fimage(n)对应的频率分布位置为:
!       0,1,.......,n/2-1,n/2,-(n/2-1),......,-1
!------------------------------------------------------------------c
    implicit none
    integer*4 n,nf,MM
    real*8 Xreal_(1:MM),Ximage_(1:MM)
    real*8,allocatable   :: Xreal(:),Ximage(:)
    integer*4 nu,n2,nu1,k,k1,k1n2,l,i,ibitr
    real*8 f,p,arg,c,s,treal,timage
    
    do i=1,100
        if(2**i>=MM)then
            n=2**i
            exit
        end if
    end do
    allocate(Xreal(1:n),Ximage(1:n))
    Xreal=0;Ximage=0
    Xreal(1:MM)=Xreal_(1:MM)
    Ximage(1:MM)=Ximage_(1:MM)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    nu=int(log(float(n))/0.693147+0.001)  !ln2=0.693147
	n2=n/2
	nu1=nu-1
	f=float((-1)**nf)
	k=0
    DO l=1,nu,1
      DO while (k.lt.n)
        do i=1,n2,1
            p=ibitr(k/2**nu1,nu)
            arg=6.2831853*p*f/float(n)
            c=cos(arg)
            s=sin(arg)
            k1=k+1
            k1n2=k1+n2
            treal=Xreal(k1n2)*c+Ximage(k1n2)*s
            timage=Ximage(k1n2)*c-Xreal(k1n2)*s
            Xreal(k1n2)=Xreal(k1)-treal
            Ximage(k1n2)=Ximage(k1)-timage
            Xreal(k1)=Xreal(k1)+treal
            Ximage(k1)=Ximage(k1)+timage
            k=k+1
        end do
        k=k+n2
      END DO
      k=0
      nu1=nu1-1
      n2=n2/2
    END DO

    DO k=1,n,1
      i=ibitr(k-1,nu)+1
      if(i.gt.k) then
        treal=Xreal(k)
        timage=Ximage(k)
        Xreal(k)=Xreal(i)
        Ximage(k)=Ximage(i)
        Xreal(i)=treal
        Ximage(i)=timage
      end if
    END DO
    
    IF(nf.ne.2) THEN
      do i=1,n,1
        Xreal(i)=Xreal(i)/float(n)
        Ximage(i)=Ximage(i)/float(n)
      end do
    END IF
    
    Xreal_(1:MM)=Xreal(1:MM)
    Ximage_(1:MM)=Ximage(1:MM)
    
    deallocate(Xreal,Ximage)
END SUBROUTINE FFT
        
FUNCTION IBITR(J,NU)
    implicit none
    integer*4 ibitr,j,nu
    integer*4 j1,itt,i,j2
	j1=j
	itt=0
    do i=1,nu,1
      j2=j1/2
      itt=itt*2+(j1-2*j2)
      ibitr=itt
      j1=j2
    end do
    END FUNCTION IBITR

module operator_i

    interface operator(.i.) !!// 自定义重载矩阵求逆操作符
    module procedure brinv !!// 实矩阵求逆
    module procedure bcinv !!// 复矩阵求逆
    end interface

contains

!!// 实矩阵求逆核心代码 trans from 徐士良《Fortran常用算法集》
!!// 作者：zuozhihua 时间：2020/8/18 地点：江西
function brinv(re) result(r)

    real*8,intent(in) :: re(:,:) !!// 原矩阵
    real*8 :: r(size(re,1),size(re,1)) !!// 逆矩阵
    integer*4 :: flag,n !!// flag判断奇异性；n是矩阵维度
    real*8 :: t,d !!// 中间变量
    integer*4 :: is(size(re,1)),js(size(re,1)) !!// 中间变量

    n=size(re,1) !!// 获取矩阵维度
    r=re !!// 复制值
    flag=1
    do k=1,n
        d=0.0
        do i=k,n
            do j=k,n
                if (abs(r(i,j)).gt.d) then
                    d=abs(r(i,j))
                    is(k)=i
                    js(k)=j
                end if
            end do
        end do
        if (d+1.0.eq.1.0) then
            flag=0
            write(*,*) "flag=0,实矩阵奇异！" 
            return !!// 矩阵奇异，退出逆矩阵求解
        end if
        do j=1,n
            t=r(k,j)
            r(k,j)=r(is(k),j)
            r(is(k),j)=t
        end do
        do i=1,n
            t=r(i,k)
            r(i,k)=r(i,js(k))
            r(i,js(k))=t
        end do
        r(k,k)=1/r(k,k)
        do j=1,n
            if (j.ne.k) r(k,j)=r(k,j)*r(k,k)
        end do
        do i=1,n
            if (i.ne.k) then
                do j=1,n
                    if (j.ne.k) r(i,j)=r(i,j)-r(i,k)*r(k,j)
                end do
            end if
        end do
        do i=1,n
            if (i.ne.k) r(i,k)=-r(i,k)*r(k,k)
        end do
    end do
    do k=n,1,-1
        do j=1,n
            t=r(k,j)
            r(k,j)=r(js(k),j)
            r(js(k),j)=t
        end do
        do i=1,n
            t=r(i,k)
            r(i,k)=r(i,is(k))
            r(i,is(k))=t
        end do
    end do
end function
!!// 复矩阵求逆核心代码 trans from 徐士良《Fortran常用算法集》
!!// 作者：zuozhihua 时间：2020/8/10 地点：江西
function bcinv(cpx)
    complex*16,intent(in) :: cpx(:,:) !!// 原矩阵
    complex*16 :: bcinv(size(cpx,1),size(cpx,2)) !!// 逆矩阵
    integer*4 :: flag,n !!// flag判断奇异性；n是矩阵维度
    real*8 :: ar(size(cpx,1),size(cpx,1)),ai(size(cpx,1),size(cpx,1)) !!// 实部矩阵ar；虚部矩阵ai
    real*8 :: d,p,t,q,s,b !!// 中间变量
    integer*4 :: is(size(cpx,1)),js(size(cpx,1)) !!// 中间变量

    n=size(cpx,1)
    forall(i=1:n,j=1:n)
        ar(i,j) = real(cpx(i,j));ai(i,j) = imag(cpx(i,j))
    end forall
    flag=1
    do k=1,n
        d=0.0
        do i=k,n
            do j=k,n
                p=ar(i,j)*ar(i,j)+ai(i,j)*ai(i,j)
                if(p.gt.d) then
                    d=p
                    is(k)=i
                    js(k)=j
                end if
            end do
        end do
        if(d+1.0.eq.1.0) then
            flag=0
            write(*,*) "flag=0,复矩阵奇异！" 
            return !!// 矩阵奇异，退出逆矩阵求解
        end if
        do j=1,n
            t=ar(k,j)
            ar(k,j)=ar(is(k),j)
            ar(is(k),j)=t
            t=ai(k,j)
            ai(k,j)=ai(is(k),j)
            ai(is(k),j)=t
        end do
        do i=1,n
            t=ar(i,k)
            ar(i,k)=ar(i,js(k))
            ar(i,js(k))=t
            t=ai(i,k)
            ai(i,k)=ai(i,js(k))
            ai(i,js(k))=t
        end do
        ar(k,k)=ar(k,k)/d
        ai(k,k)=-ai(k,k)/d
        do j=1,n
            if(j.ne.k) then
                p=ar(k,j)*ar(k,k)
                q=ai(k,j)*ai(k,k)
                s=(ar(k,j)+ai(k,j))*(ar(k,k)+ai(k,k))
                ar(k,j)=p-q
                ai(k,j)=s-p-q
            end if
        end do
        do i=1,n
            if(i.ne.k) then
                do j=1,n
                    if (j.ne.k) then
                        p=ar(k,j)*ar(i,k)
                        q=ai(k,j)*ai(i,k)
                        s=(ar(k,j)+ai(k,j))*(ar(i,k)+ai(i,k))
                        t=p-q
                        b=s-p-q
                        ar(i,j)=ar(i,j)-t
                        ai(i,j)=ai(i,j)-b
                    end if
                end do
            end if
        end do
        do i=1,n
            if (i.ne.k) then
                p=ar(i,k)*ar(k,k)
                q=ai(i,k)*ai(k,k)
                s=(ar(i,k)+ai(i,k))*(ar(k,k)+ai(k,k))
                ar(i,k)=q-p
                ai(i,k)=p+q-s
            end if
        end do
    end do
    do k=n,1,-1
        do j=1,n
            t=ar(k,j)
            ar(k,j)=ar(js(k),j)
            ar(js(k),j)=t
            t=ai(k,j)
            ai(k,j)=ai(js(k),j)
            ai(js(k),j)=t
        end do
        do i=1,n
            t=ar(i,k)
            ar(i,k)=ar(i,is(k))
            ar(i,is(k))=t
            t=ai(i,k)
            ai(i,k)=ai(i,is(k))
            ai(i,is(k))=t
        end do
    end do
    forall(i=1:n,j=1:n)
        bcinv(i,j) = cmplx(ar(i,j),ai(i,j),8)
    end forall

end function

    function eye(a)
        implicit none
        integer*4 a, i, j
        real*8   eye(a,a)
        eye=0.0
        do j=1,a
            do i=1,a
                if(i==j)  eye(i,j)=1.0
            end do
        end do
    end function eye

end module