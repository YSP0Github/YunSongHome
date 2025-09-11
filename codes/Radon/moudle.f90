!********************************************************************
!  my_module
!
!  PURPOSE:   复数之间和实数之间的叉乘计算
!********************************************************************
module my_module



    contains    
     
    !-----------递归求行列式的值--------------------
    recursive function det(A,col,row) result(D)
    Implicit None
    integer*4 row,col
    real*8 ::A(col,row),B(row-1,col-1)
    real*8 ::D
    integer*4 row_now,col_now,k,c,f
    row_now=row 
    col_now=col
    if (row_now>1) then
        D = 0.0;
        do k=1,row_now
            if(k==1)then
                B=A(2:col_now,2:row_now)
            elseif(k<row_now)then
                B(1:k-1,:)=A(1:k-1,2:row_now)
                B(k:col_now-1,:)=A(k+1:col_now,2:row_now)
            else
                B=A(1:col_now-1,2:row_now)
            endif
            c=col-1
            f=row-1
            D = D + (-1)**(1+k) * A(k,1) *&
                det(B,c,f)
        end do
    else
        D = A(1,1);
    end if
    end function det
 
 
    function delete_col_row(A,m,n,k,tar) result(b)
        !tar == 1，在 m*n 的矩阵中删去第 k 行
        !tar != 1，在 m*n 的矩阵中删去第 k 列,方便伴随矩阵的计算
        Implicit None
        integer*4::m,n,k,tar
        real*8::a(m,n)
        real*8,allocatable::B(:,:)
        if(tar==1)then
            allocate(B(m-1,n))
            if(k==1)then
                B=A(2:m,:)
            elseif(k<m)then
                B(1:k-1,:)=A(1:k-1,:)
                B(k:m-1,:)=A(k+1:m,:)
            else
                B=A(1:m-1,:)
            endif
        else
            allocate(B(m,n-1))
            if(k==1)then
                B=A(:,2:n)
            elseif(k<n)then
                B(:,1:k-1)=A(:,1:k-1)
                B(:,k:n-1)=A(:,k+1:n)
            else
                B=A(:,1:n-1)
            endif
        endif
    end function delete_col_row
 
! ------------用  A* / |A|  计算逆矩阵---------
    function inv(A,row) 
    Implicit None
    integer*4 ::row,i,j
    real*8 ::A(row,row)
    real*8 ::inv(row,row),b(row-1,row-1)
    do i=1,row
        do j=1,row
            b=delete_col_row(delete_col_row(A,row,row,i,1),row-1,row,j,2)
            inv(i,j)=(-1)**(i+j)*det(b,row-1,row-1)/det(A,row,row)
        end do
    end do
    inv=transpose(inv)
    end  function inv
    
    ! ******   实数向量叉乘函数   *********已验证
    function vec_product(v1, v2)
        implicit none
        real(8), dimension(3) :: vec_product
        real(8), dimension(3), intent(in) :: v1, v2

        vec_product(1) = v1(2) * v2(3) - v2(2) * v1(3)
        vec_product(2) = v1(3) * v2(1) - v2(3) * v1(1)
        vec_product(3) = v1(1) * v2(2) - v2(1) * v1(2)
    end function
    
    ! ******   两个复数向量叉乘函数   *********已验证
    function vec_productcc(v1, v2)
        implicit none
        Complex(16), dimension(3) :: vec_productcc
        Complex(16), dimension(3), intent(in) :: v1, v2

        vec_productcc(1) = v1(2) * v2(3) - v2(2) * v1(3)
        vec_productcc(2) = v1(3) * v2(1) - v2(3) * v1(1)
        vec_productcc(3) = v1(1) * v2(2) - v2(1) * v1(2)
    end function
    
    ! ******   复数  实数  向量叉乘函数   *********已验证
    function vec_productcr(v1, v2)
        implicit none
        Complex(16), dimension(3) :: vec_productcr
        Complex(16), dimension(3), intent(in) :: v1
        Real(8), dimension(3), intent(in) :: v2

        vec_productcr(1) = v1(2) * v2(3) - v2(2) * v1(3)
        vec_productcr(2) = v1(3) * v2(1) - v2(3) * v1(1)
        vec_productcr(3) = v1(1) * v2(2) - v2(1) * v1(2)
    end function
    ! ******   实数  复数  向量叉乘函数   *********已验证
    function vec_productrc(v1, v2)
        implicit none
        Real(8), dimension(3), intent(in) :: v1
        Complex(16), dimension(3), intent(in) :: v2
        Complex(16), dimension(3) :: vec_productrc

        vec_productrc(1) = v1(2) * v2(3) - v2(2) * v1(3)
        vec_productrc(2) = v1(3) * v2(1) - v2(3) * v1(1)
        vec_productrc(3) = v1(1) * v2(2) - v2(1) * v1(2)
    end function



end module


!program main
!
!    use operator_i
!    implicit none
!    real*8 :: r(2,2),ri(2,2)
!    integer*4 :: i,j
!    complex*16 :: c(4,4),ci(4,4)
!
!    c=reshape((/(0.2368,0.1345),(1.1161,1.2671),(0.1582,-0.2836),(0.1968,0.3576),&
!    (0.2471,0.1678),(0.1254,0.2017),(1.1675,-1.1967),(0.2071,-1.2345),&
!    (0.2568,0.1825),(0.1397,0.7024),(0.1768,0.3558),(1.2168,2.1185),&
!    (1.2671,1.1161),(0.1490,0.2721),(0.1871,-0.2078),(0.2271,0.4773)/),(/4,4/))
!    ci=.i.c
!    write(*,*) "复矩阵："
!    write(*,"(8(4x,es10.3))") ((c(i,j),j=1,4),i=1,4)
!    write(*,*) "逆矩阵："
!    write(*,"(8(4x,es10.3))") ((ci(i,j),j=1,4),i=1,4)
!    write(*,*) "校核："
!    write(*,"(8(4x,es10.3))") matmul(c,ci)
!
!    r=reshape((/1,3,2,4/),(/2,2/))
!    ri=.i.r
!    write(*,*) "实矩阵："
!    write(*,"(2(4x,es10.3))") ((r(i,j),j=1,2),i=1,2)
!    write(*,*) "逆矩阵："
!    write(*,"(2(4x,es10.3))") ((ri(i,j),j=1,2),i=1,2)
!    write(*,*) "校核："
!    write(*,"(2(4x,es10.3))") matmul(r,ri)
!
!end program