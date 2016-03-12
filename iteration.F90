MODULE iteration

IMPLICIT none
CONTAINS

SUBROUTINE iterate(eout, esout, alpha, beta, e, es, omega, smearing, prec, max_step, rank)
    INTEGER, PARAMETER :: dp = kind(1.0d0)
    REAL(kind=dp) :: omega, smearing, prec
    INTEGER :: max_step, step, rank, cnt
    COMPLEX(kind=dp) :: alpha(rank,rank), beta(rank,rank), omega_mat(rank, rank),&
     green(rank, rank), e(rank,rank), es(rank,rank)
    COMPLEX(kind=dp), INTENT(out) :: eout(rank, rank), esout(rank, rank)
    step = 0
    omega_mat = 0
    DO cnt=1, rank
        omega_mat(cnt, cnt) = 1.0_dp
    ENDDO
    omega_mat = omega_mat * (omega + (0.0_dp, 1.0_dp) * smearing)
    DO WHILE (maxval(abs(alpha)) .gt. prec .or.  maxval(abs(beta)) .gt. prec .and. step .lt. max_step)
        green = zInverse(rank, omega_mat - e)
        e = e + matmul(alpha, matmul(green, beta)) + matmul(beta, matmul(green, alpha))
        es = es + matmul(alpha, matmul(green, beta))
        alpha = matmul(alpha, matmul(green, alpha))
        beta = matmul(beta, matmul(green, beta))
        step = step + 1
    ENDDO
    eout = e
    esout = es
END SUBROUTINE iterate

FUNCTION zInverse(n, a)  RESULT(ra)
    INTEGER, PARAMETER :: dp = kind(1.0d0)
    INTEGER :: n,lda,ipiv(n),info,lwork
    COMPLEX(kind=dp) :: a(n,n),ra(n,n),work(n)
    ra=a
    lwork=n
    lda=n
    CALL zgetrf(n, n, ra, lda, ipiv, info)
    IF(info/=0) WRITE(0,*) 'Error occured in zgetrf!'
    CALL zgetri(n, ra, lda, ipiv, work, lwork, info)
    IF(info/=0) WRITE(0,*) 'Error occured in zgetri!'
END FUNCTION zInverse
END MODULE iteration
