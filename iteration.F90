MODULE iteration

IMPLICIT none
CONTAINS

SUBROUTINE iterate_k_omega(g, alpha, beta, e, es, omega, smearing, prec, max_step, rank)
    INTEGER, PARAMETER :: dp = kind(1.0d0)
    REAL(kind=dp), INTENT(in) :: omega, smearing, prec
    INTEGER, INTENT(in) :: max_step,  rank
    INTEGER :: step, cnt
    COMPLEX(kind=dp) :: alpha(rank,rank), beta(rank,rank), omega_mat(rank, rank),&
     green(rank, rank), e(rank,rank), es(rank,rank)
    COMPLEX(kind=dp), INTENT(out) :: g(rank, rank)
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
    g = zInverse(rank, omega_mat - es)
END SUBROUTINE iterate_k_omega

SUBROUTINE iterate_k(g_list, alpha, beta, e, es, omega_list, smearing, prec, max_step, rank, omega_num)
    INTEGER, PARAMETER :: dp = kind(1.0d0)
    REAL(kind=dp), INTENT(in) :: smearing, prec, omega_list(omega_num)
    REAL(kind=dp) :: omega
    INTEGER :: max_step, rank, omega_num, cnt
    COMPLEX(kind=dp), INTENT(in) :: alpha(rank,rank), beta(rank,rank), e(rank,rank), es(rank,rank)
    COMPLEX(kind=dp) :: alpha0(rank,rank), beta0(rank,rank), e0(rank,rank), es0(rank,rank)
    COMPLEX(kind=dp), INTENT(out) :: g_list(rank, rank, omega_num)
    DO cnt = 1, omega_num
        omega = omega_list(cnt)
        alpha0 = alpha
        beta0 = beta
        e0 = e
        es0 = es
        call iterate_k_omega(g_list(:, :, cnt), alpha0, beta0, e0, es0, omega, smearing, prec, max_step, rank)
    ENDDO
END SUBROUTINE iterate_k

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