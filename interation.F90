SUBROUTINE INTERATE(alpha0, beta0, epsilon0, epsilons0, omega, smearing, prec, min_step, max_step, step, rank)
 REAL(kind=8) omega, smearing, prec
 INTEGER min_step, max_step, step, rank
 COMPLEX alpha0(rank), beta0(rank), epsilon0(rank), epsilons0(rank), omega_mat(rank), green(rank)
 omega_mat = 0
 DO INTEGER i=1,rank 
      omega_mat(i,i) = 1.0
 ENDDO 
 green = 
 

END 
