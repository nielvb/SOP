# SOP
## General purpose code for sequentially dominated diagonal Sobolev inner products:
General purpose code to generate Sobolev orthonormal polynomials as described in [1]:
1) Arnoldi iteration 'Arnoldi.m'
2) Updating procedure 'updating.m'

The modified Chebyshev method and discretized Stieltjes procedure, together with code to generate quadrature rules can be found on the webpage of Walter Gautschi: https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html

The folder 'Examples/' contains several examples how to use the code to generates Sobolev orthogonal polynomials.
To run this code the following functions are required (available on https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html):
- 'Chebyshev_sob.m'
- 'Stieltjes_sob.m'
- 'gauss.m'
- 'r_jacobi.m'
- 'r_laguerre.m'
- 'sobzeros.m'

## Specialized code for Gegenbauer-Sobolev polynomials:
Folder 'GegenbauerSobolev/' contains specialized code to generate Gegenbauer-Sobolev orthogonal polynomials [2]
1) The main file 'main.m' illustrates the use of the code
2) To generate the recurrence relation for Gegenbauer-Sobolev polynomials use 'GegSob_generate.m'
3) To generate the quantities alpha (as described in [2]), use 'computerAlpha.m'


## Specialized code for Laguerre-Sobolev polynomials:
Folder 'LaguerreSobolev/' contains specialized code to compute the zeros of Laguerre-Sobolev orthogonal polynomials [3,4]
1) The two methods described in [3]  (the Ehrlich-Aberth method and the solution of a generalized eigenvalue problem) are implemented in the file 'Laguerre_Sobolev_zeros_NA25F.m'
2) An adapted version of the Ehrlich-Aberth method that avoids overflow [4] is given in 'Comput_zeros_LagSobF.m' and a multiprecision version using ADVANPIX is given in 'Comput_zeros_LagSob_mp.m'
3) 'Example2_ex_time.m' is a numerical experiment which compares the execution time of the three methods, and corresponds to Figure 1 in [4]


## Bibliography
[1] Van Buggenhout, N. "On generating Sobolev orthogonal polynomials" Numer Math 155, 415–443 (2023). https://doi.org/10.1007/s00211-023-01379-3 \
[2] Van Buggenhout, N., Laudadio, T., Mastronardi, N. et al. "Recurrence relations and zeros of Gegenbauer–Sobolev orthogonal polynomials". Adv Comput Math 52, 63 (2026). https://doi.org/10.1007/s10444-026-10338-z \
[3] Laudadio, T., Mastronardi, N., Marcellán Español, F.J. et al. "On computing the zeros of Laguerre–Sobolev polynomials". Numer Algor 100, 1507–1526 (2025). https://doi.org/10.1007/s11075-025-02021-z \
[4] Laudadio, T., Marcellán, F., Mastronardi, N. et al. "Computation of the zeros of Laguerre–Sobolev polynomials by the Ehrlich–Aberth method". (2026) Submitted
