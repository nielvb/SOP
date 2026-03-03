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


## Bibliography
[1] Van Buggenhout, N. "On generating Sobolev orthogonal polynomials" Numer. Math. 155, 415–443 (2023). https://doi.org/10.1007/s00211-023-01379-3
[2] Laudadio, T., Mastronardi, N. , Marcellán, F. , Van Buggenhout, N., and Van Dooren, P. "Recurrence relations and zeros of Gegenbauer-Sobolev orthogonal polynomials" Paper submitted (2026)
