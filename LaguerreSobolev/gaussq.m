function[d,z,a1,b1,mu1]=gaussq(ic,n,alpha,beta,kpts,endpts)
% Compute Gauss rule corresponding to orthogonal polynomials
%INPUT
%       ic = 1 (Legendre), 2 (Chebyshev type I), 3 (Chebyshev type II), 
%            4 (Hermite), 5 (Jacobi), 6 (Laguerre)            
%       n = number of nodes in quadrature rule
%       alpha, beta = parameters in orthogonal polynomials
%       kpts = 
%       endpts = 
%OUTPUT
%       d = 
%       z = 
%       a1 = 
%       b1 = 
%       mu1 = 
%
[a1,b1,mu1]=class1(ic,n,alpha,beta);
if kpts == 1,
   [gam]=solve(endpts(1),n,a1,b1);
   a1(n)=gam*b1(n-1)^2+endpts(1);
elseif kpts==2,
   [gam]=solve(endpts(1),n,a1,b1);
   [gam1]=solve(endpts(2),n,a1,b1);
   t1=((endpts(1)-endpts(2))/(gam1-gam));
   b1(n-1)=sqrt(t1);
   a1(n)=gam*t1+endpts(1);
end


[d,e,z,ier]=gausq2(n,a1,b1);
z=z.^2*mu1;



 function[d,e,z,ierr]= gausq2(n, d, e)
%
%     this subroutine is a translation of an algol procedure,
%     num. math. 12, 377-383(1968) by martin and wilkinson,
%     as modified in num. math. 15, 450(1970) by dubrulle.
%     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
%     this is a modified version of the 'eispack' routine imtql2.
%
%     this subroutine finds the eigenvalues and first components of the
%     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
%     method.
%
%     on input:
%
%        n is the order of the matrix;
%
%        d contains the diagonal elements of the input matrix;
%
%        e contains the subdiagonal elements of the input matrix
%          in its first n-1 positions.  e(n) is arbitrary;
%
%        z contains the first row of the identity matrix.
%
%      on output:
%
%        d contains the eigenvalues in ascending order.  if an
%          error exit is made, the eigenvalues are correct but
%          unordered for indices 1, 2, ..., ierr-1;
%
%        e has been destroyed;
%
%        z contains the first components of the orthonormal eigenvectors
%          of the symmetric tridiagonal matrix.  if an error exit is
%          made, z contains the eigenvectors associated with the stored
%          eigenvalues;
%
%        ierr is set to
%          zero       for normal return,
%          j          if the j-th eigenvalue has not been
%                     determined after 30 iterations.
%
%     ------------------------------------------------------------------
%
      z=zeros(1,n); z(1)=1;



%
%
      ierr = 0;
      if  n ~= 1 ,
%
      e(n) = 0;
      for l = 1 : n,
         j = 0;
%     :::::::::: look for small sub-diagonal element ::::::::::
      ivar =1;
      while ivar ==1,
        for  m = l : n,
            if m == n,
               break
            end
            if abs(e(m)) <= eps * (abs(d(m)) + abs(d(m+1))),
               break
            end
 %          end  % if  m ~= n
        end % for m
%
       p = d(l);
         if  m ~= l,

         if j == 30,
            ierr=l;
            return,
         end
         j = j + 1;
%     :::::::::: form shift ::::::::::
         g = (d(l+1) - p) / (2 * e(l));
         r = sqrt(g*g+1);
         if g  < 0,
            ggg=-abs(r);
        else
            ggg=abs(r);
        end
         g = d(m) - p + e(l) / (g + ggg);
         s = 1;
         c = 1;
         p = 0;
         mml = m - l;

%     :::::::::: for i=m-1 step -1 until l do -- ::::::::::


          for ii = 1: mml, %      200 ii = 1, mml
            i = m - ii;
            f = s * e(i);
            b = c * e(i);
            if abs(f) >= abs(g),
              c = g / f;
              r = sqrt(c*c+1.0);
              e(i+1) = f * r;
              s = 1.0 / r;
              c = c * s;
            else
              s = f / g;
              r = sqrt(s*s+1.0);
              e(i+1) = g * r;
              c = 1.0 / r;
              s = s * c;
            end % if
            g = d(i+1) - p;
              r = (d(i) - g) * s + 2.0 * c * b;
            p = s * r;
            d(i+1) = g + p;
            g = c * r - b;
%     :::::::::: form first component of vector ::::::::::
            f = z(i+1);
            z(i+1) = s * z(i) + c * f;
            z(i) = c * z(i) - s * f;
        end % for
%
         d(l) = d(l) - p;
         e(l) = g;
         e(m) = 0.0;
        else
            ivar =2;

        end % if ivar
    %   end % forse questo non ci vuuole
        end % while
      end % for        240 continue
%
%     :::::::::: order eigenvalues and eigenvectors ::::::::::
      for  ii = 2 : n,
         i = ii - 1;
         k = i;
         p = d(i);
%
         for  j = ii : n,
            if d(j) <= p,
               k = j;
               p = d(j);
            end % if
        end % for   260    continue
%
         if k ~= i,
           d(k) = d(i);
           d(i) = p;
           p = z(i);
          z(i) = z(k);
          z(k) = p;
        end % if
      end % for
%
%     :::::::::: set error -- no convergence to an
%                eigenvalue after 30 iterations ::::::::::
 end % if      1001 return
%c     :::::::::: last card of gausq2 ::::::::::


function[a,b,muzero]=class1(kind, n, alpha, beta)
%
%           this procedure supplies the coefficients a(j), b(j) of the
%        recurrence relation
%
%             b p (x) = (x - a ) p   (x) - b   p   (x)
%              j j            j   j-1       j-1 j-2
%
%        for the various classical (normalized) orthogonal polynomials,
%        and the zero-th moment
%
%             muzero = integral w(x) dx
%
%        of the given polynomial's weight function w(x).  since the
%        polynomials are orthonormalized, the tridiagonal matrix is
%        guaranteed to be symmetric.
%
%           the input parameter alpha is used only for laguerre and
%        jacobi polynomials, and the parameter beta is used only for
%        jacobi polynomials.  the laguerre and jacobi polynomials
%        require the gamma function.
%

      nm1 = n - 1 ;

      if kind ==1,
%
%              kind = 1:  legendre polynomials p(x)
%              on (-1, +1), w(x) = 1.
%
      muzero = 2.0;
      for   i = 1: nm1,
         a(i) = 0.0;
         abi = i;
         b(i) = abi/sqrt(4*abi*abi - 1.0);
      end % for
      a(n) = 0.0;

      elseif kind==2,
%
%              kind = 2:  chebyshev polynomials of the first kind t(x)
%              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
%
      muzero = pi;
      for   i = 1: nm1,
         a(i) = 0.0;
         b(i) = 0.5;
      end % for
      b(1) = sqrt(0.5) ;
      a(n) = 0.0;

      elseif kind==3,
%
%              kind = 3:  chebyshev polynomials of the second kind u(x)
%              on (-1, +1), w(x) = sqrt(1 - x*x)
%
      muzero = pi/2.0;
      for  i = 1 : nm1,
         a(i) = 0.0d0 ;
         b(i) = 0.5d0;
      end % for i
      a(n) = 0.0d0;
%

      elseif kind==4,
%              kind = 4:  hermite polynomials h(x) on (-infinity,
%              +infinity), w(x) = exp(-x**2)
%
      muzero = sqrt(pi);
      for i = 1 : nm1,
         a(i) = 0.0;
         b(i) = sqrt(i/2);
      end % for i
      a(n) = 0.0;
%
      elseif kind==5,
%              kind = 5:  jacobi polynomials p(alpha, beta)(x) on
%              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
%              beta greater than -1
%
      ab = alpha + beta;
      abi = 2.0 + ab;
      muzero = 2.0^(ab+1.)*gamma(alpha+1)*gamma(beta+1) / gamma(abi);
      a(1) = (beta - alpha)/abi;
      b(1) = sqrt(4.0*(1.0 + alpha)*(1.0 + beta)/((abi + 1.0)*abi*abi));
      a2b2 = beta*beta - alpha*alpha;
      for  i = 2: nm1,
         abi = 2.0*i + ab;
         a(i) = a2b2/((abi - 2.0)*abi);
         b(i) = sqrt(4*i*(i+alpha)*(i+beta)*(i+ab)/((abi*abi-1)*abi*abi));
      end % for
      abi = 2.0*n + ab;
      a(n) = a2b2/((abi - 2.0)*abi);

        else
%
%              kind = 6:  laguerre polynomials l(alpha)(x) on
%              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
%              than -1.
%
      muzero = gamma(alpha + 1.0);
      for   i = 1 :nm1,
         a(i) = 2.0*i - 1.0 + alpha;
         b(i) =  sqrt(i*(i + alpha));
      end % for
      a(n) = 2.0*n - 1 + alpha;
    end % if


