Subroutine Secant(root,Xold,epsi1,epsi2,maxit,numit,f,fp,flag)
!------------------------------------------------------------------------------
! Ryan Kaplan
!
!Generic Subroutine that uses newtons method of root finding for a given 
!function and its derivative.
!
!Identifiers:
!     root:     Final root estimate for f(x) (intent out)        
!     Xold:     Initial root estimate (intent in/out)  
!     epsi1:    A convergence tolerance used if root is found(intent in)
!     epsi2:    2nd tolerance for slow progress(intent in)  
!     maxit:    The max number of iterations allowed (intent in)
!     numit:    Counts the number of iterations (intent out)
!     f:        Function to find root
!     fp:       Functions derivative
!     flag      1-4 flag used to let user know if estimate is accurate(out) 
!     fxnew:    Function evaluated at (xn+1) or Xnew (local)
!     Xnew:     Estimates the root through iterations (local)
!     Xolder:   Previous x value. (local)
!------------------------------------------------------------------------------
use types
Implicit None

real(dp),intent(inout):: Xold
real(dp),intent(out):: root
integer,intent(in):: maxit  
integer,intent(out):: numit, flag
real(dp),intent(in):: epsi1,epsi2

!Local variables used in subroutine
real(dp):: fxnew,fpXnew, Xnew, Xolder
integer:: n

!Interface for the external functions used to evaluate the problem at hand.
Interface
  Function f(x)
      use types
      real(dp),intent(in):: x
      real(dp):: f
      real(dp),parameter:: R= 0.082054_dp
  end function f
  Function fp(x)
      use types
      real(dp),intent(in):: x
      real(dp):: fp
  end function
End interface

!Initialize the number of iterations, fxnew,xnew,xolder,xold:
n=0 
!Xnew = Xold- f(Xold)/fp(Xold)
Xnew = Xold - (f(Xold)*(Xolder-Xold))/(f(Xolder)-f(Old))
fxnew= f(Xnew)
Xolder= Xold
Xold= Xnew

do
  !First criteria for found root (good root).
  if(abs(fxnew).LT.epsi1)then
    root= Xnew
    flag= 1
    numit=n
    exit
  end if
  !Second criteria for root (not a good approx)
  if(abs(xold-xolder).LT.epsi2)then
    root= Xnew
    flag=2
    numit=n
    exit
  end if
  !Third exit criteria divergence....
  if(n==maxit)then
    root= Xnew
    flag=3
    numit=n
    exit
  end if
  !No exit continue the itterations
  n= n+ 1
  fpXnew= fp(Xold)
  if(fpXnew==0)then
    stop "The derivative is zero!"
  end if
  Xnew= Xold - fXnew/fpXnew
  fXnew= f(Xnew)
  Xolder= Xold
  Xold= Xnew
  fxnew= f(xnew)
end do
return
End Subroutine Secant
