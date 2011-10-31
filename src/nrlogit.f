C ===========================================================================
      subroutine nrlogit(x1,beta,p,n)

      implicit integer (I-N) 
      double precision beta(3),p(n),fx,dfx,x1,logitp
      J=0
      fx = 1.0
      DO 100 I = 1, n
         logitp = LOG(p(I)/(1-p(I)))
 10      IF((J .LT. 100) .AND. (ABS(fx) .GT. 0.000001)) THEN
            fx = beta(1)-logitp+beta(2)*x1+beta(3)*LOG(x1)
            dfx = beta(2)+beta(3)/x1
            x1 = x1 - fx/dfx
            IF(x1<0) THEN x1 = 0.001
            J = J+1
            GOTO 10
         ENDIF

         p(I) = x1
 100  END DO
      END

