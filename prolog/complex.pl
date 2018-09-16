/** <module> Complex number
 *
 * This is module to calculate imaginary number.
 *
 * The imaginary unit is represented by `i`.
 *
 * @author PiotrLi
 * @license GPL
 *
 */
:-module(complex,[
             iis/2,
             c_equals/2,
             complex_canonical/3,
             complex_exponential/3,
             complex_number/1
         ]).
:-op(700,xfx,user:iis).
:-op(700,xfx,iis).
:-op(700,xfx,user:c_equals).
:-op(700,xfx,c_equals).

%!  iis(-Number,++Expr) is det
%
%   It is is witch imaginary number.
%
%    * abs(+Z)
%    Return modulus of Z.
%
%    * phase(+Z)
%    Return phase of Z, normalized to be between -pi and pi.
%
%    * real(+Z)
%    Return real part of Z.
%
%    * imaginary(+Z)
%    Return imaginary part of Z.
%
%    * conjugate(+Z)
%    Return the conjugate of Z.
%
%    * reciprocal(+Z)
%    Return the reciprocal of Z.
%
%    * exp(+Z)
%    Return the complex exponential of Z.
%
%    * sin(+Z)
%    Return the complex sine of Z.
%
%    * cos(+Z)
%    Return the complex cosine of Z.
%
%    * tan(+Z)
%    Return the complex tangent of Z.
%
%   Examples:
%   ```prolog
%   ?- X iis i*i.
%   X = -1.
%
%   ?- X iis 5*i*6*7.
%   X = 210*i.
%
%   ?- X iis 5*i+3.
%   X = 3+5*i.
%
%   ?- X iis phase(2+3*i).
%   X = 0.982793723247329.
%
%   ?- X iis abs(5+10*i).
%   X = 11.180339887498949.
%
%   ?- X iis conjugate(1+5*i).
%   X = 1-5*i.
%
%   ?- X iis imaginary(5+2*i).
%   X = 2.
%
%   ?- X iis imaginary(sqrt(2)*exp(45*pi/180*i)).
%   X = 1.0.
%
%   ?- X iis 6*i*5/(i*2*5).
%   X = 3.
%
%   ```

%  Normal (is)/2
% ---------------
(_ iis Y):-
    \+ ground(Y),
    throw(error(instantiation_error, context(complex:(iis)/2, _))).

(X iis Y):-
    catch((Z is Y),_,(!,fail)),!,
    X = Z.

%  Functions
% -----------
(M iis abs(Z)):-
    ZZ iis Z,
    complex_exponential(ZZ,M,_),!.

(F iis phase(Z)):-
    ZZ iis Z,
    complex_exponential(ZZ,_,F),!.

(R iis real(Z)):-
    ZZ iis Z,
    complex_canonical(ZZ,R,_),!.

(I iis imaginary(Z)):-
    ZZ iis Z,
    complex_canonical(ZZ,_,I),!.

(R-I*i iis conjugate(Z)):-
    ZZ iis Z,
    complex_canonical(ZZ,R,I),!.

(X iis reciprocal(Z)):-
    ZZ iis Z,
    complex_canonical(ZZ,R,I),
    S is R*R+I*I,
    X iis R/S + I/S*i,!.

(X iis exp(Z)):-
    ZZ iis Z,
    complex_canonical(ZZ,R,I),
    X iis exp(R)*cos(I) + exp(R)*sin(I)*i,!.

(X iis sin(Z)):-
    ZZ iis Z,
    complex_canonical(ZZ,R,I),
    X iis sin(R)*cosh(I) + cos(R)*sinh(I)*i,!.

(X iis cos(Z)):-
    ZZ iis Z,
    complex_canonical(ZZ,R,I),
    X iis cos(R)*cosh(I) + sin(R)*sinh(I)*i,!.

(X iis tan(Z)):-
    ZZ iis Z,
    complex_number(ZZ),
    X iis sin(ZZ)/cos(ZZ),!.

(X iis ctan(Z)):-
    ZZ iis Z,
    complex_number(ZZ),
    X iis cos(ZZ)/sin(ZZ),!.

%  Optimize
% ----------
(-1 iis i*i):-!.

(0 iis 0*i):-!.

(0 iis 0.0*i):-!.

(i iis 1*i):-!.

(i iis 1.0*i):-!.

(0 iis O1+O2*i):-
    zero(O1),
    zero(O2),!.
(0 iis O1-O2*i):-
    zero(O1),
    zero(O2),!.

(R iis R+O*i):-
    number(R),
    zero(O),!.
(I*i iis O+I*i):-
    number(I),
    zero(O),!.
(II*i iis O-I*i):-
    number(I),
    zero(O),!,
    II is -I.

(R-II*i iis R+I*i):-
    number(R),
    number(I),
    I<0,!,
    II is -I.

(R+II*i iis R-I*i):-
    number(R),
    number(I),
    I<0,!,
    II is -I.

(X iis X):-
    complex_number(X),!.

%  Plus
% ------
(R iis R-I*i + I*i):-!.

(RR+I*i iis A+Z):-
    number(A),
    complex_number(Z),
    complex_canonical(Z,R,I),!,
    RR is A+R.

(XX iis A*i + B*i):-
    number(A),
    number(B),!,
    X is A+B,
    XX iis X*i.

(XX iis R-A*i + B*i):-
    number(R),
    number(A),
    number(B),!,
    X is B-A,
    XX iis R+X*i.

(XX iis R+A*i + B*i):-
    number(R),
    number(A),
    number(B),!,
    X is A+B,
    XX iis R+X*i.

(RR+II*i iis Z1+Z2):-
    complex_canonical(Z1,R1,I1),
    complex_canonical(Z2,R2,I2),!,
    RR is R1+R2,
    II is I1+I2.

(X iis A+B):-
    AA iis A,
    BB iis B,
    complex_number(AA),
    complex_number(BB),!,
    X iis AA+BB.

%  Minus
% -------
(0 iis I*i - I*i):-
    number(I),!.

(I*i iis R+I*i - R):-!.

(R iis R+I*i - I*i):-!.

(I*i iis I1*i - I2*i):-
    number(I1),
    number(I2),!,
    I is I1-I2.

(XX iis R-A*i - B*i):-
    number(R),
    number(A),
    number(B),!,
    X is A+B,
    XX iis R-X*i.

(XX iis R+A*i - B*i):-
    number(R),
    number(A),
    number(B),!,
    X is A-B,
    XX iis R+X*i.


(RR+I*i iis Z-A):-
    number(A),
    complex_canonical(Z,R,I),!,
    RR is R-A.

(RR+II*i iis Z1-Z2):-
    complex_canonical(Z1,R1,I1),
    complex_canonical(Z2,R2,I2),!,
    RR is R1-R2,
    II is I1-I2.

(X iis A-B):-
    AA iis A,
    BB iis B,
    complex_number(AA),
    complex_number(BB),!,
    X iis AA-BB.

%  Time
% ------
(A*i iis i*A):-!.

(AA iis A*i*i):-!,
    AA iis -A.

(X iis A*i*B*i):-!,
    X iis -A*B.

(I*i iis A*i*B):-
    number(A),
    number(B),!,
    I is A*B.

(RR+II*i iis A*Z):-
    number(A),
    complex_canonical(Z,R,I),!,
    RR is A*R,
    II is A*I.

(RR+II*i iis Z1*Z2):-
    complex_canonical(Z1,R1,I1),
    complex_canonical(Z2,R2,I2),!,
    RR is R1*R2 - I1*I2,
    II is R1*I2 + I1*R2.

(X iis A*B):-
    AA iis A,
    BB iis B,
    complex_number(AA),
    complex_number(BB),!,
    X iis AA*BB.

%  Divide
% -------
(1 iis i/i):-!.
(C*i iis i/A):-
    number(A),!,
    C is 1/A.
(C iis i/(A*i)):-
    number(A),!,
    C is 1/A.
(A iis A*i/i):-
    number(A),!.
(C iis A*i/(B*i)):-
    number(A),
    number(B),!,
    C is A/B.
(C*i iis A*i/B):-
    number(A),
    number(B),!,
    C is A/B.


(RR+II*i iis A/B):-
    number(B),
    AA iis A,
    complex_canonical(AA,R,I),!,
    RR is R/B,
    II is I/B.

(X iis A/B):-
    complex_number(A),
    complex_number(B),
    X iis A*reciprocal(B),!.

(X iis A/B):-
    AA iis A,
    BB iis B,
    complex_number(AA),
    complex_number(BB),!,
    X iis AA/BB.



%!  c_equals(@Term1,@Term2) is det
%
%   It is for (iis)/2, like (=:=)/2 for (is)/2.

c_equals(A,B):-
    Z1 iis A,
    complex_canonical(Z1,R1,I1),
    Z2 iis B,
    complex_canonical(Z2,R2,I2),
    R1==R2,
    I1==I2.

%!  complex_canonical(+Complex,-Real:number,-Imaginary:number) is semidet
%
%   Get abs and phase from complex number and test
%   if complex_number(Complex) in the same time.

complex_canonical(i,O,J):-
    zero(O),
    one(J),!.
complex_canonical(R+i,R,J):-
    number(R),
    one(J),!.
complex_canonical(R-i,R,-1):-
    number(R),!.
complex_canonical(I*i,0,I):-
    number(I),!.
complex_canonical(R+I*i,R,I):-
    number(R),
    number(I),!.
complex_canonical(R,R,0):-
    number(R),!.
complex_canonical(R-II*i,R,I):-
    number(R),
    number(II),
    I is -II,!.
complex_canonical(Z,R,I):-
    (   Z=M*exp(i*F)
    ;   Z=M*(cos(F)+i*sin(F))
    ),
    pol_rec(R,I,M,F),!.

complex_canonical(R+I*i,R,I):-
    number(R),
    number(I),!.

%!  complex_exponential(+Complex,-Abs:number,-Phase:number) is semidet
%
%   Get abs and phase from complex number and test
%   if complex_number(Complex) in the same time.

complex_exponential(i,O,J):-
    zero(O),
    one(J),!.
complex_exponential(M,M,O):-
    number(M),
    zero(O),!.
complex_exponential(M*i,M,F):-
    number(M),
    F is pi/2,!.
complex_exponential(R+i,M,F):-
    number(R),
    pol_rec(R,1,M,F),!.
complex_exponential(R-i,M,F):-
    number(R),
    pol_rec(R,-1,M,F),!.
complex_exponential(R+I*i,M,F):-
    number(R),
    number(I),
    pol_rec(R,I,M,F),!.
complex_exponential(R-I*i,M,F):-
    number(R),
    number(I),
    II is I,
    pol_rec(R,II,M,F),!.
complex_exponential(M*exp(i*F),M,F):-
    number(M),
    number(F),!.
complex_exponential(M*(cos(F)+i*sin(F)),M,F):-
    number(M),
    number(F),!.

pol_rec(R,I,M,F):-
    number(R),
    number(I),
    M is sqrt(R*R+I*I),
    F is atan2(I,R),!.
pol_rec(R,I,M,F):-
    number(M),
    number(F),
    R is M*cos(F),
    I is M*sin(F),!.

%!  complex_number(@Complex) is semidet
%
%   True if Term currently is a complex number in a form like:
%    * canonical
%    ```
%    4, i, 2*i, 3+2*i, 8-4*i, 3+ -2*i
%    ```
%
%    * exponential
%    ```
%    12*exp(i*5)
%    ```
%
%    * trigonometric
%    ```
%    10*(cos(3)+i*sin(3))
%    ```

complex_number(A):-
    var(A),!,fail.

complex_number(i):-!.
complex_number(R):-
    number(R),!.
complex_number(R + I*i):-
    number(R),
    number(I),!.
complex_number(R - I*i):-
    number(R),
    number(I),!.
complex_number(R - i):-
    number(R),!.
complex_number(M*exp(i*F)):-
    number(M),
    number(F),!.
complex_number(M*(cos(F)+i*sin(F))):-
    number(M),
    number(F),!.
complex_number(I):-
    imaginary_number(I,_),!.

imaginary_number(i,1):-!.
imaginary_number(I*i,I):-
    number(I),!.

zero(0):-!.
zero(0.0):-!.

one(1):-!.
one(1.0):-!.
