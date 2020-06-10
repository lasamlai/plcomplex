/** <module> Cooley–Tukey FFT algorithm
 *
 *  Compute the FFT and inverse FFT of a length n complex sequence
 *  using the radix 2 Cooley-Tukey algorithm.
 *
 *  Bare bones implementation that runs in O(n log n) time.
 *
 *  Limitations
 *  -----------
 *   -  assumes length is a power of 2
 *   -  not the most memory efficient algorithm (because it uses
 *      an object type for representing complex numbers and because
 *      it re-allocates memory for the subarray, instead of doing
 *      in-place or reusing a single temporary array)
 *
 *  @see https://introcs.cs.princeton.edu/java/97data/FFT.java.html
 *  @author PiotrLi
 *  @license GPL
 */
:-module(ct_fft,[
             ct_fft/2,
             ct_ifft/2,
             ct_cconvolve/3,
             ct_convolve/3
         ]).

:-use_module(complex).

%!  ct_fft(+X:list,-Y:list) is det
%
%   Compute the FFT of X, assuming its length is a power of 2.

ct_fft([],[]):-!.
ct_fft([X],[X]):-!.
ct_fft(X,Y):-
    length(X,N),
    length(Y,N),

    deal(X, Event, Odd),

    ct_fft(Event,Q),
    ct_fft(Odd,R),
    Th is -2 * pi / N,

    divide(Y, Y1, Y2),

    % combine
    combine(Q,R,Th,0,Y1,Y2).

combine([],[],_,_,_,_).
combine([Q|QT],[R|RT],Th,K,[Yk|Y1],[Ykn2|Y2]):-
    WkR iis R*exp(i*(K*Th)),
    Yk iis Q+WkR,
    Ykn2 iis Q-WkR,
    succ(K,KK),
    combine(QT,RT,Th,KK,Y1,Y2).

deal([], [], []) :- !.
deal([A], [A], []) :- !.
deal([A,F,B,G,C,H,D,I,E,J|T],[A,B,C,D,E|L],[F,G,H,I,J|R]):-
    !,
    deal(T,L,R).
deal([A,F,B,G,C,H,D,I|T],[A,B,C,D|L],[F,G,H,I|R]):-
    !,
    deal(T,L,R).
deal([A,F,B,G,C,H|T],[A,B,C|L], [F,G,H|R]):-
    !,
    deal(T,L,R).
deal([A,F,B,G|T],[A,B|L], [F,G|R]):-
    !,
    deal(T,L,R).
deal([A,F|T],[A|L],[F|R]):-
    !,
    deal(T,L,R).

divide(A,B,C):-
    divide_(A, B, A, C).

divide_([], [], C, C).
divide_([_,_|T],[X|L],[X|LL], C):-
    divide_(T, L, LL, C).

%!  ct_ifft(+X:list,-Y:list) is det
%
%   Compute the inverse FFT of X, assuming its length is a power of 2.

ct_ifft([],[]):-!.
ct_ifft(X,Y):-
    length(X,N),

    % take conjugate
    conjugate(X,Y1),

    % compute forward FFT
    ct_fft(Y1,Y2),

    % take conjugate again
    conjugate(Y2,Y3),

    % divide by n
    maplist({N}/[AA,BB]>>(BB iis AA/N),Y3,Y),!.

conjugate([],[]):-!.
conjugate([A|AT],[B|BT]):-
    B iis conjugate(A),
    conjugate(AT,BT),!.

%!  ct_cconvolve(+X:list,+Y:list,-Z:list) is det
%
%   Compute the circular convolution of X and Y.

ct_cconvolve(X,Y,Z):-

    % compute FFT of each sequence
    ct_fft(X,A),
    ct_fft(Y,B),

    % point-wise multiply
    maplist([AA,BB,CC]>>(CC iis AA*BB),A,B,C),

    % compute inverse FFT
    ct_ifft(C,Z).

%!  ct_convolve(+X:list,+Y:list,-Z:list) is det
%
%   Compute the linear convolution of Y and Y.

ct_convolve(X,Y,Z):-
    length(X,N),
    N2 is N*2,
    length(A,N2),
    length(B,N2),
    length(ZZ,N),
    maplist([0]>>true,ZZ),
    append(X,ZZ,A),
    append(Y,ZZ,B),
    ct_cconvolve(A,B,Z).
