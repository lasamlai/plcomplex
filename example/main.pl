:-use_module(library(ct_fft)).

show(X,T):-
    write(T),nl,
    write('-------------------'),nl,
    show_(X).

show_([]).
show_([A|T]):-
    writeln(A),
    show_(T).

main(N):-
    length(X,N),

    % original data
    maplist([A]>>random(A),X),
    show(X,'X'),

    nl,
    % FFT of original data
    time(ct_fft(X,Y)),
    show(Y,'fft(Y,X)'),

    nl,
    % take inverse FFT
    time(ct_ifft(Y,Z)),
    show(Z,'ifft(Y,Z)'),

    nl,
    % circular convolution of x with itself
    time(ct_cconvolve(X,X,C)),
    show(C,'cconvolve(X, X, C)'),

    nl,
    % linear convolution of x with itself
    time(ct_convolve(X,X,CC)),
    show(CC,'convolve(X, X, C)').
