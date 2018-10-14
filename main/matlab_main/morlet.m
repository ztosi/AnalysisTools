function [ wvlt ] = morlet( freq, tempPres )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    sig = 6/(2*pi*freq);
    T = (-5*sig):tempPres:(5*sig);
    sine_wv = exp(2*pi*1i*freq.*T);
    gauss_win = exp(-T.^2./(2*sig^2));
    wvlt = sine_wv .* gauss_win;
end

