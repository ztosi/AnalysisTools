function [ wvlt ] = morelet( freq, tempPres )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    T = -t:tempPres:1;
    sine_wv = real(exp(2*pi*1i*freq.*T));
    sig = 6/(2*pi*freq);
    gauss_win = exp(-T.^2./(2*s^2));
    wvlt = sine_wv * gauss_win;
end

