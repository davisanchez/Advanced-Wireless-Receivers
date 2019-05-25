function  [EbN0Range]  = EbN0shift(P)
%
% The function takes the parameter variable P as a input
% and outputs the shifted Eb/N0 range corresponding to the 
% given SNR range. 
%
 
G = length(P.Sequence);

shift = 10*log10(G);

EbN0Range = P.SNRRange + shift;