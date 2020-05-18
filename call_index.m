function ci = call_index(L); 

% Compute call indices for blue and fin whale calls from 1-minute LTSA.
%   Blue whale B call, 3rd harmonic
%   Fin whale 20 Hz calls
%
% Input: L, a 1-minute resolution LTSA {time, frequency, spectrum level}
%
% Algorithm: For both blue whale B calls and fin whale 20 Hz calls, the
% algorithm is the calibrated signal:noise ratio.
% * signal from the mean of adjacent 1 Hz peak bands for call
% * noise from the mean of separate 1 Hz bands in the quiet background
%
% Output: structured array with index values

% algorithm specification: methods and frequencies
cif.blue = [37 43 44 50]; cif.fin = [12 20 21 34]; 

% blue
[q,ia,ib] = intersect(cif.blue,L.freq);    
ci.blue = mean(L.ltsa(ib([2 3]),:)) ./ mean(L.ltsa(ib([1 4]),:));

% fin
[q,ia,ib] = intersect(cif.fin,L.freq);    
ci.fin = mean(L.ltsa(ib([2 3]),:)) ./ mean(L.ltsa(ib([1 4]),:));

end
