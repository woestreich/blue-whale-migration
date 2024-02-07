function L = calibratedLTSA(X,P)

% Compute calibrated LTSA.  called by batch_process_LTSA.m

% input
%   X = input audio, structured array {time, x, FS}
%   P = parameters of input audio and LTSA parameters
%     {tave, fs, dfreq, nfft, noverlap, sa, nspec, cal}
%
% output
%   L = LTSA output (time, frequency, power), units are dB re uPa^2/Hz

X.x = X.x * (6.0 / 2.0); % voltage scaling 

window = hann(P.nfft);  % Hanning window


for N = 1:P.nspec
    cs = (N-1)*P.sa+1; ce = cs+P.sa - 1; cdat = X.x(cs:ce);
    [pwr,freq] = pwelch(cdat,window,P.noverlap,P.nfft,P.fs);
    % power in volts
    pwr = 10*log10(pwr); 
    % correct for the Hann window
    pwr = pwr - 10*log10(1.5); 
    % power spectral density in volts, offset is zero if np = FS
    %   where FS = sampling frequency and np = ???
    psd = pwr - 10*log10(X.FS/P.nfft); 
    % add to LTSA
    ltsa(:,N) = pwr(:);
end

% Correct for hydrophone sensitivity
%xp = P.cal(:, 1); yp = P.cal(:, 2);
%s = interp1(xp,yp,freq,'spline');
%s = repmat(s,1,floor(P.nspec)); 
%upa = ltsa - s;
% currently not using hydrophone calibration
upa = ltsa;

L.time = [0:P.nspec-1]*P.tave; + 0.5*P.tave;  % seconds
L.freq = freq;
L.ltsa = upa;

