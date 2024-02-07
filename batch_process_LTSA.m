% Process time-series of wav files files into calibrated LTSA files.
%
% MBARI specific processing -
% Input is daily 1 kHz wav files.
% Output is daily calibrated power spectral density at 1 minute resolution.
%
% Because all decimated wav files are padded with zeros for missing data
% to maintain accurate timekeeping, LTSAs for the zero periods become -Inf.
% You could exclude these and associated time reference before saving the
% LTSA so that there is no wasted storage. Currently, you are filtering out
% the missing LTSA periods downstream using the -Inf as indicator.
%
% 30-Nov-2018   John Ryan, MBARI, ryjo@mbari.org
% 02-Jan-2023   modified by Will Oestreich, MBARI

clear all;
close all;

dpath = '/Users/woestreich/Dropbox/Documents/MBARI_postdoc/BWO_2022/PAM_data/300m/1kHz/';
flist = dir([dpath,'*.wav']); nfiles = numel(flist);
    
for F = 1:nfiles
    
    disp(['Processing ' flist(F).name]);
    disp(['File ', num2str(F), ' of ', num2str(nfiles)]);

    % Read
    fn = [flist(F).folder '/' flist(F).name]; disp(['Reading ' fn]);
    [X.x,X.FS] = audioread(fn);
            
    % File start time (for calibration and time reference)
    [filepath,filename,ext] = fileparts(fn);
    filestart = datenum(filename(6:17),'yymmddHHMMSS');
            
    % LTSA parameterization
    P.tave = 60;  % output LTSA averaging time bin [seconds] ****
    P.fs = X.FS;  % input sample rate [Hz]
    P.dfreq = 1;  % frequency bin size [Hz]
    P.nfft = P.fs/P.dfreq;  % number of ffts
    overlap = 50;  % percent
    P.noverlap = round((overlap/100)*P.nfft);
    P.sa = P.tave * P.fs;  % samples per LTSA spectral average
    P.nspec = numel(X.x)/P.sa;
            
    % Calibration
    % P.cal = get_icListen_cal(filestart);
            
    % LTSA
    L = calibratedLTSA(X,P);
                                    
    % subset to primary frequency range of interest
    % idx = find(L.freq >= flim(1) & L.freq <= flim(2));
    % L.freq = L.freq(idx); L.ltsa = L.ltsa(idx,:);
    L.time = filestart + L.time/86400;
                        
    % save
    ofile = ['/Users/woestreich/Dropbox/Documents/MBARI_postdoc/BWO_2022/PAM_data/300m/1kHzLTSA/' filename(6:17) '.mat'];
    disp(['Saving ' ofile]); disp(' '); eval(['save ' ofile ' L']);

end

