function funcdat = my_spm_batch_filt_funcdata(funcdat, ltrh, htrh, tr)

Fs = 1/tr; %Sample frequency
Ny = Fs/2; %Nyquist frequency

if ltrh==0 || ltrh<0.0, ltrh=-1; end
if htrh==Ny || htrh>Ny, htrh=-1; end

s = size(funcdat);

mean_funcdat = mean(funcdat,2);
funcdat = funcdat - repmat(mean_funcdat,[1,s(end)]);

funcdat = brant_Filter_FFT_Butterworth(funcdat, ltrh, htrh, Fs, 0);

funcdat = funcdat + repmat(mean_funcdat,[1,s(end)]);

clear mean_funcdat