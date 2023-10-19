function [modf, peak_times, freq]= audmod (file, frame_rate, cutoff, full)
%
% calculate audio modulation function using MFCCs
%
% Louis Goldstein
% 1 March 2022
% Input arguments
% file:        audio file or cell structure containing the waveform in the first cell and the sampling rate in the second  
% frame_rate   frame rate of MFCC analysis (Hz, (def. = 200 Hz)
% cutoff       filter cutoff (Hz, def. = 12 Hz)
% Output arguments:
% modf          modulation function
% peak_times    time of peaks (secs)
% freq          mean pulse frequency
% full          if equal to one statistics are computed and results are plotted(added by L. Lancia)

% 11/10/2023 L. Lancia modified to accept cell structure as first argument and added the possibility to skip the computation of statistics and plots

close all;

if nargin<2; frame_rate = 200; end;
if nargin<3; cutoff=12; end;
if nargin<4; full=0; end;


% get audio
if iscell(file) % added by L. Lancia to process a waveform and its sampling rate as first argument
    s=file{1};
    sr=file{2};
else
    [s, sr] = audioread(file);
end
% calculate filter parameters

Wn = cutoff/(frame_rate/2);
N = 9;
[B,A] = butter(N,Wn);

 %%
 

% get mfccs and filter them
[mfccs FBEs] = get_mfcc(s, sr, 1000/frame_rate); 

[Ncoef, Nframes]= size(mfccs);

for i=1:size(mfccs,1)
    mfccs(i,:) = filtfilt(B,A,mfccs(i,:));
end

% calculate modulation function ignoring mfcc 1
for i = 2:size(mfccs,2)
    mdist(i) = sum((mfccs (2:Ncoef,i) - mfccs (2:Ncoef,i-1)).^2);
end  
 
modf = filtfilt(B,A,mdist);
%% 
if full ==1 % added by L. Lancia to skip statistics and plots
    % calculate peak locations as zero crossings of the derivative of modf
    dmodf = diff(modf);
    dmodf_shift = [0 dmodf(1:end-1)];
    peaks = find ((dmodf_shift>0)&(dmodf<0));
    peak_locs = ((dmodf_shift>0)&(dmodf<0));



    % get mean freq of peaks

    peak_times = peaks/frame_rate;
    ipi = peak_times(2:end) - peak_times(1:end-1);
    freq = 1/mean(ipi);

    %%

    % PLOTS

    t = [0 : Nframes-1].* (1000/frame_rate);
    figure (1)
    h=subplot (2,1,1);
    Nmsec = make_spect2 (s,sr, 5, 6000, 5);
    title (file);
    xlabel ('Time in milliseconds',  'fontsize', 14);
    hold on
    stem (t(1:end-1), peak_locs*6000);
    hold off

    subplot (2,1,2)
    plot (t, modf);
    xlim ([1 Nmsec]);
    ylabel ('MFCC modulatiom',  'fontsize', 14);
    hold on
    stem (t(1:end-1), peak_locs*max(modf));
    title ([' FREQ: ' num2str(freq)]);
    xlabel ('Time in milliseconds',  'fontsize', 14);
    hold off
else
    peak_times =[];
    freq=[];
end
