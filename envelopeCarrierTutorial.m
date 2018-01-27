% envelopeCarrierTutorial
%
% Little program to help understand envelope/carrier effects
%
% 11/29/13  dhb  Wrote it.

%% Clear
clear; close all;

%% Parameters
carrierFs = [5 10 20 30 40];
envelopeF = 1;
nSeconds = 5;
samplesSec = 5000;

% Envelope case
%  - us
%  - stockman
envelopeCase = 'stockman';

% Nonlinerity case
%  - 'log'
%  - 'sqrt'
%  - 'fourthroot';
%  - 'halfrect'
%  - 'fullrect'
%  - 'sqrtfullrect'
%  - 'sqrtmod'
nonlinCase = 'sqrtfullrect';

% Early filter case
%  - 'gaussianconv';
earlyfilterCase = 'gaussianconv';
kernalSd = 0.01;

% String defs
outputDir = sprintf('%s_%s_%g_%g_',envelopeCase,nonlinCase,envelopeF,round(1000*kernalSd));
if (~exist(outputDir,'file'))
    mkdir(outputDir);
end
mtfOutFileRoot = outputDir;
mtfTitleStr = sprintf('Version: %s, nonlin: %s, envelope %g Hz, kernalSd = %g seconds',envelopeCase,nonlinCase,envelopeF,kernalSd);

%% Loop over carrier frequencies
for c = 1:length(carrierFs)
    close all;
    
    carrierF = carrierFs(c);
    carrierOutFileRoot = sprintf('%s_%s_%g_%g_%g_',envelopeCase,nonlinCase,carrierF,envelopeF,round(1000*kernalSd));
    carrierTitleStr = sprintf('Version: %s, nonlin: %s, carrier %g Hz, envelope %g Hz, kernalSd = %g seconds\n',envelopeCase,nonlinCase,carrierF,envelopeF,kernalSd);
    
    %% Generate time and frequency bases
    %
    % Use odd number of samples so that DC comes out
    % where I expect it after fftshift.
    nSamples = samplesSec*nSeconds+1;
    theTime = linspace(0,nSeconds,nSamples);
    theFreqStep = 1/nSeconds;
    maxFreq = (nSamples/2)*theFreqStep;
    theFrequencies = linspace(-maxFreq,maxFreq,nSamples);
    
    %% Generate signal
    switch (envelopeCase)
        case 'us'
            theEnvelope = sin(2*pi*envelopeF*theTime);
        case 'stockman'
            theEnvelope = (0.5 + 0.5*cos(2*pi*envelopeF*theTime));
        otherwise
            error('Leave the jungle, leave the amazonis.  You will not find the signal here.')
    end
    theCarrier = sin(2*pi*carrierF*theTime);
    theSignal = 1 + theEnvelope.*theCarrier;
    
    %% Apply early filter to signal
    switch (earlyfilterCase)
        case 'gaussianconv'
            % Compute a Gaussian convolution kernal, normalized to unity.
            % Compute with an odd number of samples, for convenience in interpretting FFT freqs.
            theKernalRangeSecs = 1;
            nKernalSamples = round(2*theKernalRangeSecs*samplesSec);
            if (round(nKernalSamples/2) == nKernalSamples/2)
                nKernalSamples = nKernalSamples+1;
            end
            theKernalTime = linspace(-theKernalRangeSecs,theKernalRangeSecs,nKernalSamples);
            theKernal = exp(-pi*((theKernalTime/kernalSd).^2));
            theKernal = theKernal/sum(theKernal(:));
            
            % Convolve with signal
            theFilteredSignal = conv(theSignal,theKernal,'same');
            
            % Get predicted MTF at carrier freq
            theKernalFFT = fftshift(fft(theKernal));
            theKernalFreqStep = 1/(2*theKernalRangeSecs);
            maxKernalFreq = (nKernalSamples/2)*theKernalFreqStep;
            theKernalFrequencies = linspace(-maxKernalFreq,maxKernalFreq,nKernalSamples);
            [nil,theKernalFreqIndexAtCarrier] = min(abs(theKernalFrequencies-carrierF));
            theKernalMTFFreqs(c) = theKernalFrequencies(theKernalFreqIndexAtCarrier(1));
            theKernalMTF(c) = abs(theKernalFFT(theKernalFreqIndexAtCarrier(1)));
            theKernalMTFAnalytic(c) = exp(-pi*((kernalSd*carrierF)^2));
        otherwise
            error('Leave the jungle, leave the amazonis.  You will not find the filter here.')
    end
    
    
    %% Nonlinearity
    switch (nonlinCase)
        case 'log'
            theNLSignal = log10(theFilteredSignal);
        case 'sqrt'
            theNLSignal = sqrt(theFilteredSignal);
        case 'fourthroot'
            theNLSignal = sqrt(sqrt(theFilteredSignal));
        case 'halfrect'
            tempSignal = theFilteredSignal-mean(theFilteredSignal);
            tempSignal(tempSignal < 0) = 0;
            theNLSignal = tempSignal + mean(theFilteredSignal);
        case 'fullrect';
            tempSignal = theFilteredSignal-mean(theFilteredSignal);
            tempSignal = abs(tempSignal);
            theNLSignal = tempSignal + mean(theFilteredSignal);
        case 'sqrtfullrect';
            tempSignal = theFilteredSignal-mean(theFilteredSignal);
            tempSignal = sqrt(abs(tempSignal));
            theNLSignal = tempSignal + mean(theFilteredSignal);
        case 'sqrtmod'
            tempSignal = theFilteredSignal-mean(theFilteredSignal);
            minTemp = min(tempSignal)
            tempSignal = sqrt(tempSignal-minTemp)+minTemp;
            theNLSignal = tempSignal + mean(theFilteredSignal);
        otherwise
            error('Leave the jungle, leave the amazonis.  You will not find the nonlinearity here.')
    end
    
    
    %% Plot the time signal   
    theTimeFig = figure; clf;
    set(gcf,'Position',[100 100 1500 1500]);
    subplot(5,1,1);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theTime,theCarrier,'r');
    xlim([0 nSeconds]); ylim([-1 1]);
    xlabel('Time (seconds)')
    ylabel('Carrier');
    title(carrierTitleStr);
    subplot(5,1,2);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theTime,theEnvelope,'r');
    xlim([0 nSeconds]); ylim([-1 1]);
    xlabel('Time (seconds)')
    ylabel('Envelope');
    subplot(5,1,3);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theTime,theSignal,'r');
    xlim([0 nSeconds]); ylim([0 2]);
    xlabel('Time (seconds)')
    ylabel('Signal');
    subplot(5,1,4);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theTime,theFilteredSignal,'r');
    xlim([0 nSeconds]); ylim([0 2]);
    xlabel('Time (seconds)')
    ylabel('Filtered Signal');
    subplot(5,1,5);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theTime,theNLSignal,'r');
    xlim([0 nSeconds]); ylim([0 2]);
    xlabel('Time (seconds)')
    ylabel('Nonlinear Signal');
    savefig(fullfile(outputDir,[carrierOutFileRoot 'Time']),theTimeFig,'png');
    
    %% Frequency domain
    theCarrierFFT = fftshift(fft(theCarrier));
    theEnvelopeFFT = fftshift(fft(theEnvelope));
    theSignalFFT = fftshift(fft(theSignal));
    theFilteredSignalFFT = fftshift(fft(theFilteredSignal));
    theNLSignalFFT = fftshift(fft(theNLSignal));
    plotMaxFreq = 16;
    
    %% Plot the frequency domain signal
    theFreqFig = figure; clf;
    set(gcf,'Position',[100 100 1500 1200]);
    subplot(5,1,1);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theFrequencies,abs(theCarrierFFT),'r');
    set(gca,'XTick',-plotMaxFreq:2:plotMaxFreq);
    xlim([-plotMaxFreq plotMaxFreq]);
    xlabel('Frequency (cycles/sec)')
    ylabel('Carrier');
    title(carrierTitleStr);
    subplot(5,1,2);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theFrequencies,abs(theEnvelopeFFT),'r');
    set(gca,'XTick',-plotMaxFreq:2:plotMaxFreq);
    xlim([-plotMaxFreq plotMaxFreq]);xlabel('Time (seconds)')
    ylabel('Envelope');
    subplot(5,1,3);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theFrequencies,log10(abs(theSignalFFT)),'r');
    set(gca,'XTick',-plotMaxFreq:2:plotMaxFreq);
    xlim([-plotMaxFreq plotMaxFreq]);
    xlabel('Time (seconds)')
    ylabel('Log10 Signal');
    subplot(5,1,4);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theFrequencies,log10(abs(theFilteredSignalFFT)),'r');
    set(gca,'XTick',-plotMaxFreq:2:plotMaxFreq);
    xlim([-plotMaxFreq plotMaxFreq]);
    xlabel('Time (seconds)')
    ylabel('Log10 Fltered Signal');
    subplot(5,1,5);
    set(gca,'FontName','Helvetica','FontSize',18);
    plot(theFrequencies,log10(abs(theNLSignalFFT)),'r');
    set(gca,'XTick',-plotMaxFreq:2:plotMaxFreq);
    xlim([-plotMaxFreq plotMaxFreq]);
    xlabel('Time (seconds)')
    ylabel('Log10 Nonlinear Signal');
    savefig(fullfile(outputDir,[carrierOutFileRoot 'Freq']),theFreqFig,'png');
    
    %% Extract power at envelope and twice envelope
    [nil,theZeroFreqIndex] = min(abs(theFrequencies));
    dcAmps(c) = abs(theNLSignalFFT(theZeroFreqIndex));
    [nil,theEnvFreqIndex] = min(abs(theFrequencies-envelopeF));
    envelopeAmpFreqs(c) = theFrequencies(theEnvFreqIndex(1));
    envelopeAmps(c) = abs(theNLSignalFFT(theEnvFreqIndex(1)));
    [nil,theEnv2FreqIndex] = min(abs(theFrequencies-2*envelopeF));
    envelope2AmpFreqs(c) = theFrequencies(theEnv2FreqIndex(1));
    envelope2Amps(c) = abs(theNLSignalFFT(theEnv2FreqIndex(1)));
    fprintf('Carrier %g Hz, envelope %g Hz\n',carrierF,envelopeF);
    fprintf('\tKernal mtf freq %0.1f\n',theKernalMTFFreqs(c));
    fprintf('\tEnvelope amp freq %0.1f, 2*envelope amp freq %0.1f\n',envelopeAmpFreqs(c),envelope2AmpFreqs(c));
    fprintf('\tDc amp, %0.3g, envelope amp %0.3g, 2*envelope amp %0.3g\n',dcAmps(c),envelopeAmps(c),envelope2Amps(c));
    fprintf('\tKernal amp %0.3g, analytic %0.3g\n',theKernalMTF(c),theKernalMTFAnalytic(c));
end

%% Rolloff figure
mtfFig = figure; clf; hold on
set(gcf,'Position',[300 300 800 800]);
set(gca,'FontName','Helvetica','FontSize',18);

% Normalize extracted amplitude MTF for max for each extracted freq
plotEnvelopeAmps = envelopeAmps/max([envelopeAmps(:) ; envelope2Amps(:)]);
plotEnvelope2Amps = envelope2Amps/max([envelopeAmps(:) ; envelope2Amps(:)]);

% Regress predicted MTF on what we're goint to plot, to match overall scale.
kernalEnvAmps = theKernalMTF(:)*(theKernalMTF(:)\plotEnvelopeAmps(:));
kernalEnv2Amps = theKernalMTF(:)*(theKernalMTF(:)\plotEnvelope2Amps(:));
kernalEnvAmpsAnalytic = theKernalMTFAnalytic(:)*(theKernalMTFAnalytic(:)\plotEnvelopeAmps(:));
kernalEnv2AmpsAnalytic = theKernalMTFAnalytic(:)*(theKernalMTFAnalytic(:)\plotEnvelope2Amps(:));
plot(carrierFs,plotEnvelopeAmps,'ro','MarkerSize',8,'MarkerFaceColor','r');
plot(carrierFs,plotEnvelope2Amps,'go','MarkerSize',8,'MarkerFaceColor','g');
plot(carrierFs,kernalEnvAmps,'k');
plot(carrierFs,kernalEnv2Amps,'k');
plot(carrierFs,kernalEnvAmpsAnalytic,'r:','LineWidth',2);
plot(carrierFs,kernalEnv2AmpsAnalytic,'g:','LineWidth',2);
xlim([0 max(carrierFs)]); ylim([0 1]);
xlabel('Carrier Frequency (Hz)');
ylabel('Normalized Recovered Early Filter');
title(mtfTitleStr);
legend({'Envelope','Twice Envelope'},4);
savefig(fullfile(outputDir,mtfOutFileRoot),mtfFig,'png');


