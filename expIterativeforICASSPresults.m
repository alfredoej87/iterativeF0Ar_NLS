addpath(genpath('\\CREATE.AAU.DK\Users\aeja\Desktop\WorkMatlab\fastF0Nls-master\fastF0Nls-master\matlab'));
%clear
load generalArctic32_256.mat 
load cb_noise_size16.mat
numberSpeechVectors = size(GeneralSpeechDictionaryOrder12,2); 
numberNoiseVectors_2 = size(noisex92Codebook2,2); 
FullDictionary_2 = [GeneralSpeechDictionaryOrder12 noisex92Codebook2];
snr1 = -5:10:15; 
nSnr = length(snr1); 
samplingFreqPitch = 100; 
Filenames = {'f1nw0000','m1nw0000','f2nw0000','m2nw0000','f3nw0000',...
   'm3nw0000','f4nw0000','m4nw0000', 'f5nw0000','m5nw0000'}; 
nFiles = length(Filenames); 
lpc_prewh_order = 25; 
p1 = lpc_prewh_order; 
noiseTypes = 2; 
fs = 8000; 
pitch_segment_duration_sec = 0.025625; 
N1 = pitch_segment_duration_sec*fs; 
shifting_time = 0.010; 
M1 = shifting_time*fs; 
Lmax = 27; 
f0BoundsHz = [60 400]; 
f0Bounds = f0BoundsHz/fs; 
f0Estimator = fastF0Nls(N1,Lmax,f0Bounds); 
%f0Estimator2 = fastF0Nls(N1,Lmax,f0Bounds); 
%%Parameters for SHRP method 
% frame_length_ms = N1/(fs/1000); % timestep = 10; % SHR_Thresh = 0.35; 
% ceiling = 1250; % med_smooth = 0;

%%pre-whitening parameters 
N2 = N1; M2=M1; 

for nF = 56:60%%52:60%%1:10%%nFiles
    filename = cell2mat(Filenames(nF-50));
    [speechKele,fs_Keele] = KeeleRead([filename '.pes']);
    s = resample(speechKele,2,5);  %%Resampled to 8kHz   
    lengthClean = length(s); 
    signalPower = s'*s;
    disp(['simulating for file ' num2str(filename)]);    
    %%Read peak values of correlation in the database    
    peakCorrValues = textread([filename '.pev'],'%s');
    if strcmp(filename,'m1nw0000')                  %extract pitch information
        peakCorrValues = str2double(peakCorrValues(41:end-2));
    else
        peakCorrValues = str2double(peakCorrValues(41:end-1));
    end
    total_pitch_frames = floor(lengthClean/M2)-2; 
    %%For GER only consider positive values, ignore 0 ones 
    idx_for_GER = find(peakCorrValues>0); 
    idx_unvoiced = find(peakCorrValues==0); 
    peaks_voicedSegments = peakCorrValues(idx_for_GER); 
    true_voicedFrames_Pitches = fs_Keele./peaks_voicedSegments;     
    %%For V/UV decision errors, we need also the 0 values 
    idx_for_VUV = find(peakCorrValues>=0); 
    peaks_allCertainSegments = peakCorrValues(idx_for_VUV); 
    true_allFrames_Pitches = fs_Keele./peaks_allCertainSegments; 
    true_allFrames_Pitches(isinf(true_allFrames_Pitches))=0; 
    
    for nT = 1:3% 1:3%%1:3%%1:noiseTypes 
        
        if nT==1, noise_name='babble', noisefilename='babbleTest_8KHz.wav'; end 
        if nT ==2, noise_name = 'factory', noisefilename = 'factoryTest_8KHz.wav'; end
        if nT ==3, noise_name = 'f16', noisefilename = 'f16Test_8KHz.wav'; end
        
        disp(['simulating for noise number ' num2str(nT)]);        
        
        for nS = 1:3%%1:nSnr
            snrdB = snr1(nS); 
            disp(['simulating for an snr' num2str(snrdB)]);      
            noiseVar = 10^(-snrdB/10)*signalPower/lengthClean;
            z1 = audioread(noisefilename);
            rndmStart = randi(150000); 
            z1 = z1(rndmStart:rndmStart+lengthClean-1); 
            z1 = sqrt(noiseVar)*z1/sqrt(z1'*z1/lengthClean);        
            y1 = s+z1; 
            
            %%Estimate for Cepstral, SHRP, RAPT and YIN 
            [pitchCEPST,nfrm] = cepstral(y1, fs, 60, 400, 0.01,2.0);
            pitchCEPST(pitchCEPST==0)=NaN;
%             [~,pitchSHRP,~,~]=shrp(y1,fs,f0BoundsHz,frame_length_ms,...
%                 timestep,SHR_Thresh, ceiling,med_smooth,1);
%            pitchSHRP(pitchSHRP==0)=NaN;
           [pitchRAPT,ttrapt] = estimateRAPT(y1,fs,'u');%%Include unvoiced frames
           periodYIN = yin2NEW(y1); 
           pitchYIN= fs./periodYIN.prd; 
            [~,vuvYIN,~,~]=SRH_PitchTracking(y1,fs,f0BoundsHz(1),f0BoundsHz(2)); 
            pitchYIN(vuvYIN==0)=NaN; %%to NaN since for GER
             
             
            total_prewht_frames = floor(lengthClean/M1);
            window = rectwin(N1); %%window for pitch-estimation
            clear lpc_mmse; clear g_mmse; clear lpc_nmf; clear g_nmf; 
            spectrum = fft(enframe(y1,window,M1),N1,2); 
            periodogram = spectrum.*conj(spectrum)/N1; 
            
            mmsePSD = estnoiseg(periodogram,M1/fs);
            mmseCov=ifft(mmsePSD');
               for m = 1:size(mmseCov,2) 
                   [lpc_mmse(:,m),g_mmse(m)]=levinson(mmseCov(:,m),lpc_prewh_order);
               end
             nmfPSD = computeNMF_PSD(periodogram,FullDictionary_2,...
                   numberSpeechVectors,lpc_mmse,numberNoiseVectors_2); 
              nmfCov = ifft(nmfPSD');
              for m = 1:size(nmfCov,2) 
                  [lpc_nmf(:,m),g_nmf(:,m)]=levinson(nmfCov(:,m),lpc_prewh_order);
              end
              y1_nmfPrew = prewhiteningTimeVariant(y1,N2,M2,lpc_nmf,g_nmf);
              idx = 1:N2;  
              for idf = 1:total_pitch_frames-4
                %  disp(['estimating frame ' num2str(idf) 'of ' num2str(...
                 %     total_pitch_frames)]);
                   [estf0_nmf(1,nF,nT,nS,idf),estOrd_nmf(1,nF,nT,nS,idf)]=...
                      f0Estimator.estimate(y1_nmfPrew(idx));
                  estf0_nmf(1,nF,nT,nS,idf) = estf0_nmf(1,nF,nT,nS,idf)*fs;
%                   [estf0_nmf(2,nF,nT,nS,idf),it(1,nF,nT,nS,idf)]=...
%                       reEstimatePitchConv(y1(idx),estf0_nmf(1,nF,nT,nS,idf)/fs,...
%                       estOrd_nmf(1,nF,nT,nS,idf),N1,M1,lpc_nmf(:,idf),...
%                       g_nmf(idf),y1_nmfPrew(idx),f0Estimator);
%                   estf0_nmf(2,nF,nT,nS,idf) = estf0_nmf(2,nF,nT,nS,idf)*fs;
                  [estf0_nmf(2,nF,nT,nS,idf),it(1,nF,nT,nS,idf)]=...
                      reEstimatePitchChoose(y1(idx),estf0_nmf(1,nF,nT,nS,idf)/fs,...
                      estOrd_nmf(1,nF,nT,nS,idf),N1,M1,lpc_nmf(:,idf),...
                      g_nmf(idf),y1_nmfPrew(idx),f0Estimator,1);
                  [estf0_nmf(3,nF,nT,nS,idf),it(2,nF,nT,nS,idf)]=...
                      reEstimatePitchChoose(y1(idx),estf0_nmf(1,nF,nT,nS,idf)/fs,...
                      estOrd_nmf(1,nF,nT,nS,idf),N1,M1,lpc_nmf(:,idf),...
                      g_nmf(idf),y1_nmfPrew(idx),f0Estimator,0);
                  estf0_nmf(2,nF,nT,nS,idf) = estf0_nmf(2,nF,nT,nS,idf)*fs;
                  estf0_nmf(3,nF,nT,nS,idf) = estf0_nmf(3,nF,nT,nS,idf)*fs;
                  idx = idx+M2;
              end          
                              
              %%Evaluate performance measures
              for kp = 1:3
                  [GER(kp,nF,nT,nS),u(kp)]=computeGPEAND(...
                      true_voicedFrames_Pitches(1:end-6),squeeze(estf0_nmf(kp,nF,...
                      nT,nS,idx_for_GER(1:end-6))),0.20,length(true_allFrames_Pitches));
              end
              [GER(4,nF,nT,nS),u(4)]=computeGPEAND(true_voicedFrames_Pitches(1:end-6),...
                  pitchCEPST(idx_for_GER(1:end-6)),0.20,length(true_allFrames_Pitches));
              [GER(5,nF,nT,nS),u(5)]=computeGPEAND(true_voicedFrames_Pitches,...
                 pitchRAPT(idx_for_GER),0.20,length(true_allFrames_Pitches));
               [GER(6,nF,nT,nS),u(6)]=computeGPEAND(true_voicedFrames_Pitches',...
                  pitchYIN(idx_for_GER(1:end)),0.20,length(true_allFrames_Pitches));
              
               estf0_nmf(isnan(estf0_nmf))=0; 
               pitchCEPST(isnan(pitchCEPST))=0;
               pitchRAPT(isnan(pitchRAPT))=0;    
               pitchYIN(isnan(pitchYIN))=0; 
               
               for kp = 1:3
                   VUV(kp,nF,nT,nS)=computeVUVErrorRatio(...
                       true_allFrames_Pitches(1:end-6),squeeze(estf0_nmf(...
                       kp,nF,nT,nS,idx_for_VUV(1:end-6))));
               end
               VUV(4,nF,nT,nS)=computeVUVErrorRatio(true_allFrames_Pitches(...
                   1:end-6),pitchCEPST(idx_for_VUV(1:end-6)));
               VUV(5,nF,nT,nS)=computeVUVErrorRatio(true_allFrames_Pitches(...
                   1:end-3),pitchRAPT(idx_for_VUV(1:end-3)'));
               VUV(6,nF,nT,nS)=computeVUVErrorRatio(true_allFrames_Pitches(...
             1:end-3),pitchYIN(idx_for_VUV(1:end-3))');
              
              for kp = 1:6
              FFE(kp,nF,nT,nS) = u(kp)+VUV(kp,nF,nT,nS)';
              end
              FFE(:,:,:,1)
              
        end
    end
end

FFE_all = mean(FFE,2); GER_all = mean(GER,2); VUV_all = mean(VUV,2); 
%%reshape(VUV_all(:,1,2,:),6,3)'

GER_dev = std(GER,0,2); VUV_dev = std(VUV,0,2); FFE_dev = std(FFE,0,2); 
ts = tinv([0.025 0.975],nF-1); 
GER_ConfInt = ts(2)*GER_dev/sqrt((nF)); 
VUV_ConfInt = ts(2)*VUV_dev/sqrt((nF)); 
FFE_ConfInt = ts(2)*FFE_dev/sqrt((nF)); 

%%Now time to plot: 
% for nT = 1:3 
%     if nT == 1
%         title = 'Babble noise'; 
%     end
%     if nT==2
%         title = 'Factory noise'; 
%     end
%     if nT==3
%         title = 'F16 noise';
%     end
%     plotbarErrors(reshape(GER_ConfInt(:,1,nT,:),6,3)',...
%         reshape(GER_all(:,1,nT,:),6,3)',title,'iSNR',...
%         {'-5 dB' '5 dB' '15 dB'},'GER(%)')
% end
%     
%%Nice plots for the ICASSP paper 
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/3 scrsz(3)/4 scrsz(4)/3]);
for nT = 1:3
    subplot(3,3,nT);
    if nT == 1        title = 'Babble noise';     end
    if nT==2        title = 'Factory noise';     end
    if nT==3        title = 'F16 noise';    end
    plotbarErrors(reshape(GER_ConfInt(:,1,nT,:),6,3)',...
        reshape(GER_all(:,1,nT,:),6,3)',title,'iSNR',...
        {'-5 dB' '5 dB' '15 dB','Fontsize', 12, 'interpreter', 'latex' },'GER(\%)')
    ylim([0 66])
    subplot(3,3,nT+3)
    if nT == 1        title = 'Babble noise';     end
    if nT==2        title = 'Factory noise';     end
    if nT==3        title = 'F16 noise';    end
    plotbarErrors(reshape(VUV_ConfInt(:,1,nT,:),6,3)',...
        reshape(VUV_all(:,1,nT,:),6,3)',title,'iSNR',...
        {'-5 dB' '5 dB' '15 dB','Fontsize', 12, 'interpreter', 'latex' },'VDE(\%)')
     ylim([0 48])
    subplot(3,3,nT+6)
    if nT == 1        title = 'Babble noise';     end
    if nT==2        title = 'Factory noise';     end
    if nT==3        title = 'F16 noise';    end
    plotbarErrors(reshape(FFE_ConfInt(:,1,nT,:),6,3)',...
        reshape(FFE_all(:,1,nT,:),6,3)',title,'iSNR',...
        {'-5 dB' '5 dB' '15 dB','Fontsize', 12, 'interpreter', 'latex' },'FFE(\%)')
    ylim([0 78])
  end
  leg = legend('NLS-NMF','NLS-NMF Iter1','NLS-NMF Iter2','Cepstrum','RAPT','YIN',...
        'orientation','horizontal');
    set(leg,'interpreter','latex', 'fontsize', 11)
    set(gcf,'color','w');
    
    

% scrsz = get(0,'ScreenSize');
% figure('Position',[scrsz(3)/4 scrsz(4)/3 scrsz(3)/4 scrsz(4)/3]);
%    subplot(2,1,1); plot(true_allFrames_Pitches(1300:1950),'or',...
%       'linewidth',2); hold; plot(squeeze(estf0_nmf(1,1,2,3,idx_for_VUV(...
%       1300:1950))),'*k','linewidth',1.5); xlim([0 651])
%    xlabel('frame index','Fontsize', 11)
%    ylabel('fund. freq (Hz)','Fontsize', 11)
%       subplot(2,1,2); plot(true_allFrames_Pitches(1300:1950),'or',...
%       'linewidth',2); hold; plot(squeeze(estf0_nmf(3,1,2,3,idx_for_VUV(...
%       1300:1950))),'*k','linewidth',1.5); xlim([0 651])
%     xlabel('frame index','Fontsize', 11)
%       ylabel('fund. freq (Hz)','Fontsize', 11)
%       leg = legend('ground truth','estimated fundamental frequency',...
%           'orientation','horizontal');
%  % set(leg,'interpreter','latex', 'fontsize', 9);
%   set(gcf,'color','w');
%  