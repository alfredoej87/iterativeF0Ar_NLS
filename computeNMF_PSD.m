function [NMFPSDhat, NMFPSDhat2, speechFrameActiv, noiseFrameActiv] = ...
    computeNMF_PSD(NoisyPSD,inputDict,numSpeechVct,MMSEfromLPC,numNoiseVct)

TotalFrames = size(NoisyPSD,1);
segmentLength = size(NoisyPSD,2);
alpha = 0.875;

for tt=1:TotalFrames
    currentNoisyPSD = NoisyPSD(tt,:)';
    MMSELPC = MMSEfromLPC(:,tt);
    MMSEPSD = computeArPsd(segmentLength,MMSELPC,1);
       currDict = [inputDict MMSEPSD];
    if tt ==1
    all_actC(tt,:) = ActivationFunction(currentNoisyPSD,currDict);
      else 
        all_actC(tt,:) = ActivationFunction(currentNoisyPSD,currDict,...
            all_actC(tt-1,:)'); 
    end    
    sp_actC = all_actC(tt,1:numSpeechVct);
    noise_actC = all_actC(tt,numSpeechVct+1:end);
    speechPower = currDict(:,1:numSpeechVct)*sp_actC';
    noisePower = currDict(:,numSpeechVct+1:end)*noise_actC';
   % aprioriSNR = speechPower./noisePower; 
    aprioriSNR = speechPower./max(noisePower,1e-12); 
    %%%%%aprioriSNR = max(10^(-15/10),aprioriSNR); 
    wghtPSD = (((1./(1+aprioriSNR))).^2).*currentNoisyPSD;
    wghtNoiseVar = (aprioriSNR./(1+aprioriSNR)).*noisePower;
    NMFPSDhat(tt,:) = wghtPSD+wghtNoiseVar;   
    NMFPSDhat2(tt,:) = NMFPSDhat(tt,:); 
  if tt >1
     NMFPSDhat2(tt,:)= alpha*NMFPSDhat2(tt-1,:)+(1-alpha)*NMFPSDhat2(tt,:);
   end
    if nargout>2
    speechFrameActiv(tt,:) = sp_actC;
    noiseFrameActiv(tt,:) = noise_actC;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ActivationCoeff = ActivationFunction(noisyPSD,AR_PSD_BasisMtx,...
    prev_Frame_Actv_Coeff)

if nargin<3
ActivationCoeff = rand(size(AR_PSD_BasisMtx,2),1);
numberIterations = 35;
else
   numberIterations = 35; %%30; 
    ActivationCoeff = rand(size(AR_PSD_BasisMtx,2),1); %%prev_Frame_Actv_Coeff; 
    %%ActivationCoeff = prev_Frame_Actv_Coeff; 
end 
for times = 1:numberIterations
    Factor=AR_PSD_BasisMtx'*(((AR_PSD_BasisMtx*ActivationCoeff).^(-2)).*...
        noisyPSD)./(AR_PSD_BasisMtx'*(AR_PSD_BasisMtx*...
        ActivationCoeff).^(-1));
    ActivationCoeff = ActivationCoeff.*Factor;    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function arPsd = computeArPsd(nData, arParameters, exVar)
    arPsd = exVar./abs(fft([arParameters(:)], nData)).^2;
end
    