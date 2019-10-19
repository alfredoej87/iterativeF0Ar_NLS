function yy = prewhiteningTimeVariant(FullSignal,...,
    segmentLength,nShift,estimatedParametersMtx,estimVars); 

%window = rectwin(segmentLength); 
window = sqrt(hanning(segmentLength,'periodic')); 
% window = hanning(segmentLength,'periodic'); 
%window = window/sqrt(2*mean(window.^2)); %%why is the 2 necessary?????  
STFT_noisy = fft(enframe(FullSignal,window,nShift),segmentLength,2); 

for ab = 1:size(estimatedParametersMtx,2) 
STFT_segment = STFT_noisy(ab,:);
PrewFilter_STFT = (fft(estimatedParametersMtx(:,ab),segmentLength)); %%%/sqrt(estimVars(ab)); 
%%PrewFilter_STFT = ones(segmentLength,1); %%to verity output is the same
STFT_filtered(ab,:) = STFT_segment.*PrewFilter_STFT.';
end

prewhitenedLongSignal_Temp= overlapadd(ifft(STFT_filtered,segmentLength,...
    2),window,nShift); 
yy = prewhitenedLongSignal_Temp; 
%%%Need to apply the same normalization before and after pre-whitening 
