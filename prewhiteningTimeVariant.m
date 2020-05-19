function yy = prewhiteningTimeVariant(FullSignal,...,
    segmentLength,nShift,estimatedParametersMtx,estimVars); 

window = sqrt(hanning(segmentLength,'periodic')); 
STFT_noisy = fft(enframe(FullSignal,window,nShift),segmentLength,2); 

for ab = 1:size(estimatedParametersMtx,2) 
STFT_segment = STFT_noisy(ab,:);
PrewFilter_STFT = (fft(estimatedParametersMtx(:,ab),segmentLength));  
STFT_filtered(ab,:) = STFT_segment.*PrewFilter_STFT.';
end

prewhitenedLongSignal_Temp= overlapadd(ifft(STFT_filtered,segmentLength,...
    2),window,nShift); 
yy = prewhitenedLongSignal_Temp; 

