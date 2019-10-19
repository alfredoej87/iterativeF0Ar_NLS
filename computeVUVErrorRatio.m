
%%function to compute the classification errors

function VoiceDetectionError = computeVUVErrorRatio(trueValues,...
    estimatedValues) 

 positions = find((trueValues==0 & estimatedValues>0)|...
     (trueValues>0 & estimatedValues==0)); 
 VoiceDetectionError = (length(positions)/length(trueValues))*100; 
end
