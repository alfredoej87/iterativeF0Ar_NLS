function [GPE,scaledGPE]=computeGPEAND(trueValues,estimatedValues, percent,N)
    % find the indices where both signals are free of NaNs
    %%See Matlab Mathworks definition: 
    %%https://se.mathworks.com/help/audio/ref/pitch.html
  nanIndicator = ~(isnan(trueValues)) & ~(isnan(estimatedValues));
   %nanIndicator = ~(isnan(trueValues) & isnan(estimatedValues));
    %nanIndicator = ~isnan(trueValues); 
    scaledGPE= length(find(abs(estimatedValues(nanIndicator)-...
        trueValues(nanIndicator))>trueValues(nanIndicator)*percent))/N*100; 
    GPE = mean(abs(estimatedValues(nanIndicator)-...
        trueValues(nanIndicator))>trueValues(nanIndicator)*percent)*100; %
    %%%%%%%%scaledGPE = proportion*GPE; 
end