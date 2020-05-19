function [GPE,scaledGPE]=computeGPEAND(trueValues,estimatedValues, percent,N)

  nanIndicator = ~(isnan(trueValues)) & ~(isnan(estimatedValues));
  scaledGPE= length(find(abs(estimatedValues(nanIndicator)-...
        trueValues(nanIndicator))>trueValues(nanIndicator)*percent))/N*100; 
    GPE = mean(abs(estimatedValues(nanIndicator)-...
        trueValues(nanIndicator))>trueValues(nanIndicator)*percent)*100; %
  end
