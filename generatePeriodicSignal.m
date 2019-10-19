function periodicSignal = generatePeriodicSignal(nData, f0, modelOrder, ...
        amps, phases,Dnls)
    if modelOrder~=0
    %complexSinusMtx = exp(1i*2*pi*f0*(0:nData-1)'*(1:modelOrder));
    complexSinusMtx = ZDc(2*pi*f0,nData,modelOrder,Dnls); 
    complexAmplitudes = amps.*exp(1i*phases);
    periodicSignal = real(complexSinusMtx*complexAmplitudes);
    else
        periodicSignal = zeros(nData,1); 
    end