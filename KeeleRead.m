% reads the microphone/laryngograph waveform from keele database
% By calling KeeleRead you will have a lot of fun.
function [vSignal, dFs] = KeeleRead(szFile, varargin)

% File   :  KeeleRead.m
% Author :  Robert Rehr <r.rehr AT uni-oldenburg.de>
% Date   :  24.10.2012
%
% Updates:

% try to open audio file (matlab assumes binary file by default)
[fWaveFile, szError] = fopen(szFile, 'r+');

% in case file could not be opened display error message
if fWaveFile == -1
    % return error message
    error(['Could not open wave file (Reason: ' szError ')'])
end

if any(strcmpi('size', varargin))
    % get only size of file
    fseek(fWaveFile, 0, 'eof');
    vSignal = [ftell(fWaveFile)/2 1];
    
    % close file handle
    fclose(fWaveFile);
else
    % read samples as 16 bit integers coded as little endian 
    % (keele documentation mentions intel format)
    vSignal = fread(fWaveFile, inf, 'int16', 0, 'l');
    % close file handle
    fclose(fWaveFile);

    % normalize signal 2^15 = 32768 = maximum value of a sample
    vSignal = vSignal ./ 32768;
end

% sample rate is always 20 kHz
dFs = 2e4;

% End of KeeleRead.m
