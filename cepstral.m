function [f0_cepst,nfrm] = cepstral(x, fs, f0min, f0max, timestep, pthr1)
%
% Cesptrum pitch detection with automatic voice (music) marker detection.
%
% Input arguments:
% x -> the audio data
% fs -> sample rate for the audio data (unit: Hz)
% f0min -> lower bound of pitch (unit: Hz)
% f0max -> higher bound of pitch (unit: Hz)
% timestep -> the time offset of the detected pitch (unit: sec)
% pthr1 -> amplitude ratio of the first two highest Cepstrum peaks, used for
% voiced/unvoiced detection. A higher value gives a more stringent threshold
% for voiced frame detection. If pthr1=1, all frames are regarded as voiced.
% pthr1 is a value greater than or equal to 1.
%
% Ouput:
% f0_cepst -> the detected Cepstrum pitch for each frame.
% nfrm -> number of frames.
%
% See the Bridge project website at Wireless Communication and Networking
% Group (WCNG), University of Rochester for more information: 
% http://www.ece.rochester.edu/projects/wcng/project_bridge.html
%
% Code modified by He Ba and Na Yang, University of Rochester.
% Originally provided by Lawrence R. Rabiner, Rutgers University and
% University of California at Santa Barbara.
% Copyright (c) 2013 University of Rochester.
% Version March 2013.
%
% Permission to use, copy, modify, and distribute this software without 
% fee is hereby granted FOR RESEARCH PURPOSES only, provided that this
% copyright notice appears in all copies and in all supporting 
% documentation, and that the software is not redistributed for any fee.
%
% For any other uses of this software, in original or modified form, 
% including but not limited to consulting, production or distribution
% in whole or in part, specific prior permission must be obtained from WCNG.
% Algorithms implemented by this software may be claimed by patents owned 
% by WCNG, University of Rochester.
%
% The WCNG makes no representations about the suitability of this 
% software for any purpose. It is provided "as is" without express
% or implied warranty.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Default input parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default ratio of the amplitude of the first two highest Cepstrum peaks.
if nargin <6
    pthr1 = 1.1; 
end

% default time offset of detected pitch
if nargin <5
    timestep = 0.01; 
end

% default higher bound of pitch
if nargin <4
    f0max = 600; 
end

% default lower bound of pitch
if nargin <3
    f0min = 50; 
end

if nargin <2
    error('Lack of input arguments!') % there should be at least 2 arguments (audio data and fs) 
end

if (size(x,1)<size(x,2))
    x=x';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Parameter settings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_per_window = 0.06;                             % frame length in second
Lmed = 5;                                           % median filter window length

point_per_window = round(time_per_window*fs);       %points per window
point_per_timestep=round(timestep*fs);              %point_per_timestep setting 

ls=length(x);                                       %length of input audio data
frame_center=ceil(point_per_window/2);              %the center point of current frame
nfrm = floor((ls-point_per_window)/point_per_timestep)+1; % number of frames
f0_cepst = zeros(nfrm,1); % cepstrum pitch decision

p1=[]; % highest pitch level
p2=[]; % 2nd highest pitch level
pd1=[]; % highest pitch level's pitch period
pd2=[]; % 2nd highest pitch level's pitch period

i=1;                                                %frame count
while (frame_center+point_per_window/2 <= ls)   
    xframe_orig=x(max(frame_center-ceil(point_per_window/2),1):frame_center+ceil(point_per_window/2),1); %original frame  
    xframe = xframe_orig-mean(xframe_orig); %zero mean
    xframe1 = filter([1 -1],1,xframe); % high-pass filter
    frame_fft = abs(fft(xframe1.*hann(length(xframe1)),2^16)); %2^16 points fft

    ppdlow=round(fs/f0max);
    ppdhigh=round(fs/f0min);
    xframe1t=log(frame_fft);
    xframec=ifft(xframe1t,2^16);
       
    % initialize local frame cepstrum over valid range from ppdlow to ppdhigh
    indexlow=ppdlow+1;
    indexhigh=ppdhigh+1;
    loghsp=xframec(indexlow:indexhigh);

    % find cepstral peak location (ploc) and level (pmax) and save results in
    % pd1 (for ploc) and p1 (for pmax)
    pmax=max(loghsp);
    ploc1=find(loghsp == pmax);
    ploc=ploc1+ppdlow-1;
    f0_cepst(i,1) = fs/ploc;
    p1=[p1 pmax];
    pd1=[pd1 ploc];

    % eliminate strongest peak in order to find highest secondary peak 
    % which is spaced away from primary peak
    % save secondary peak in pd2 and secondary level in p2
    n1=max(1,ploc1-4);
    n2=min(ploc1+4,length(loghsp));
    loghsp2=loghsp;
    loghsp2(n1:n2)=0;
    pmax2=max(loghsp2);
    p2=[p2 pmax2];
    ploc2=find(loghsp2 == pmax2);
    ploc2=ploc2+ppdlow-1;
    pd2=[pd2 ploc2];
    
    i=i+1;
    frame_center = frame_center + point_per_timestep;
end

% voice marker based on ratio of first two highest cepstrum peaks 
ppdf=smoothpitch(pd1,pd2,p1,p2,pthr1); % cepstrum pitch period
%p1m=medf(ppdf,Lmed,length(pd1)); % median filter cepstrum
for myf = 1:length(ppdf)%%(p1m)
    if ppdf(myf)==0%%%p1m(myf) == 0
        f0_cepst(myf) = 0;
    else
        f0_cepst(myf) = fs/ppdf(myf); %fs/p1m(myf); % cepstrum pitch
    end
end

% for myf = 1:length(pd1)
%     if pd1(myf) == 0
%         f0_cepst(myf) = 0;
%     else
%         f0_cepst(myf) = fs/pd1(myf); % cepstrum pitch
%     end
% end
