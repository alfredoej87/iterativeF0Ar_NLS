function [estNewf0,i]=reEstimatePitchChoose(y,estF0,L,N,M,lpcPrewParam,gPrew,...
    prew_y,f0Estimator,flag)

fs = 8000; 
estNewf0 = estF0; %%for the 1st iteration 
LNew = L; %%for the 1st iteration
diff = Inf;
i = 0; 
noMaxIter = 8; 
p1 = 25; %%pre-whitening order; 
small_delta = 1e-5; %%regularization for APES 
%f0Estimator = fastF0Nls(N,27,[65 370]/8000); 

%%Uncomment when considering inharmonicities
%  D_range = 25; %%(Hz), since the smallest value of f0 is 57 Hz
%  D_range_rel = D_range/(4000)*2*pi;
%  stepsize = 1; %%(Hz);
%  stepsize_rel = stepsize./(4000)*2*pi;
% noIterFib = 15; 

while (i<noMaxIter) & (diff>1e-6) %repeat regardless 
  if flag ==1 
     if LNew>0
        Dnls = zeros(1,LNew); %%search of perturbations from pre-whitened data
        expMtx = ZDreal(2*pi*estNewf0,N,LNew,Dnls);        
        J_old=real((prew_y'*expMtx)*inv(expMtx'*expMtx)*(expMtx'*prew_y));        
        ampsLS = inv(expMtx'*expMtx)*expMtx'*y;
        ampsLSReal = 2*abs(ampsLS(1:LNew));
        phasesPlus = angle(ampsLS(1:LNew));
        sHat=generatePeriodicSignal(N,estNewf0,LNew,ampsLSReal,phasesPlus,...
            Dnls);
    else
        sHat=zeros(1,N);
        diff = 0; 
        J_old = Inf; %%to expect if the frame was unvoiced it remains:) !!! 
    end
    else %%Odd iterations do APES-filtering
        if LNew>0
        lpcPrewParam = [lpcPrewParam(:);zeros(N-size(lpcPrewParam(:),1),1)];
        Dnls = zeros(1,LNew);
       %  Dnls = D_nls_real(prew_y,2*pi*estNewf0,D_range_rel,...
        %             stepsize_rel,LNew,noIterFib);
        expMtx = ZDreal(2*pi*estNewf0,N,LNew,Dnls);       
        J_old=real((prew_y'*expMtx)*inv(expMtx'*expMtx)*(expMtx'*prew_y));      
        toep1 = tril(toeplitz(lpcPrewParam'));
        toep2 = tril(toeplitz([0 lpcPrewParam(end:-1:2)']));
        invQ_apes = (toep1*toep1'-toep2*toep2')/gPrew;
        H_apes = real(invQ_apes*expMtx*inv(expMtx'*invQ_apes*...
            expMtx+small_delta*eye(2*LNew))*expMtx');
        sHat=H_apes'*y;
        else
        sHat=zeros(1,N);
        diff = 0; 
        J_old = Inf; %%to expect if the frame was unvoiced it remains:) !!! 
        end
    end
    zHat = y-sHat(:); 
    lpcPrewParam = lpc(zHat,p1);
    prew_y = filter(lpcPrewParam,1,y); 
    [estNewf0,LNew] = f0Estimator.estimate(prew_y); 
    %%The next step should be to estimate in-harmonicities 
    Dnls = zeros(1,LNew);     
  %  Dnls=D_nls_real(prew_y,2*pi*estNewf0,D_range_rel,stepsize_rel,...
   %         LNew,noIterFib);
    %%Compute cost function  (only if frame was detected as voiced)
    if LNew>0
        expMtx = ZDreal(2*pi*estNewf0,N,LNew,Dnls);
        J_new = real((prew_y'*expMtx)*inv(expMtx'*expMtx)*(expMtx'*prew_y));
    else
        J_new = Inf;
    end    
    diff = abs(J_new-J_old)/J_old;   
i = i+1;
end
