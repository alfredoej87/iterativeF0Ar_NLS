
function Zk = ZDreal(omega_k,M,L,Delta_l)
m = (0:M-1)';
Zk = zeros(M,2*L);
if length(Delta_l) == L
    psi = omega_k*(1:L)+Delta_l;%this to be used for positive and negative
else
    error('Length of Delta_l does not conform to model order')
end

for ll = 1:L 
    Zk(:,ll) = exp(1j*psi(ll)*m); 
    Zk(:,ll+L) = exp(-1j*psi(ll)*m); 
end
end

% if strcmp(dir,'back')
%     for ll = 1:L;
%         Zk(:,ll) = exp(-1j*psi(ll)*m);
%     end
% else
%     for ll = 1:L;
%         Zk(:,ll) = exp(1j*psi(ll)*m);
%     end
% end
% end
% 
% 
% for ll = 1:L;
%     Zk(:,ll) = exp(1j*ll*(omega_k+khalf*m).*m);
%     Zk(:,ll+L) = exp(-1j*ll*(omega_k+khalf*m).*m);
% end
