function tauRR = tau_recurrence_rate(RP)
% TAU_RECURRENCE_RATE computes the tau recurrence rate of the given
% recurrence plot RP

N = size(RP,1);
M = size(RP,2);
assert(N==M, "Input needs to be a square matrix.")
tauRR = zeros(1,N);
cnt = 1;
for k = 0:N-1
    d = diag(RP,k);
    tauRR(cnt) = length(find(d))/length(d);
    clear d
    cnt = cnt + 1;
end

end