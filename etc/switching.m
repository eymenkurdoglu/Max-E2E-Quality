dbstop if error
close all
clear all
clc

hP = 1;
k = zeros(32,1);
if hP
    k(1:4:end) = 25; k(1) = 70;
    k(3:4:end) = 25;
    k(2:2:end) = 25;
else
    k(1:end) = 25; %k(1) = 70;
end
k = k(:);
eps = 0.05;
f = round( 2 * sum(k) * eps );

[allocs, scores] = greedy_fec_search2( f, k, eps );

cw_lengths = allocs + repmat(k,1,f);
fec_rates = allocs./cw_lengths;

figure; stairs(allocs(:,end))
figure; plot(scores(:))
figure; stairs(fec_rates(:,end))

% F = allocs(:,end);

phi = zeros(32,f);
for j = 1:f
    F = allocs(:,j);
    for i = 1:32
        phi(i,j) = sum( binopdf(0:F(i),F(i)+k(i),eps) );
    end
end
figure; stem(phi(:,end))
% pr = ones(32,length(f)); pr(1) = phi(1);
% for i = 2:32
%     pr(i) = pr(i-1) * phi(i);
% end
% figure; plot(pr)

% figure;
% j=1;
% for i = 1:4:32
%    subplot(10,1,j)
%    stairs(fec_rates(i,:))
%    ylim( [0 0.15] )
%    j = j+1;
% end
% subplot(10,1,j)
% stairs(sum(allocs(3:4:end,:))./sum( cw_lengths(3:4:end,:) ))
% subplot(10,1,j+1)
% stairs(sum(allocs(2:2:end,:))./sum( cw_lengths(2:2:end,:) ))

F = zeros(32,1);
F(1:4:end) = 2; F(1) = 4;
% F(1) = 5; F(5) = 2; F(9) = 1; F(3) = 1;%F(13) = 1;
phi_ = zeros(32,1);
for i = 1:32
    phi_(i) = sum( binopdf(0:F(i),F(i)+k(i),eps) );
end
calc_score( phi_, 1 )