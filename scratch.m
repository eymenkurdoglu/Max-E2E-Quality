close all
clc
videos = {'CREW','CITY','FOREMAN','HARBOUR','ICE','AKIYO'};

figure; hold all;
xlabel('reference distance')
ylabel('P-frame size / I-frame size')

for v = 1:length(videos)

path = ['~/Dropbox/Matlab/3_OptimalQuality/data/',videos{v},'-352x288-30-32-3/'];
l = []; qp = []; rates = [];
for file = dir( strcat(path,'*.txt') )'
    contents = dlmread( strcat(path,file.name) );
    l  = [l contents(:,1)];
    qp = [qp contents(:,2)];
    dum = strsplit(file.name,'.');
    rates = [rates str2double(dum{1})];
end
[rates,dum] = sort(rates);
qp = qp(:,dum);
l = l(:,dum);

L = zeros(4,length(rates));
for r = 1:length(rates)
    count = zeros(4,1);
    for i = 1:size(l,1)
        if mod(i,32) == 1
            L(1,r) = L(1,r)+l(i,r);
            count(1) = count(1)+1;
        elseif mod(i,4) == 1
            L(2,r) = L(2,r)+l(i,r);
            count(2) = count(2)+1;
        elseif mod(i,4) == 3
            L(3,r) = L(3,r)+l(i,r);
            count(3) = count(3)+1;
        else
            L(4,r) = L(4,r)+l(i,r);
            count(4) = count(4)+1;
        end
    end
end
L(1,:) = L(1,:)/count(1);
L(2,:) = L(2,:)/count(2);
L(3,:) = L(3,:)/count(3);
L(4,:) = L(4,:)/count(4);

for r = 1:length(rates)
    L(:,r) = L(:,r)/L(1,r);
end
% bar3(L)

%     figure; hold all;
%     plot(rates,L(2,:))
%     plot(rates,L(3,:))
%     plot(rates,L(4,:))
%     xlabel('video rate')
%     ylabel('P-frame size / I-frame size')
%     title(videos{v})
%     legend('Base layer', 'Layer 2', 'Layer 3')

    plot( [4, 2, 1], mean(L(2:4,5:end),2) )

end
legend(videos,'Location','NorthWest')