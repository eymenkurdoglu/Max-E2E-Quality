dbstop if error

% load the encoded frame size matrix (crew) with x264 using our settings
load frame.mat

% for demonstration, we choose the 4th intra-period
intra_period = 4;
packet_size = 200; %bytes
bw_index = 12;

% extract the frame sizes in this intra-period as payloads of packet_size bytes
% should be predicted
frames_in_ip = frame(32*intra_period-31 : 32*intra_period, :); % frame lengths (in bytes) for each source rate
sr_in_ip = 0.008*sum(frames_in_ip)/(16/15); % source rate in this intra-period (kbps), could be higher than the intended source rate 
frames_in_ip = ceil(frames_in_ip/packet_size); % frame lengths (in packets)
frames_in_ip = frames_in_ip(:,1:bw_index);
abw = max(sum(frames_in_ip));

% number of FEC blocks available
f = abw-sum(frames_in_ip);

scores = zeros(bw_index,1);
SCORES = zeros(bw_index,6);

alloc = zeros(32,bw_index);

% average qp values for each source rate (sadly averaged over all intra-periods), should be predicted
qp = [43.6871875;
      36.95625;
      33.9309375;
      31.8509375;
      30.2390625;
      28.9309375;
      27.79625;
      26.8228125;
      25.986875];

h1 = figure; hold all;
h2 = figure; hold all;  
tic
for alpha = [0.001 0.01 0.1]

    beta = 9*alpha;

    % search!
    for i = 1:bw_index % for each source rate
            [best_alloc, best_obj] = greedy_fec_search3( f(i), frames_in_ip(:,i), alpha, beta );
            scores(i) = best_obj;
            alloc(:,i) = best_alloc;
    end

    Q = mnqq(qp,4.51,20).*mnqt(scores*15/16,3.09,30);

    figure(h1);
    % SCORES(:,1+p/0.05) = scores(:);
    plot(100*sum(frames_in_ip)/abw,scores);

    figure(h2);
    plot(100*sum(frames_in_ip)/abw,Q);

end
toc
ylabel('Quality (MNQQ*MNQT)')
xlabel('Source rate percentage (%)')