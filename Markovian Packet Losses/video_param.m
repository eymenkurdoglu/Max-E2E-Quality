function Param = video_param( video, ipr, numTempLayers, N )

if      strcmp(video, 'CREW')
    beta_q=1.061165; beta_t=0.707289; Rmax=1188.956800; qmin=14.100638;
    alpha_q = 4.51; alpha_t = 3.09;
elseif  strcmp(video, 'CITY')
    beta_q=1.142239; beta_t=0.471515; Rmax=1150.621600; qmin=10.662216;
    alpha_q = 7.25; alpha_t = 4.10;
elseif  strcmp(video, 'HARBOUR')
    beta_q=1.320098; beta_t=0.584034; Rmax=1196.192000; qmin=21.765459;
    alpha_q = 9.65; alpha_t = 2.83;
elseif  strcmp(video, 'ICE')
    beta_q=0.868087; beta_t=0.653480; Rmax=943.812800; qmin=5.318390;
    alpha_q = 5.61; alpha_t = 3.00;
elseif  strcmp(video, 'FOREMAN') %aq=1.087917, at=0.558064, Rmax=1171.495200, qmin=9.878532
    beta_q=1.087917; beta_t=0.558064; Rmax=1171.495200; qmin=9.878532;
    alpha_q = 4.57; alpha_t = 3.80;
elseif  strcmp(video, 'FG') %aq=1.155468, at=0.570933, Rmax=1201.996000, qmin=22.806485
    beta_q = 1.155468; beta_t = 0.570933; Rmax = 1201.996; qmin = 22.806485;
    alpha_q = 10.68; alpha_t = 2.8;    
end
tmax = 30;

numDescendants = zeros(1,N);
G = 2^(numTempLayers-1);
ix_so_far = 1:G:N;

for i = ix_so_far
    numDescendants(i) = N-i;
end

while length(ix_so_far) < N
    G = G/2;
    numDescendants( ix_so_far + G ) = G-1;
    temp = [ix_so_far, ix_so_far + G];
    ix_so_far = temp;
end

layerMap = zeros(N,1);
for layer = numTempLayers : -1 : 1   
    layerMap(mod(0:N-1,2^(numTempLayers-layer))==0) = layer;
end

referenceMap = zeros(N,1);
for i = 1:N
    if mod(i,4) == 1
        referenceMap(i) = i-4;
    elseif mod(i,4) == 2 || mod(i,4) == 0
        referenceMap(i) = i-1;
    else
        referenceMap(i) = i-2;
    end
end

Param = struct('beta_q',beta_q,'beta_t',beta_t,'Rmax',Rmax,'qmin',qmin,...
    'alpha_q',alpha_q,'alpha_t',alpha_t,'tmax',tmax,...
    'ipr',ipr,'numTempLayers',numTempLayers,'numDescendants',numDescendants,'layerMap',layerMap,...
    'referenceMap',referenceMap);

return