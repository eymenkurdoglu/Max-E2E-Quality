function P = star_param( video, ipr )

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
elseif  strcmp(video, 'FOREMAN')    
    beta_q=1.087917; beta_t=0.558064; Rmax=1171.495200; qmin=9.878532;
    alpha_q = 4.57; alpha_t = 3.80;
end
tmax = 30;

P = struct('beta_q',beta_q,'beta_t',beta_t,'Rmax',Rmax,'qmin',qmin,...
    'alpha_q',alpha_q,'alpha_t',alpha_t,'tmax',tmax,'ipr',ipr);

return