function [FSz, QPs, R] = parse_log_files( path ) 
% parse_log_files       Parse the output of the x264 debug info (already
% processed with awk)
%  This function reads all files with extension .txt on the given path. All
%  such files should contain, for each compresssed frame, the size in
%  bytes on the left column and the mean QP on the right column. The
%  target bitrate is part of the file name. 
%  FSz & QPs: rows <-> frames, columns <-> target bitrates
%  R: target bitrates
% 
%  Example: [size,qp,target] = parse_log_files( '~/Matlab/CITY-352x288-30.0-32-3/' )

files = dir( strcat(path,'*.txt') )';

numFiles = length(files);
numLines = size(dlmread( strcat(path,files(1).name) ),1);

FSz = zeros(numLines,numFiles);
QPs = zeros(numLines,numFiles);
R = zeros(1,numFiles);

for i = 1 : numFiles
    file = files(i);
    contents = dlmread( strcat(path,file.name) );
    FSz(:,i)  = contents(:,1);
    QPs(:,i) = contents(:,2);
    dum = strsplit(file.name,'.');
    R(i) = str2double(dum{1});
end

[R,ix] = sort(R);
FSz = FSz(:,ix);
QPs = QPs(:,ix);

return