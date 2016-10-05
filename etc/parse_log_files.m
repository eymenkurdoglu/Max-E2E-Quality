function [Lengths, QPs, Target_BitRates] = parse_log_files( path ) 

Lengths = [];
QPs = [];
Target_BitRates = [];

for file = dir( strcat(path,'*.txt') )'
    contents = dlmread( strcat(path,file.name) );
    Lengths  = [Lengths contents(:,1)];
    QPs = [QPs contents(:,2)];
    dum = strsplit(file.name,'.');
    Target_BitRates = [Target_BitRates str2double(dum{1})];
end

[Target_BitRates,ix] = sort(Target_BitRates);
Lengths = Lengths(:,ix);
QPs = QPs(:,ix);

return