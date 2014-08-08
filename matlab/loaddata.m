function H=loaddata(filename)

fileID = fopen(filename);
if fileID == -1
    error(['Cannot open file "' filename '"'])
end
m = fread(fileID, 1, 'uint64');
n = fread(fileID, 1, 'uint64');
H = fread(fileID, [m, n], 'double');
fclose(fileID);

H = H';
