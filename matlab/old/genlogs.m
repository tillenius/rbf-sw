function genlogs(str)

%  system('mex CXXOPTIMFLAGS='-DNDEBUG -mavx -Wall -Ofast -ffast-math'  evalCartRhs_fd_mex.cpp -largeArrayDims -I/home/martin/git/tl2/tasklib')

for numthreads=[0 3]
for chunksize = [6400 3200 2130  1600  1280 1066 914 800 640 400]
for i=0:4
  name = ['results4/cpp-2d' str '-' num2str(chunksize) '-' num2str(numthreads) '-' num2str(i) '.dat'];
  disp(name)
  mt
  system(['cp rbf-cpp.log ' name])
end
end
end


