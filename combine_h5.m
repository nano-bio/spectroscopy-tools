% this file reads all *.h5 files in a folder and combines them to a single
% one, where the writes are the sum-spectra of the original files

% new file to use
newh5 = 'new.h5';

% output file
outputfile = 'filelist_864nm.txt';

% directory to combine
directory = '864nm';
h5files = dir([directory '/*.h5']);

newbuffer = 1;
firstrun = 1;

outputfilehandle = fopen(outputfile,'a'); 

for filenumber = 1:length(h5files)
    fn = h5files(filenumber).name;
    ['using: ' fn]
    %signal = h5read(fullfile(directory,fn),'/FullSpectra/SumSpectrum');
    timebins = getnumberofinstancesinh5(fullfile(directory,fn), 'timebins');
    writes = getnumberofinstancesinh5(fullfile(directory,fn), 'writes');
    buffers = getnumberofinstancesinh5(fullfile(directory,fn), 'buffers');
    fprintf(outputfilehandle, '%s\n', fn);
    
    for buffer = 1:buffers
       for write = 1:writes
           signal = h5read(fullfile(directory,fn),'/FullSpectra/TofData', [1 1 buffer write], [timebins 1 1 1]);
           if firstrun == 1
               signalsum = zeros(size(signal));
               firstrun = 0;
           end
           signalsum = signalsum + signal;
           h5write(newh5, '/FullSpectra/TofData', signal, [1 1 newbuffer 1], [timebins 1 1 1]);
           newbuffer = newbuffer + 1;
       end
    end
end
timebinsinnewh5 = getnumberofinstancesinh5(newh5, 'timebins');
timebindifference = timebinsinnewh5 - timebins;
signalsum(end+1:end+timebindifference) = 0;
h5write(newh5, '/FullSpectra/SumSpectrum', signalsum);
fclose(outputfilehandle);