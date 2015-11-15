% this file reads all *.h5 files in a folder and combines them to a single
% one, where the writes are the sum-spectra of the original files

% new file to use
newh5 = 'new3.h5';

% output file
outputfile = 'filelist_bereich_3_test2.txt';

% directory to combine
directory = 'bereich_3_2_fine_960nm';
h5files = dir([directory '/*.h5']);

newbuffer = 1;
firstrun = 1;
newwrite = 1;

outputfilehandle = fopen(outputfile,'a'); 

timebinsinnewh5 = getnumberofinstancesinh5(newh5, 'timebins');
writesnew = getnumberofinstancesinh5(newh5, 'writes');
buffersnew = getnumberofinstancesinh5(newh5, 'buffers');
signalsum = zeros(timebinsinnewh5, 1);

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
            if timebins >= timebinsinnewh5
                signal = h5read(fullfile(directory,fn),'/FullSpectra/TofData', [1 1 buffer write], [timebinsinnewh5 1 1 1]);
            else
                signal = h5read(fullfile(directory,fn),'/FullSpectra/TofData', [1 1 buffer write], [timebins 1 1 1]);
                timebindifference = timebinsinnewh5 - timebins;
                signal(end+1:end+timebindifference) = 0;
            end
            
            newbuffer
            newwrite
            signalsum = signalsum + signal;
			h5write(newh5, '/FullSpectra/TofData', signal, [1 1 newbuffer newwrite], [timebinsinnewh5 1 1 1]);
			
		    if newbuffer < buffersnew
			    newbuffer = newbuffer + 1;
            elseif newbuffer == buffersnew
				newbuffer = 1;
				newwrite = newwrite + 1;
			end
		end
    end
end

finalsize = size(h5read(newh5, '/FullSpectra/SumSpectrum'), 1);
sumsize = size(signalsum, 1);
if finalsize > sumsize
    diff = finalsize - sumsize;
    signalsum(end+1:end+diff) = 0;
elseif finalsize < sumsize
    signalsum = signalsum(1:finalsize);
end

h5write(newh5, '/FullSpectra/SumSpectrum', signalsum);
fclose(outputfilehandle);