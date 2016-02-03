% this script uses one large h5 file with a lot of buffers and averages
% over a certain amount of spectra (n). the new averaged ones are written
% to buffers in a new file

% new file to use (place one of the templates there)
newh5 = 'Z:\Experiments\Clustof\C60He + Spektroskopie\ID 1\Abschnitt 1\Bereich 3 Abschnitt 1 fein 960nm averaged.h5';

% old file to use
oldh5 = 'Z:\Experiments\Clustof\C60He + Spektroskopie\ID 1\Abschnitt 1\Bereich 3 Abschnitt 1 fein 960nm.h5';

oldbuffers = 2380;
number_to_average = 10;

timebinsinnewh5 = getnumberofinstancesinh5(newh5, 'timebins');
writesnew = getnumberofinstancesinh5(newh5, 'writes');
buffersnew = getnumberofinstancesinh5(newh5, 'buffers');
signalsum = zeros(timebinsinnewh5, 1);
totalsignal = zeros(timebinsinnewh5, 1);

timebinsinoldh5 = getnumberofinstancesinh5(oldh5, 'timebins');

if oldbuffers/number_to_average > buffersnew
    fprint('Not enough buffers in the new file!');
end

j=1;
for i=1:oldbuffers
    if timebinsinoldh5 >= timebinsinnewh5
        signal = h5read(oldh5,'/FullSpectra/TofData', [1 1 i 1], [timebinsinnewh5 1 1 1]);
    else
        signal = h5read(oldh5,'/FullSpectra/TofData', [1 1 i 1], [timebinsinoldh5 1 1 1]);
        timebindifference = timebinsinnewh5 - timebinsinoldh5;
        signal(end+1:end+timebindifference) = 0;
    end
    
    signalsum = signalsum + signal;
    totalsignal = totalsignal + signal;
    
    if mod(i,number_to_average) == 0
        h5write(newh5, '/FullSpectra/TofData', signalsum, [1 1 j 1], [timebinsinnewh5 1 1 1]);
        signalsum = zeros(timebinsinnewh5, 1);
        j=j+1;
        disp(sprintf('Writing buffer: %d',j));
    end
end

finalsize = size(h5read(newh5, '/FullSpectra/SumSpectrum'), 1);
sumsize = size(totalsignal, 1);
if finalsize > sumsize
    diff = finalsize - sumsize;
    totalsignal(end+1:end+diff) = 0;
elseif finalsize < sumsize
    totalsignal = totalsignal(1:finalsize);
end

h5write(newh5, '/FullSpectra/SumSpectrum', totalsignal);