% takes a wavelengths file and averages over n lines - returning a file
% with 1/n times the number of lines, but two columns - the second one
% being the error

wlfile = 'Z:\Experiments\Clustof\C60He + Spektroskopie\ID 1\Abschnitt 1\Bereich 3 Abschnitt 1 fein 960nm - Wellenlaengen.txt';
newfile = 'Z:\Experiments\Clustof\C60He + Spektroskopie\ID 1\Abschnitt 1\Bereich 3 Abschnitt 1 fein 960nm - Wellenlaengen averaged.txt';

% average over x lines
average = 10;

data = dlmread(wlfile);
oldlines = size(data, 1);

datanew = reshape(data, [average, oldlines/average])';
wl = mean(datanew,2);
err = std(datanew,0,2);

dlmwrite(newfile, horzcat(wl, err),'delimiter','\t','precision',8);