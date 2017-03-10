filename = 'Curves/ZERO_DET_low_P.csv';

zdet = csvread(filename);

f = transpose(zdet);

csvwrite(filename, f);