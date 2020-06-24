startup_mtex
prompt1 = 'How many random orientations are needed: ';
prompt2 = 'How many grains are there?  ';
n = input(prompt1);
grains = input(prompt2);
cs = crystalSymmetry('cubic');
% ss = specimenSymmetry('cubic');
x = n;
while x > 0
    ori = orientation.rand(grains,cs);
    filename = ['cubic_', num2str(grains), '_', num2str(x), '.txt'] 
    fname = fullfile(mtexDataPath, 'ODF', filename);
    export(ori, fname, 'Bunge')
    x = x-1;
end
