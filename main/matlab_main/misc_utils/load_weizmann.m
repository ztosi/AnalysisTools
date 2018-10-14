actions = {'bend', 'jack', 'jump', 'pjump', 'run', 'side', 'skip', ...
    'walk', 'wave1', 'wave2'};
names = {'daria', 'denis', 'eli', 'ido', 'ira', 'lena', 'lyova', 'moshe', ...
    'shahar'};

actlkp = 1:length(actions);
 [bx, by, p, hx, hy, out] = gen_ret_mat(20, 180, 144, 3, 0.01, 0.25);

neups = cell(length(actions)*length(names),1);
l=1;
for i=1:length(actions)
    for j=1:length(names)
        v = VideoReader(['./' actions{i} '/' names{j} '_' actions{i} '.avi']);
        fs = zeros(v.Height * v.Width, v.NumberOfFrames);
        v = VideoReader(['./' actions{i} '/' names{j} '_' actions{i} '.avi']);
        k=1;
        while hasFrame(v)
            gs = double(rgb2gray(v.readFrame))./255;
            fs(:, k) = 1-gs(:);
            k=k+1;
        end
        neups{l} = out*fs;
        l
        l=l+1;
    end
end