% Extract position information from .xyz file

Rep = 5;
L1 = 24; N1 = 625;
T_vec = {210:10:250};
fname = {'6A2B2'};
mode = 0;

for sqr = length(T_vec)
    T = T_vec{sqr}; fn = fname{sqr};
    filepath = ['L24/ROUT_', fn, '_run2/'];
    outfilepath = ['L24/XYZFILES_', fn,'_run2/']; mkdir(outfilepath)
    for rr = 1:length(T)
        Temp = T(rr);
        for jj = 1:Rep
            switch mode
                case 0 
                    tempFILE = ['L_' num2str(L1) ...
                        '_N_' num2str(N1) '_6A2B2_T' num2str(Temp) '_Rep' num2str(jj)];
                case 1
                    tempFILE = ['L_' num2str(L1) ...
                        '_N_' num2str(N1) '_A6B6_T' num2str(Temp) '_Rep' num2str(jj)];
                case 2
                    tempFILE = ['L_' num2str(L1) ...
                        '_N_' num2str(N1) '_2A3B3_T' num2str(Temp) '_Rep' num2str(jj)];
                case 3
                    tempFILE = ['L_' num2str(L1) ...
                        '_N_' num2str(N1) '_3A2B2_T' num2str(Temp) '_Rep' num2str(jj)];
                case 4
                    tempFILE = ['L_' num2str(L1) ...
                        '_N_' num2str(N1) '_6AB_T' num2str(Temp) '_Rep' num2str(jj)];
            end
            % import data file
            A = regexp(fileread([filepath tempFILE '_pos.xyz']),'\n','split');
            token = 'ITEM: ATOMS id x y z';
            tokenline = find(contains(A,token)); % find token
            L = length(tokenline); % number of timeslices
            N = L1*N1; % number of atoms
            XYZ = zeros(N,4,L); % preallocate
            for cc = 1:L % toss the values into the matrix
                inds = (tokenline(cc)+1):(tokenline(cc)+N);
                for ii = 1:length(inds)
                    temp = strsplit(cell2mat(A(inds(ii))));
                    for kk = 1:length(temp)-1
                        XYZ(ii,kk,cc) = sscanf(sprintf(' %s',temp{kk}),'%f/%f',[1,Inf]);
                    end
                end
            end
            save([outfilepath tempFILE '_POS.mat'], 'XYZ'); % save the file
        end
    end
end