clear; clc

L1 = 24; 
N1 = 625;
mode = 5;
RepNum = 1; A = 4;
T = 280; 

currFILE = generate_filenames(L1, N1, T, 1);
for jj = 1:length(currFILE)
    tempFILE = currFILE{jj};
    % output location of the in.file 
    InFolder = '/Restart_Files/'; mkdir(InFolder)
    % output location of the restart.file
    ReadFolder = '/Outputs/';
    ReadFilename = [tempFILE '1.restart'];
    % output file for dump
    OutFolder = '/Restart_Outputs/'; mkdir(OutFolder)
    for rc = 1:RepNum
        InFilename = ['Runners' num2str(crc) '.in'];
        OutFilename = [tempFILE num2str(rc)];
        write_IN(InFolder, InFilename, ReadFolder, ReadFilename, OutFolder, OutFilename, A, Temp(kk));
        crc = crc + 1;
    end
end

function currFILE = generate_filenames(L1, N1, Temp, mode)
    L = length(N1);
    currFILE = cell(L,1);
    cc = 1;
    for jj = 1:L
        switch mode
            case 1
                tempFILE = ['L_' num2str(L1) ...
                    '_N_' num2str(N1) '_A' ...
                    num2str(round(L1/2)) 'B' num2str(round(L1/2)) '_T' num2str(Temp) '_Rep'];
        end
        currFILE{cc} = tempFILE;
        cc = cc + 1;
    end
end
function write_IN(InFolder, InFilename, ReadFolder, ReadFilename, OutFolder, OutFilename, A, Temp)
    kBT_O = (1.38*10^-2)*Temp;
    kBT = (1.38*10^-2)*300; % set to 300
    BeadSize = 1;
    Af = -A*kBT; 
    R0 = 5; Damp = 1;
    dt = Damp/100;
    K = 0.56*kBT_O;
    Tsim = 1e8;
    % write in.file
    fid = fopen([InFolder InFilename],'w');
    fprintf(fid, 'units nano \n');
    fprintf(fid, 'boundary p p p \n');
    fprintf(fid, 'atom_style bond\n\n');
    fprintf(fid, 'processors 1 * *\n\n');
    % load restart file
    fprintf(fid, ['read_restart ' ReadFolder ReadFilename '\n\n']);
    % pair style
    fprintf(fid, ['pair_style hybrid lj/cut ' num2str(BeadSize*2^(1/6)) ' soft ' num2str(BeadSize/2) '\n']);
    fprintf(fid, ['pair_coeff 1 1 lj/cut ' num2str(kBT) ' ' num2str(BeadSize) ' ' num2str(BeadSize*2^(1/6)) '\n']);
    fprintf(fid, ['pair_coeff 2 2 lj/cut ' num2str(kBT) ' ' num2str(BeadSize) ' ' num2str(BeadSize*2^(1/6)) '\n']);
    fprintf(fid, ['pair_coeff 1 2 soft ' num2str(Af) ' ' num2str(BeadSize/2) '\n']);
    fprintf(fid, 'pair_modify shift yes\n');
    fprintf(fid, 'special_bonds lj/coul 1.0 1.0 1.0\n\n');
    % bond style
    fprintf(fid, 'bond_style fene/expand\n');    
    fprintf(fid, ['bond_coeff * ' num2str(K) ' ' num2str(R0) ' ' num2str(0) ' ' num2str(0) ' ' num2str(0) '\n\n']);
    % sim params
    fprintf(fid, ['neighbor ' num2str(R0+BeadSize-BeadSize*2^(1/6)) ' bin \n']);
    fprintf(fid, 'neigh_modify every 1 delay 0\n\n');
    % fixes
    fprintf(fid, 'fix 1 all nve\n');
    fprintf(fid, ['fix 2 all langevin ' num2str(Temp) ' ' num2str(Temp) ' ' num2str(Damp) ' ' num2str(randi(10^7)) ' zero no\n']);
    % simulation
    fprintf(fid, ['thermo ' num2str(Tsim) '\n']);
    fprintf(fid, ['timestep ' num2str(dt) '\n\n']);
    % dumps 
    fprintf(fid, ['dump 1 all movie ' num2str(round(Tsim/250)) ' ' OutFolder OutFilename '.mpeg type type zoom 4.5 box yes 0.01 view 85 85 size 1000 400 shiny 0.5\n']);
    fprintf(fid, 'dump_modify 1 acolor 1 blue\n'); 
    fprintf(fid, 'dump_modify 1 acolor 2 orange\n'); 
    fprintf(fid, 'dump_modify 1 adiam 1 3\n'); 
    fprintf(fid, 'dump_modify 1 adiam 2 3\n\n');

    fprintf(fid, ['fix 3 all ave/time 1 ' num2str(round(Tsim/100)) ' ' num2str(Tsim/100) ' c_thermo_press[1] c_thermo_press[2] c_thermo_press[3] file ' OutFolder OutFilename '_GloPress' '.dump\n\n']);
    fprintf(fid, ['dump 2 all custom ' num2str(round(Tsim/100))  ' ' [OutFolder OutFilename '_pos'] '.xyz' ' id' ' x y z\n']); 
    fprintf(fid, ['run ' num2str(Tsim) '\n']);
    fclose(fid); 
end