clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mass = attograms = 10^-18 gram
% distance = nanometers
% time = nanoseconds
% energy = attogram-nanometer^2/nanosecond^2
% velocity = nanometers/nanosecond
% force = attogram-nanometer/nanosecond^2
% torque = attogram-nanometer^2/nanosecond^2
% temperature = Kelvin
% pressure = attogram/(nanometer-nanosecond^2)
% dynamic viscosity = attogram/(nanometer-nanosecond)
% charge = multiple of electron charge (1.0 is a proton)
% dipole = charge-nanometer
% electric field = volt/nanometer
% density = attograms/nanometer^dim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters
BeadSize = 1; 
WaterEta=1; %ag/nm/ns %0.001 kg/m/s
BeadCsi = 6*pi*WaterEta*BeadSize/2; %ag/ns
% temperature
Temp = 280; 
index = 1; Rep_num = 1;
kB=1.38*10^-2; %ag*nm^2/ns^2/K 1.38*10^-23 Kg*m^2/s^2/K
kBT=kB*Temp;
% diffusivity
Damp=1; % velocity relaxation time (~gamma)
D=kBT/BeadCsi; BeadMass=Damp*BeadCsi;
% filepath and save
SaveFolder='Parameters/'; mkdir(SaveFolder);
% save([SaveFolder 'Parameter.mat'],'BeadSize','BeadCsi','Damp','Temp','kBT');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Make the initial condition file and write .in file %%%%%%%
% create directory to store initial systems and set parameters
BoxSize = [250, 30, 30]; 
L1 = 24; % polymer lengths
NP1 = 625; % number of polymers
NP = 625*L1; % number of monomers
A = 4; % prefactor in soft potential
Rep = 1;
tempFILE = ['L_' num2str(L1) ...
            '_N_' num2str(NP1) '_A' ...
            num2str(round(L1/2)) 'B' num2str(round(L1/2)) '_T' num2str(Temp)];
ini_cond(L1,NP,SaveFolder,BoxSize,tempFILE, mode)
% output location of the initial file
ReadFolder = 'InitialFiles/'; mkdir(ReadFolder);
ReadFilename = [tempFILE '.initial'];
% output file for dump
OutFolder = 'Outputs/'; mkdir(OutFolder);
makefile_output(SaveFolder, ReadFolder,tempFILE, BoxSize, BeadMass);
% output location of the in.file
InFolder = 'InFiles/'; mkdir(InFolder);
InFilename = ['Runners' num2str(index) '.in'];
OutFilename = [tempFILE '_Rep1' num2str(Rep)];
fid = fopen('Index_ID.txt','a');
fprintf(fid,['Index ' num2str(index) ': ' tempFILE '_Rep' num2str(Rep) '\n']);
fclose(fid);
write_IN(InFolder, InFilename, ReadFolder, ReadFilename, OutFolder, OutFilename, BeadSize, kBT, Temp, Damp, A);
index = index + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% subfunctions %%%%%%%
% initial setting generator
function ini_cond(L1, NP, Folder, BoxSize, tempFILE, mode)
% Ratio: monomer ratio, NP: number of monomers, Folder: initial filepath
    NP1= round(NP/L1); % number of polymers type 1
    % allocate polymer stack dimensions - need to expand if NP1+NP2>2*NB^2
    NB = 25; % number of polymers to stack in y and z
    x = [-L1/2 L1/2]; y=(1:NB)-mean(1:NB); z=(1:NB)-mean(1:NB);    
    [X,Y,Z]=meshgrid(x,y,z);
    X=reshape(X,[],1); Y=reshape(Y,[],1); Z=reshape(Z,[],1);
    nx=length(x); ny=length(y); nz=length(z); % number polymers in each dim
    Polymer1=zeros(3,L1,NP1); % (xyz, monomers in polymer, polymer number)
    for np=1:nx*ny*nz % stack polymers - align along x-dir
        Polymer1(1,:,np)=(0:L1-1)+X(np,1)-L1/2;
        Polymer1(2,:,np)= Y(np,1); Polymer1(3,:,np)= Z(np,1);
    end
    figure(); hold on
    for np=1:NP1
        plot3(Polymer1(1,:,np),Polymer1(2,:,np),Polymer1(3,:,np),'k.-');
    end
    axis equal; grid on
    axis([-BoxSize(1)/2 BoxSize(1)/2 -BoxSize(2)/2 BoxSize(2)/2 -BoxSize(3)/2 BoxSize(3)/2])
    % Assign bonds
    NM = NP; % number of monomers - recalculated since rounded above
    Monomer = zeros(3,NM); % pre-allocate monomer position 
    NB = (L1-1)*NP1; % number of bonds (intra-chain only)
    Bond = zeros(1,NB); % preallocate bonds
    Atype = zeros(1,NM); % preallocate bond type
    Btype = zeros(1,NB);
    nm = 0; nb = 0;
    atype = [ones(L1/2,1); 2*ones(L1/2,1)]; % atom type in polymer
    for np=1:NP1 % run through polymer
        for l=1:L1 % run through each monomer
            nm=nm+1;
            Monomer(:,nm)=Polymer1(:,l,np); % monomer position
            Atype(1,nm) = atype(l); % bond type 1 
            if l~=L1 % non-terminal monomer
                nb=nb+1;
                Bond(1,nb)=nm; Bond(2,nb)=nm+1; % intra-chain bond to adjacent monomer
                Btype(1,nb) = 1;
            end
        end
    end
    % save initial data
    save([Folder tempFILE '.mat'],'Monomer','Bond','Atype','Btype', 'BoxSize');
    % output: .mat files of Atype, Bond, BoxSize and Monomer position
end
% make initial file 
function makefile_output(SaveFolder, ReadFolder,tempFILE, BoxSize, BeadMass)
    load([SaveFolder tempFILE '.mat']);
    Filename = [tempFILE '.initial'];
    % .initial file
    Natom_type=2; %number of atom types;
    Natom=length(Monomer);
    Nbond_type=1; %number of bond types;
    Nbond=length(Bond); %number of DNA bonds;
    Nangl_type=0;  %number of angle types;
    Nangl=0; %number of angles;
    Ndihe_type=0; %number of dihedral types;
    Ndihe=0; %number of dihedrals;
    Nimpr_type=0; %number of improper types;
    Nimpr=0; %number of impropers;

    xlo=-BoxSize(1)/2; xhi=BoxSize(1)/2; %x boundary
    ylo=-BoxSize(2)/2; yhi=BoxSize(2)/2; %y boundary
    zlo=-BoxSize(3)/2; zhi=BoxSize(3)/2; %z boundary

    V=zeros(3,Natom);

    fid=fopen([ReadFolder Filename],'w');
    fprintf(fid,'LAMMPS chain data file\n\n');
    fprintf(fid,'%d atoms\n', Natom);
    fprintf(fid,'%d bonds\n', Nbond);
    fprintf(fid,'%d angles\n', Nangl);
    fprintf(fid,'%d dihedrals\n', Ndihe);
    fprintf(fid,'%d impropers\n\n', Nimpr);
    fprintf(fid,'%d atom types\n', Natom_type);
    fprintf(fid,'%d bond types\n', Nbond_type);
    fprintf(fid,'%d angle types\n', Nangl_type);
    fprintf(fid,'%d dihedral types\n', Ndihe_type);
    fprintf(fid,'%d improper types\n\n', Nimpr_type);
    fprintf(fid,'%8.5f %8.5f xlo xhi\n', xlo, xhi);
    fprintf(fid,'%8.5f %8.5f ylo yhi\n', ylo, yhi);
    fprintf(fid,'%8.5f %8.5f zlo zhi\n\n', zlo, zhi);

    fprintf(fid,'Masses\n\n');
    fprintf(fid,'%d %8.5f\n',1,BeadMass);
    fprintf(fid,'%d %8.5f\n\n',2,BeadMass);

    fprintf(fid,'Atoms\n\n');
    for i=1:Natom
        Mole_type=1;
        Atom_type=Atype(i);
        fprintf(fid,[num2str(i) ' ' num2str(Mole_type) ' ' num2str(Atom_type) ' ' ...
                     num2str(Monomer(1,i)) ' ' num2str(Monomer(2,i)) ' ' num2str(Monomer(3,i)) '\n']);
    end
    fprintf(fid,'\n');

    fprintf(fid,'Velocities\n\n');
    for i=1:Natom
        fprintf(fid,[num2str(i) ' ' num2str(V(1,i)) ' ' num2str(V(2,i)) ' ' num2str(V(3,i)) '\n']);
    end
    fprintf(fid,'\n');

    fprintf(fid,'Bonds\n\n');
    for i=1:Nbond
        fprintf(fid,'%d %d %d %d\n',...
                i,Btype(i),Bond(1,i),Bond(2,i));
    end
    fprintf(fid,'\n');
    fclose(fid); 
end
% write the LAMMPS .in file
function write_IN(InFolder, InFilename, ReadFolder, ReadFilename, OutFolder, OutFilename, BeadSize, kBT, Temp, Damp, A)
    kBT_O = kBT;    
    kBT = kBT*300/Temp; % set to 300
    Ai = -0*kBT; Af = -A*kBT; 
    K = 0.56*kBT_O; % FENE scales with kBT
    Rb = 40; R0 = 5;
    dt = Damp/100;
    Tsim = 1e8;
    % write in.file
    fid = fopen([InFolder InFilename],'w');
    fprintf(fid, 'units nano \n');
    fprintf(fid, 'boundary p p p \n');
    fprintf(fid, 'atom_style bond\n\n');
    fprintf(fid, 'processors 1 * *\n\n');
    % load restart file
    fprintf(fid, ['read_data ' ReadFolder ReadFilename  '\n\n']);
    % pair style
    fprintf(fid, ['pair_style hybrid lj/cut ' num2str(BeadSize*2^(1/6)) ' soft ' num2str(BeadSize/2) '\n']);
    fprintf(fid, ['pair_coeff 1 1 lj/cut ' num2str(kBT) ' ' num2str(BeadSize) ' ' num2str(BeadSize*2^(1/6)) '\n']);
    fprintf(fid, ['pair_coeff 2 2 lj/cut ' num2str(kBT) ' ' num2str(BeadSize) ' ' num2str(BeadSize*2^(1/6)) '\n']);
    fprintf(fid, ['pair_coeff 1 2 soft ' num2str(Ai) ' ' num2str(BeadSize/2) '\n']);
    fprintf(fid, 'pair_modify shift yes\n');
    fprintf(fid, 'special_bonds lj/coul 1.0 1.0 1.0\n\n');
    % bond style
    fprintf(fid, 'bond_style fene/expand\n');
    fprintf(fid, ['bond_coeff * ' num2str(K) ' ' num2str(R0) ' ' num2str(0) ' ' num2str(0) ' ' num2str(0) '\n\n']);
    fprintf(fid, ['variable A equal "ramp(' num2str(Ai) ',' num2str(Af) ')"\n\n']);
    fprintf(fid, ['region wallx block -' num2str(Rb) ' ' num2str(Rb) ...
                                ' -' num2str(100) ' ' num2str(100) ...
                                ' -' num2str(100) ' ' num2str(100) ' open 3 open 4 open 5 open 6\n\n']);
    % sim params
    fprintf(fid, ['neighbor ' num2str(R0+BeadSize-BeadSize*2^(1/6)) ' bin \n']);
    fprintf(fid, 'neigh_modify every 1 delay 0\n\n');
    % fixes
    fprintf(fid, 'fix 1 all nve\n');
    fprintf(fid, ['fix 2 all langevin ' num2str(Temp) ' ' num2str(Temp) ' ' num2str(Damp) ' ' num2str(randi(10^7)) ' zero no\n']);
    fprintf(fid, 'fix 3 all adapt 1 pair soft a 1 2 v_A\n'); % anneal
    fprintf(fid, ['fix 4 all wall/region wallx lj126 ' num2str(kBT) ' ' num2str(BeadSize) ' ' num2str(BeadSize*2^(1/6)) '\n\n']);
    % simulation
    fprintf(fid, ['thermo ' num2str(Tsim) '\n']);
    fprintf(fid, ['timestep ' num2str(dt) '\n\n']);
    fprintf(fid, ['run ' num2str(Tsim/10) '\n']);
    fprintf(fid, 'unfix 3\n');
    fprintf(fid, ['run ' num2str(Tsim/10) '\n']);
    fprintf(fid, 'unfix 4\n');
    fprintf(fid, ['run ' num2str(Tsim) '\n\n']);
    % dumps 
    fprintf(fid, ['dump 1 all movie ' num2str(round(Tsim/100)) ' ' OutFolder OutFilename '.mpeg type type zoom 4.5 box yes 0.01 view 85 85 size 1000 400 shiny 0.5\n']);
    fprintf(fid, 'dump_modify 1 acolor 1 blue\n'); 
    fprintf(fid, 'dump_modify 1 acolor 2 orange\n'); 
    fprintf(fid, 'dump_modify 1 adiam 1 3\n'); 
    fprintf(fid, 'dump_modify 1 adiam 2 3\n\n');
    fprintf(fid, ['dump 2 all custom ' num2str(round(Tsim/400))  ' ' [OutFolder OutFilename '_pos'] '.xyz' ' id' ' x y z\n']); 
    fprintf(fid, ['run ' num2str(Tsim) '\n']);
    fprintf(fid, ['write_restart ' OutFolder OutFilename '.restart\n']);
    fclose(fid);
end






















