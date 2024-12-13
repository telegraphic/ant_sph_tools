function [Q,j,j3D,mmax,nmax] = ReadSphModes(flnm)
% This function reads in spherical mode data in spherical mode Q coefficients
% stored in TIRCA .sph format (See GRASP documentation). 
% The modes are stored in the single index compressed mode index scheme
% (see FEKO documentation). 

% DB Davidson, Feb 2024.

% Status - appears to be working. Storing correct amount of data in
% correct locations, according to the compressed storage scheme.


modedata=importfile_sph(flnm); % First two lines skipped when reading in file
headerlines=6; 
nmax=modedata{1,3};
mmax=modedata{1,4};

mode_counter=1;

lines = headerlines+1;
modecheck=modedata{lines,1};
lines = lines+1;
m = 0;
if modecheck ~= m 
    error(['Mode ',num2str(m),' incorrectly read'])
end
for n=1:nmax
    for s=1:2 % TE then TM modes
        j(mode_counter)  = 2 *(n*(n+1)+ m-1) + s; % compressed mode index
        j3D(j(mode_counter)) = string([num2str(s),' ',num2str(m),' ',num2str(n)]);
        if s == 1
            Q(j(mode_counter)) = complex(modedata{lines,1},modedata{lines,2});
        else
            Q(j(mode_counter)) = complex(modedata{lines,3},modedata{lines,4});
        end
        mode_counter = mode_counter +1;
    end
    lines = lines+1;
end

for mmode = 1:mmax
    modecheck=modedata{lines,1};
    if modecheck ~= mmode
        error(['Mode ',num2str(mmode),' incorrectly read'])
    end
    lines = lines+1; %Advance to next mode.
    for n=abs(mmode):nmax
        for m = -mmode:2*mmode:mmode
            for s=1:2 % TE then TM modes
                j(mode_counter) = 2 *(n*(n+1)+ m-1) + s; % compressed mode index
                j3D(j(mode_counter)) = string([num2str(s),' ',num2str(m),' ',num2str(n)]);
                if s == 1
                    Q(j(mode_counter)) = complex(modedata{lines,1},modedata{lines,2});
                else
                    Q(j(mode_counter)) = complex(modedata{lines,3},modedata{lines,4});
                end
                mode_counter = mode_counter +1;
            end
            lines = lines+1;
        end
    end
end