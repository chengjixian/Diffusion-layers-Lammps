clc;close all; clear all;
% MSD slab


bincount=10;

water=192001:384000; %water atoms index 

fid1=fopen('MSD_oxy_bin.txt','wt');
fid3=fopen('MSD_film_avg_oxy.txt','wt');

T_max=4000000;

%getting Tmin from path address automatically
path=cd;


k=strfind(path,'/t_');
x=length('/_');

if path(end)=='/'
path(end)=[];
end

for i=1:length(path)-(k+x)
Tmin_str(i)=path(k+x+i);
end

T_min=str2num(Tmin_str)

cd ..
cd ..

fid=fopen('traj_UW_stage2.lammpstrj','r');


t=0; cnt_steps=0;
while ~feof(fid)
    
    
    
    id=fgetl(fid);
    if strcmp(id,'ITEM: TIMESTEP')==1
        id=fgetl(fid);
        
        Timestep=str2num(id)
        
        if Timestep>T_max 
            break;
        end
        
    end
    
    cnt_if=0;
    if Timestep>=T_min && Timestep<=T_max && cnt_if<1
       
        cnt_if=cnt_if+1;
        t=t+1;
        id=fgetl(fid);
        if strcmp(id,'ITEM: NUMBER OF ATOMS')==1
            
            id=fgetl(fid);
            id=fgetl(fid);
        end
        
        if strcmp(id,'ITEM: BOX BOUNDS pp pp pp')==1 || strcmp(id,'ITEM: BOX BOUNDS pp pp ff')
            cnt_box=10;
            id=fgetl(fid);id=fgetl(fid);id=fgetl(fid);
            if t==1
                
                Z_box=str2num(id);
                bin_size=(Z_box(2)-Z_box(1))/bincount;
            end
            id=fgetl(fid);
        end
        
        if strcmp(id,'ITEM: ATOMS id type xu yu zu ix iy iz ')==1 || strcmp(id,'ITEM: ATOMS id type x y z ix iy iz ')==1
            
            cnt_steps=cnt_steps+1;
            cnt_steps
            %skip lines for all atoms before the first water atom
            for i=1:water(1)-1
                id=fgetl(fid);
            end
            
            water_count=0;
            
            %read the data from a particular timestep
            for i=1:length(water) %loop for the number of water molecules
                
                water_count=water_count+1; % count of water molecules
                id=fgetl(fid);
                ts(t).water_mol(water_count).atom_data=str2num(id);
                
                
                if t==1
                    ts(t).water_mol(water_count).bin_oxy= ceil((ts(t).water_mol(water_count).atom_data(5)-Z_box(1))/bin_size); %bin count for a water molecule from the position of oxygen atom
                end
                
               
            end
            
            % initialize each bin. each bin structure will have all three msd's
            % from all the atoms in the bin
            
            for k=1:bincount
                
                bin(k).MSD_oxy=[0 0 0];
           end
            
            if cnt_steps > 1 %compute msd
                
                %compute msd for all water beads
                for mol=1:length(water)
                    
                    %oxygen_MSD
                    MSD_oxy= (ts(t).water_mol(mol).atom_data(3:5)-ts(1).water_mol(mol).atom_data(3:5)).^2;
                    
                    bin_number=ts(1).water_mol(mol).bin_oxy;
                    
                    bin(bin_number).MSD_oxy=[bin(bin_number).MSD_oxy;MSD_oxy];
                    
                end
                
                
                
                %average binwise
                
                fprintf(fid1,'%i\n',Timestep);
             
                
                film_avg_MSD_oxy=[0 0 0] ;
                
                for k=1:bincount
                    
                    %sum msd in each bin and divide by total moecules in each
                    %bin
                    bin(k).Avg_MSD_oxy= sum(bin(k).MSD_oxy)/(size(bin(k).MSD_oxy,1)-1);

                    
                    %sum over msd's for all molecules in all bins
                    film_avg_MSD_oxy=film_avg_MSD_oxy + sum(bin(k).MSD_oxy);

                    
                    %print msd-x msd-y msd-z and msd-total
                    fprintf(fid1,'%10.3f', [bin(k).Avg_MSD_oxy sum(bin(k).Avg_MSD_oxy)]);
                    fprintf(fid1,'\n');
                    
                end
                
                %divide by total water molecules to get average
                film_avg_MSD_oxy=film_avg_MSD_oxy/(length(water));

                
                fprintf(fid3,'%i\t\t',Timestep);

                
                %print msd-x msd-y msd-z and msd-total
                fprintf(fid3,'%10.3f', [film_avg_MSD_oxy sum(film_avg_MSD_oxy)]);
                fprintf(fid3,'\n');

                
                if t>1
                    t=1; % rewrite the next timestep as t=2
                    ts(2)=[];
                    clear bin;
                end
                
            end
            
        end
        
    end
end
             
fclose(fid);fclose(fid1);fclose(fid3);       
