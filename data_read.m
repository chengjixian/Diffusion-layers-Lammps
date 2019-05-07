clc;close all;clear all;
Nbins=10;

path='MSD_oxy_bin.txt';
fid=fopen(path,'r');
ts=0;

while ~feof(fid)
for j=1:Nbins+1
    id=fgetl(fid);
    n=str2num(id);
    if j==1
      timestep=n;
      ts=ts+1;
    else
      layer(j-1).data(ts,1:5)=timestep;
      layer(j-1).data(ts,2:end)=n;
    end
end
end

h=figure;
cnt=0;
for j=1:3:Nbins
    cnt=cnt+1;
    X=layer(j).data(:,1);
    X=X-X(1);
    y=layer(j).data(:,2:3);
    Y=sum(y,2);
    
    hold;
    loglog(X,Y,'LineWidth',3);
    s(cnt)={strcat('bin',num2str(j))};
    hold;
    
end
legend(string(s));
xlabel('Time');ylabel('MSD');
savefig(h,'msd_bin.fig');
        