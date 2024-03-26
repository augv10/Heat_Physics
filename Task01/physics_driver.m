% Read the land_sea mask


era5_file='G:/ERA5_PHYSICS/geopot_land_sea_mask.grib';

info = georasterinfo(era5_file);
metadata = info.Metadata;
band_lsm=find(info.Metadata.Element=='LSM');
band_geopot=find(info.Metadata.Element=='Z');


[LSM,R] = readgeoraster(era5_file,Bands=band_lsm);
[ORO,R] = readgeoraster(era5_file,Bands=band_geopot);
ORO=ORO/9.81;


xax_era5=R.LongitudeLimits(1):R.SampleSpacingInLongitude:R.LongitudeLimits(end);
yax_era5=R.LatitudeLimits(end):-R.SampleSpacingInLatitude:R.LatitudeLimits(1);

n_filter=numel(LSM(LSM>=0.5));
xy_filter=zeros(n_filter,2);
xy=zeros(numel(yax_era5)*numel(xax_era5),2);
ij=zeros(numel(yax_era5)*numel(xax_era5),2);
ij_filter=zeros(n_filter,2);
ig=0;ig1=0;
for j=1:numel(yax_era5)
    for i=1:numel(xax_era5)
        ig=ig+1;
        xy(ig,1)=yax_era5(j);
        xy(ig,2)=xax_era5(i);
        ij(ig,1)=j;
        ij(ig,2)=i;
        if LSM(j,i)>=0.5
            ig1=ig1+1;
            xy_filter(ig1,1)=yax_era5(j);
            xy_filter(ig1,2)=xax_era5(i);
            ij_filter(ig1,1)=j;
            ij_filter(ig1,2)=i;
        end
    end
end




T2M1D=[];
D2M1D=[];
iy=2006;
for im=1:12
clear cal
cal=datetime(iy,im,1,0,0,0):hours(1):datetime(iy,im,eomday(iy,im),23,0,0);

% First read T2M D2M

era5_file=['D:/ERA5/',num2str(iy),'/era5_',num2str(iy),sprintf('%.2d',im),'.grib'];

info = georasterinfo(era5_file);
metadata = info.Metadata;
band_t2m=find(info.Metadata.Element=='2T');
band_d2m=find(info.Metadata.Element=='2D');
band_sp=find(info.Metadata.Element=='SP');


[T2M,R] = readgeoraster(era5_file,Bands=band_t2m);
%[D2M,R] = readgeoraster(era5_file,Bands=band_d2m);

% filter
tmp_1=zeros(n_filter,numel(cal));
tmp_2=zeros(n_filter,numel(cal));

for ig1=1:n_filter
    tmp_1(ig1,:)=T2M(ij_filter(ig1,1),ij_filter(ig1,2),:)-273.15;
   % tmp_2(ig1,:)=D2M(ij_filter(ig1,1),ij_filter(ig1,2),:)-273.15;
end
clear T2M D2M
T2M1D=[T2M1D ; reshape(tmp_1,n_filter*numel(cal),1)]; 
end %im


SSRD1D=[];
iy=2006;
for im=1:12
clear cal
cal=datetime(iy,im,1,0,0,0):hours(1):datetime(iy,im,eomday(iy,im),23,0,0);
era5_file=['G:/ERA5_PHYSICS/',num2str(iy),'/era5_',num2str(iy),sprintf('%.2d',im),'.grib'];



info = georasterinfo(era5_file);
metadata = info.Metadata;
band_slhf=find(info.Metadata.Element=='SLHF'); % Surface latent heat flux
band_ssr=find(info.Metadata.Element=='SSR'); % Surface solar radiation
band_str=find(info.Metadata.Element=='STR'); % Surface thermal radiation
band_sshf=find(info.Metadata.Element=='SSHF'); % Surface sensible heat flux
band_ssrd=find(info.Metadata.Element=='SSRD'); % Surface solar radiation downwards


[SSRD,R] = readgeoraster(era5_file,Bands=band_ssrd);
tmp_1=zeros(n_filter,numel(cal));
for ig1=1:n_filter
    tmp_1(ig1,:)=SSRD(ij_filter(ig1,1),ij_filter(ig1,2),:)/3600;
end
clear SSRD
SSRD1D=[SSRD1D ; reshape(tmp_1,n_filter*numel(cal),1)]; 
end %im


cut_off_t2m=quantile(T2M1D(SSRD1D>20),0.95);
clear oo
oo=find(SSRD1D>20 & T2M1D>cut_off_t2m);
TRGT_DAY=zeros(numel(SSRD1D),1);
TRGT_DAY(oo)=1;

clear qq
qq=reshape(TRGT_DAY,n_filter,yeardays(iy)*24);
for ig1=1:n_filter
    for t=1:yeardays(iy)*24
    if qq(ig1,t)==1
    plot(xy_filter(:,2),xy_filter(:,1),'.')
    end
    end
end


for ig1=1:numel(SSRD1D)
    if SSRD1D(ig1) > 




[SLHF,R] = readgeoraster(era5_file,Bands=band_slhf);
SLHF1D=zeros(n_filter,numel(cal));

for ig1=1:n_filter
    SLHF1D(ig1,:)=SLHF(ij_filter(ig1,1),ij_filter(ig1,2),:)/3600;
end
clear SLHF

[SSR,R] = readgeoraster(era5_file,Bands=band_ssr);
SSR1D=zeros(n_filter,numel(cal));

for ig1=1:n_filter
    SSR1D(ig1,:)=SSR(ij_filter(ig1,1),ij_filter(ig1,2),:)/3600;
end
clear SSR

[STR,R] = readgeoraster(era5_file,Bands=band_str);
STR1D=zeros(n_filter,numel(cal));

for ig1=1:n_filter
    STR1D(ig1,:)=STR(ij_filter(ig1,1),ij_filter(ig1,2),:)/3600;
end
clear STR

[SSHF,R] = readgeoraster(era5_file,Bands=band_sshf);
SSHF1D=zeros(n_filter,numel(cal));
for ig1=1:n_filter
    SSHF1D(ig1,:)=SSHF(ij_filter(ig1,1),ij_filter(ig1,2),:)/3600;
end
clear SSHF


T=table(SLHF1D',SSR1D',STR1D',SSHF1D',SSRD1D','VariableNames',{'SLHF','SSR','STR','SSHF','SSRD'});
T.Properties.VariableDescriptions={'Surface latent heat flux','Surface solar radiation', ...
'Surface thermal radiation','Surface sensible heat flux','Surface solar radiation downwards'};

t=1; % time
v=5; % variable

V=table2array(T(t,v));
