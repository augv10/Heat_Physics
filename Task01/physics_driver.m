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

xy=zeros(numel(yax_era5)*numel(xax_era5),2);
ij=zeros(numel(yax_era5)*numel(xax_era5),2);
ig=0;
for j=1:numel(yax_era5)
    for i=1:numel(xax_era5)
        ig=ig+1;
        xy(ig,1)=xax_era5(i);
        xy(ig,2)=yax_era5(j);
        ij(ig,1)=i;
        ij(ig,2)=j;
    end
end

LSM1D=reshape(LSM,numel(yax_era5)*numel(xax_era5),1);
clear oo_filter
oo_filter=find(LSM1D<0.5);


% Testing with parquet

era5_file='G:/ERA5_PHYSICS/2006/era5_200601.grib';

info = georasterinfo(era5_file);
metadata = info.Metadata;
band_slhf=find(info.Metadata.Element=='SLHF'); % Surface latent heat flux
band_ssr=find(info.Metadata.Element=='SSR'); % Surface solar radiation
band_str=find(info.Metadata.Element=='STR'); % Surface thermal radiation
band_sshf=find(info.Metadata.Element=='SSHF'); % Surface sensible heat flux
band_ssrd=find(info.Metadata.Element=='SSRD'); % Surface solar radiation downwards


[SLHF,R] = readgeoraster(era5_file,Bands=band_slhf);
[SSR,R] = readgeoraster(era5_file,Bands=band_ssr);
[STR,R] = readgeoraster(era5_file,Bands=band_str);
[SSHF,R] = readgeoraster(era5_file,Bands=band_sshf);
[SSRD,R] = readgeoraster(era5_file,Bands=band_ssrd);

n_filter=numel(yax_era5)*numel(xax_era5)-numel(oo_filter);
SLHF_1=reshape(SLHF,numel(yax_era5)*numel(xax_era5),size(SLHF,3));
SLHF_1(oo_filter,:)=[];
SLHF1D=reshape(SLHF_1,n_filter*size(SLHF,3),1);


parquetwrite('outagesDefault.parquet',T)
