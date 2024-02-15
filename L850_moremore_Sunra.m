clear
close all
clc

frequency = 300;
nn = 9;
dep = 14;
rs = 50;
roi = 50;
nroi = 50;
if nn==3
    roi = 20;
    nroi = 20;
end
path = ['/data/eggebrecht/data1/Weihao/CSF/semi_homo/new/voxel_1mm/',num2str(nn),'/',num2str(frequency)];
% path0 = ['/data/eggebrecht/data1/Weihao/CSF/semi_homo/Beta/',num2str(nn)];
load([path,'/A_',num2str(nn),'_',num2str(frequency),'.mat'])
load('/data/eggebrecht/data1/matlab_codes/neurodot_dev/NeuroDOT_Beta/Support_Files/Spectroscopy/Ecoeff_Prahl.mat','prahlEC')
keep1 = (info.pairs.WL==1 & info.pairs.r3d<=rs);
keep2 = (info.pairs.WL==2 & info.pairs.r3d<=rs);
Extco(1,:) = prahlEC(ismember(prahlEC(:,1),690),2:3).*[2.303,2.303]./10000000;
Extco(2,:) = prahlEC(ismember(prahlEC(:,1),850),2:3).*[2.303,2.303]./10000000;

l1 = [0.1,0.06,0.01,0.008,0.006,0.003,0.001,0.0006];
l2 = [0.6,0.3,0.1,0.06,0.03,0.01,0.006,0.003,0.001,0.0006,0.0003,0.0001,0.00006,0.00003];
psfside = 10;
psfnum = psfside^2;
ale = zeros(psfnum,length(l1),length(l2));
afw = zeros(psfnum,length(l1),length(l2));
afv = zeros(psfnum,length(l1),length(l2));
aer = zeros(psfnum,length(l1),length(l2));
asnr = zeros(psfnum,length(l1),length(l2));
asignal = zeros(psfnum,length(l1),length(l2));
anoise = zeros(psfnum,length(l1),length(l2));
anorm = zeros(psfnum,length(l1),length(l2));
afov = zeros(6,length(l1),length(l2));
a690 = squeeze(A(keep1,:));
a850 = squeeze(A(keep2,:));
if frequency>0
    a690 = [real(a690); imag(a690)];
    a850 = [real(a850); imag(a850)];
end
ax = ceil(size(t1,1)/2);
ay = ceil(size(t1,2)/2);
tic
for i = 4:5
    for j = 1:length(l2)
%         iA690=Tikhonov_invert_Amat_test(a690,l1(i),l2(j));
        iA850=Tikhonov_invert_Amat_test(a850,l1(i),l2(j));
%         iA690=smooth_Amat(iA690,info.tissue.dim,3);
        iA850=smooth_Amat(iA850,info.tissue.dim,3);
        ffr=makeFlatFieldRecon(a850,iA850);

        fooV=Good_Vox2vol(ffr,info.tissue.dim);
        fooV=fooV./max(fooV(:));
        foV = fooV(ax-roi:ax+roi,ay-roi:ay+roi,:);
        afov(1,i,j) = sum(foV(:)>0.5);
        afov(2,i,j) = sum(foV(:)>0.4);
        afov(3,i,j) = sum(foV(:)>0.3);
        afov(4,i,j) = sum(foV(:)>0.2);
        afov(5,i,j) = sum(foV(:)>0.1);
        afov(6,i,j) = sum(foV(:)>0.01);
        
        %%% PSF simulation
        for k=1:psfnum
            dim1 = floor(size(t1,1)/2);
            dim2 = floor(size(t1,2)/2);
            dim3 = size(t1,3)-dep;
            if mod(k,psfside)==0
                dim1 = dim1-ceil(psfside/2)+floor(k/psfside);
                dim2 = dim2+floor(psfside/2);
            else
                dim1 = dim1-floor(psfside/2)+floor(k/psfside);
                dim2 = dim2-ceil(psfside/2)+mod(k,psfside);
            end
            perturb_point = sub2ind([info.tissue.dim.nVx info.tissue.dim.nVy info.tissue.dim.nVz], dim1, dim2, dim3);
            Xsim = zeros(size(a690,2),1);
            Xsim(ismember(info.tissue.dim.Good_Vox, perturb_point),1) =1;
            rsd690 = squeeze(info.pairs.r3d(keep1));
            rsd850 = squeeze(info.pairs.r3d(keep2));
            rand1 = randn(sum(keep1),1);
            rand2 = randn(sum(keep2),1);
            noise_real690 = 0.01*(0.2502*exp(0.02913*rsd690)+4.625*10^(-6)*exp(0.2128*rsd690))*10^(6.769e-4*(frequency-140)).*rand1;
            noise_img690 = (3.933*10^(-11)*exp(0.4161*rsd690)+0.0105*exp(0.05585*rsd690))*10^(0.0013*(frequency-140))*pi.*rand1/180;
            noise_real850 = 0.01*(0.6019*exp(0.01052*rsd850)+9.685*10^(-5)*exp(0.1382*rsd850))*10^(6.785e-4*(frequency-140)).*rand2;
            noise_img850 = (1.917*10^(-10)*exp(0.3708*rsd850)+0.03573*exp(0.02002*rsd850))*10^(0.0013*(frequency-140))*pi.*rand2/180;
            if frequency == 0
                noise690 = noise_real690;
                noise850 = noise_real850;
            else
                noise690 = [noise_real690;noise_img690];
                noise850 = [noise_real850;noise_img850];
            end
            % noise PSF
            Ysim = 5*a850*Xsim+noise850;
            Xrecon = iA850*Ysim;
            % noise Recon          
            N_recon = iA850*noise850;
            % PSF
            psf0 = Good_Vox2vol(Xrecon,info.tissue.dim);
            psf = psf0./max(psf0(:)); 
            % perturbation
            perturb = zeros(size(t1,1),size(t1,2),size(t1,3));
            perturb(dim1,dim2,dim3) = 1;
            
            error = reshape((psf-perturb),size(t1,1)*size(t1,2)*size(t1,3),1);
            anorm(k,i,j) = norm(error);
            psf = StrongestContiguousRegionbeta(psf,0.5);
            noise_recon = Good_Vox2vol(N_recon,info.tissue.dim);         
            imagenoise = noise_recon(ax-nroi:ax+nroi,ay-nroi:ay+nroi,:);
            asignal(k,i,j) = mean(psf0(psf>0));
            anoise(k,i,j) = std(abs(imagenoise(:)));
            asnr(k,i,j) = asignal(k,i,j)/anoise(k,i,j);
            [aFW, aFV, aLE, aER] = Cal_metrics_dev(dim1,dim2,dim3,psf);
            ale(k,i,j) = info.tissue.flags.voxmm*aLE;
            afw(k,i,j) = info.tissue.flags.voxmm*aFW;
            afv(k,i,j) = info.tissue.flags.voxmm*aFV;
            aer(k,i,j) = info.tissue.flags.voxmm*aER;
        end
        i,j
   
    end
end
toc

% psf0 = psf0./max(psf0(:));
psf(dim1,dim2,dim3) = 0.5;
psf0(dim1,dim2,dim3) = -0.1;
pA_psf.slices_type='idx'; pA_psf.slices=[dim1 dim2 dim3];
pA_psf.PD = 0;pA_psf.Scale=max(psf0(:));pA_psf.Th.P=0;pA_psf.Th.N=-pA_psf.Th.P;pA_psf.CH=0;
PlotSlices(t1,info.tissue.dim,pA_psf, psf0)   
pA_psf.PD = 0;pA_psf.Scale=1;pA_psf.Th.P=0;pA_psf.Th.N=-pA_psf.Th.P;pA_psf.CH=0;
PlotSlices(t1,info.tissue.dim,pA_psf, psf)
pA_psf.PD = 0;pA_psf.Scale=max(noise_recon(:));pA_psf.Th.P=0;pA_psf.Th.N=-pA_psf.Th.P;pA_psf.CH=0;
PlotSlices(t1,info.tissue.dim,pA_psf, noise_recon)

fw = squeeze(median(afw,1));
le = squeeze(median(ale,1));
er = squeeze(median(aer,1));
snr = squeeze(median(asnr,1));
signal = squeeze(median(asignal,1));
noise = squeeze(median(anoise,1));
norme = squeeze(median(anorm,1));
fov = afov;

% load([path,'/lambda850_noise_',num2str(nn),'_',num2str(frequency),'MHz.mat'])
% fw1(:,1:2) = squeeze(median(afw(:,:,1:2),1));
% fw1(:,3:10) = fw;
% fw1(:,11:14) = squeeze(median(afw(:,:,3:end),1));
% le1(:,1:2) = squeeze(median(ale(:,1:2),1));
% le1(:,3:10) = le;
% le1(:,11:14) = squeeze(median(ale(:,:,3:end),1));
% er1(:,1:2) = squeeze(median(aer(:,1:2),1));
% er1(:,3:10) = er;
% er1(:,11:14) = squeeze(median(aer(:,:,3:end),1));
% snr1(:,1:2) = squeeze(median(asnr(:,1:2),1));
% snr1(:,3:10) = snr;
% snr1(:,11:14) = squeeze(median(asnr(:,:,3:end),1));
% signal1(:,1:2) = squeeze(median(asignal(:,1:2),1));
% signal1(:,3:10) = signal;
% signal1(:,11:14) = squeeze(median(asignal(:,:,3:end),1));
% noise1(:,1:2) = squeeze(median(anoise(:,1:2),1));
% noise1(:,3:10) = noise;
% noise1(:,11:14) = squeeze(median(anoise(:,:,3:end),1));
% norme1(:,1:2) = squeeze(median(anorme(:,1:2),1));
% norme1(:,3:10) = norme;
% norme1(:,11:14) = squeeze(median(anorme(:,:,3:end),1));
% fov1(:,:,1:2) = squeeze(median(afov(:,:,1:2),1));
% fov1(:,:,3:10) = fov;
% fov1(:,:,11:14) = squeeze(median(afov(:,:,3:end),1));
save([path,'/lambda850_more_noise_',num2str(nn),'_',num2str(frequency),'MHz.mat'],'fw','le','er','snr','signal','noise','norme','fov','-v7.3')
