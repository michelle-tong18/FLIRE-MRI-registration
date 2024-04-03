function [vol_corr M_eddy]= QD_Eddy(vol,qmat,bvals,pe_dim,niter)
%

% Phase encode must be along second dimension
switch pe_dim
  case 1
    vol = permute(vol,[1 2 3 4]);
  case 2
    vol = permute(vol,[2 1 3 4]);
  case 3
    vol = permute(vol,[3 1 2 4]);
end

[nx,ny,ns,nf] = size(vol);
[rowvol,colvol,slicevol] = ndgrid(single(1:nx),single(1:ny),single(1:ns));
coordvol = cat(4,rowvol,colvol,slicevol,ones(size(rowvol)));
clear rowvol colvol slicevol

paramvec = 0;
if ~exist('niter','var'); niter = 5;end

% Compute unit deformation field for each deformation parameter
bmax = max(bvals);
qdiff = zeros(nf,3);
for i = 1:3 % Cardinal axes
    qdiff(:,i) = qmat(:,i)/sqrt(bmax);
end

% PUT IN DEV Version!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta qmat terms
%deltaqdiff = [0 0 0;diff(qdiff,1)];
%qdiff = cat(2,qdiff,deltaqdiff);

% Cross-terms (should check if these are needed)
%for i = 1:3 
%    for j = i:3
%        colnum = colnum+1;
%        qdiff(:,colnum) = (qmat(:,i).*qmat(:,j))/bmax;
%    end
%end

%%%% add drift correction term %%%
%ramp = linspace(0,1,nf)';
%qdiff = cat(2,qdiff,ramp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nparams = size(qdiff,2);

vol_corr = vol; % Start with uncorrected data
cost = 1e10;
for iter = 0:niter
    % synthesis volume based on diffusion tensor
    volDT = QD_Fit_Tensor(vol_corr,qmat,bvals,1);
    vol_err = QD_Synth_Diffusion_Data(volDT.D,qmat,bvals); 

    % calculate error
    vol_err = vol_corr - vol_err;
    newcost = sum(vol_err(:).^2)/sum(vol_corr(:).^2);
    if newcost>cost,break;end % TODO: fix numerical errors in gradient and hessian when porting to c!! 
    cost = newcost;
    fprintf('%s -- %s.m:    iter=%d: cost=%0.6e...\n',datestr(now),mfilename,iter,cost)

    dTdy = (cat(1,zeros(1,ny,ns,nf),diff(vol_corr,1,1))+cat(1,diff(vol_corr,1,1),zeros(1,ny,ns,nf)))/2;
    nbetas = nparams*size(coordvol,4); % total number of parameters to solve for
    g = 0;
    H = 0;
    for k = 1:nf
        colnum = 0;
        tmp0 = zeros(nx,ny,ns,nbetas);
        tmp1 = zeros(nx,ny,ns,nbetas);
        for i = 1:size(coordvol,4)
            for j = 1:nparams
                colnum = colnum+1;
                tmp1(:,:,:,colnum) = squeeze(dTdy(:,:,:,k).*coordvol(:,:,:,i)*qdiff(k,j));
                tmp0(:,:,:,colnum) = squeeze(vol_err(:,:,:,k)).*tmp1(:,:,:,colnum);
            end
        end
        tmp0 = reshape(tmp0,[nx*ny*ns nbetas]);
        tmp1 = reshape(tmp1,[nx*ny*ns nbetas]);
        g = g + sum(tmp0,1)';
        H = H + tmp1'*tmp1;
    end
    clear dTdy vol_err
    betahat = -H\g;
    paramvec = paramvec + betahat;
    dymat = zeros(nx,ny,ns,nf);
    colnum = 0;
    for i = 1:size(coordvol,4) 
        for j = 1:nparams
            colnum = colnum+1;
            for k = 1:nf
                dymat(:,:,:,k) = dymat(:,:,:,k)+coordvol(:,:,:,i)*qdiff(k,j)*paramvec(colnum);
            end
        end
    end
    for k = 1:nf
        vol_corr(:,:,:,k) = QD_interp3(vol(:,:,:,k),dymat(:,:,:,k),zeros(size(dymat(:,:,:,k))),zeros(size(dymat(:,:,:,k))),'linear');
    end
end

% put corrected data back into original space
switch pe_dim
  case 1
    vol_corr = permute(vol_corr,[1 2 3 4]);
  case 2
    vol_corr = permute(vol_corr,[2 1 3 4]);
  case 2
    vol_corr = permute(vol_corr,[3 1 2 4]);
end


% compute equivalent affine matrices (M_eddy) for given paramvec
parammat = [paramvec(1:3) paramvec(4:6) paramvec(7:9) paramvec(10:12)]';
M_eddy = zeros(4,4,nf);
for i = 1:nf
    M_eddy(:,:,i) = eye(4);
    rowtmp = (parammat*qdiff(i,1:3)')';
    M_eddy(1,:,i) = M_eddy(1,:,i)+rowtmp;
end
