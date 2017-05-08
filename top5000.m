%Copyright (c) 2013, 2014, 2015 Imperial College London
%              2016, 2017       Technical University of Denmark
%All rights reserved.
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%a. Redistributions of source code must retain the above copyright notice,
%this list of conditions and the following disclaimer.
%b. Redistributions in binary form must reproduce the above copyright
%notice, this list of conditions and the following disclaimer in the
%documentation and/or other materials provided with the distribution.
%c. Neither the name of PyFR nor the names of its contributors
%may be used to endorse or promote products derived from this software
%without specific prior written permission.
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
%LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
%OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%DAMAGE.


%TOP5000 by Kristian Ejlebjerg Jensen 2017-01-30
%Department of Micro- and Nanotechnology, Technical University of Denmark (DTU)
%Ørsteds Plads, DK-2800 Kgs. Lyngby
%top5000(nan,5e-3,0.5,false,[],[],0.1,[],1/600,300,1e-3,pi/5,'fig7a',1,1.025,false)
%top5000(4e-3,5e-3,0.5,false,[],[],0.1,[],1/20,424,1e-3,pi/5,'fig7b')
%top5000(2e-3,5e-3,0.5,false,[],[],0.1,[],1/20,300,1e-3,pi/5,'fig7c')
%top5000(1e-3,5e-3,0.5,false,[],[],0.1,[],1/20,300,1e-3,pi/5,'fig7d')
%top5000(0.8,7e-2,0.1,true,0.75,0.7,0.5,1.5,1/10,200,1e-3,pi/4,'fig9a')
%top5000(0.4,7e-2,0.1,true,0.75,0.7,0.5,1.5,1/10,283,1e-3,pi/4,'fig9b')
%top5000(0.2,7e-2,0.1,true,0.75,0.7,0.5,1.5,1/10,400,1e-3,pi/4,'fig9c')
%top5000(4e-3,5e-3,0.5,false,[],[],0.1,[],1/20,200,1e-3,pi/5,'fig10a',0)
%top5000(0.8,7e-2,0.1,true,0.75,0.7,0.5,1.5,1/10,212,1e-3,pi/4,'fig10b',0)
%top5000(nan,7e-2,0.1,true,0.75,0.7,0.5,1.5,1/40,400,1e-3,pi/4,'fig11a',1,1.025,false)
%top5000(0.4,7e-2,0.1,true,0.75,0.7,0.5,1.5,1/10,283,2e-4,pi/4,'fig12b')
%top5000(0.4,7e-2,0.2,true,0.75,0.7,0.5,1.5,1/10,283,1e-3,pi/4,'fig13a')
%top5000(0.2,7e-3,0.1,true,0.75,0.7,0.5,1.5,1/10,400,1e-3,pi/8,'fig14b')

function top5000(eta, Lmin, Tvol, do3D, z1, rw, r1, r2, meshszI, itertotal, kmin, theta, outnm, strg, mvlimit, doadapt)
%%%%% INPUT %%%%%
%eta      : scaling factor for the metric_pnorm, can also be an integer controlling the desired number of nodes
%Lmin     : filter length
%Tvol     : Volume fraction
%do3D     : toggles between the 2D and 3D (true) problem
%z1       : geometric parameter for the 3D problem
%rw       : geometric parameter for the 3D problem
%r1       : geometric parameter 
%r2       : geometric parameter for the 3D problem
%meshszI  : initial mesh size
%itertotal: total number of optimization iterations
%kmin     : ratio of minimum to maximum thermal conductivity
%theta    : geometric parameter 
%outnm    : string specifying the output directory
%strg     : toggle between transfering design variables and sensitivities using a naive hack or real interpolation (true)
%mvlimit  : maximum change of design variables between iterations (optional)
%doadapt  : optinal boolean argument for disabling adaptivity, doadapt==kk triggers adapt every kk iteration
if nargin < 14
 strg = 1;
end;
if nargin < 15
 mvlimit = 1.025;
end;
if nargin < 16
 doadapt = 1;
end;
frzz20 = false;
simpP = 1.;
%initialize mesh, geometry and boundary conditions 
if do3D
  [tri,xy,bndmesh,geomfunc,options] = mesh_cake(r1,r2,z1,rw,round(1./meshszI),theta);
  bcs{1} = 7; bcs2{1} = 12; %options.area = ((1+r2)^2-(1-r1)^2)*z1*theta;
  eta = eta/4.^(1./3.); %Nmetric = Nmetric*4.^(1./3.);
  if isnan(eta)
   Nmetric = metric_uniform(xy, [meshszI meshszI meshszI]); options.innerit = 40;
   [tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
  end;
else
  [tri,xy,bndmesh,geomfunc,options] = mesh_pizza(r1,round(1./meshszI),theta);
  bcs{1} = 3; bcs2 = []; %options.area = (1-r1^2)*theta;
end;
X = fem_getX(tri,xy); [I,dem,demN] = elem_inv(tri,xy);
simpF = ((3+do3D)/simpP)^(2./(itertotal-20*frzz20	));
%export initial design
system(['rm -rf ' outnm ' && mkdir ' outnm]); 
rho = ones(size(demN))*Tvol; flname = [outnm '/rho.pvd']; export_vtk_binary(flname,tri,xy,rho,1,{'rho'}); 
sclr = zeros(itertotal,17);
options.verbose = -1; 
for ii=1:itertotal
 sclr(ii,7) = now();
 rhof = fem_filter(tri,xy,bndmesh,Lmin,rho,bcs2,options);
 rhof = mean(rhof(tri),2); rhof(rhof<0) = 0.;
 sclr(ii,8) = now();
 T = fem_heat(tri,xy,bndmesh,rhof,kmin,simpP,bcs,options,X);
 sclr(ii,1) = sum((kmin+(1.-kmin)*rhof.^simpP).*fem_sq_grad(tri,xy,T,X).*dem)/(6+18*do3D); %obj
 sclr(ii,9) = now();
 dOdrhof = -simpP*(1.-kmin)*rhof.^(simpP-1).*fem_sq_grad(tri,xy,T,X);
 sclr(ii,10) = now();
 dOdrho  = fem_filter(tri,xy,bndmesh,Lmin,dOdrhof,bcs2,options,X); %inconsistent
 %dOdrho = fvm_filter(tri,xy,bndmesh,Lmin,dOdrhof,bcs2,options);
 %dOdrho = full(sparse(tri,ones(size(tri)),repmat(dOdrho.*dem/(2+4*do3D),1,size(tri,2)),size(xy,1),1))./demN;
 sclr(ii,11) = now();
 %adapt
 if ii < itertotal-frzz20*20 && doadapt && mod(ii,doadapt) == 0
 sclr(ii,12) = now();
 Nmetric = metric_pnorm(tri,xy,dOdrho,eta,2,options,X); 
 sclr(ii,13) = now();
 Nmetric = [Nmetric rho dOdrho T];
 if mod(eta,1) == 0
  Nmetric = metric_scale(Nmetric,eta,demN); end;
 if strg
  bndmesh.xyold = xy; end; %triggers exact interpolation
 sclr(ii,14) = now();
 [tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
 sclr(ii,15) = now();
 X = fem_getX(tri,xy); [I,dem,demN] = elem_inv(tri,xy);
 rho = Nmetric(:,end-2); dOdrho = Nmetric(:,end-1); T = Nmetric(:,end); %(potentially) hacky interpolation
 rho(rho < 0.) = 0; rho(1. < rho) = 1.; dOdrho(0. < dOdrho) = 0.; %extrapolation fix
 end;
 %Optimality criteria
 dOdrho = dOdrho.*demN;
 sclr(ii,16) = now();
 rho = optC(rho, dOdrho, demN./sum(demN), Tvol, mvlimit);
 sclr(ii,17) = now();
 %Output
 sclr(ii,6) = simpP;
 sclr(ii,2) = sum(rho.*demN)/sum(demN);
 sclr(ii,3) = sum(sqrt(max(fem_sq_grad(tri,xy,rho,X),zeros(size(tri,1),1))).*dem)/(2.+4.*do3D); %complexity
 sclr(ii,4) = get_ND(tri,rho,dem);
 sclr(ii,5) = numel(rho);
 
 export_vtk_binary(flname,tri,xy,[rho T dOdrho./demN],ii+1,{'rho','T','sensitivity'});
 disp(sprintf('%0d, O = %2.3f, nodes=%d, P=%2.2f, V=%2.2f, C=%2.2f, ND=%0.0f %%',...
 ii,sclr(ii,1),sclr(ii,5),simpP,sclr(ii,2),sclr(ii,3),sclr(ii,4)*100));
 if mod(ii,10) == 0 || ii==itertotal %save every 10 iterations only
  if ii ~= 10
   system(['rm ' outnm '/sclr.mat']);
  end;
  save([outnm '/sclr.mat'],'sclr','-ascii');
 end;
 %update simpP
 if simpP < 3.+do3D
  simpP = simpP*simpF;
  if 3.+do3D < simpP
   simpP = 3.+do3D;
  end;
 end;
end; %end for
%do timings
Ttotal = sum(sclr(:,17)-sclr(:,7))/100; disp(sprintf('solve took %1.1f %%, metric took %1.1f %% adapt took %1.1f %%. The total time was %2.2fh',sum(sclr(:,9)-sclr(:,8))/Ttotal,sum(sclr(:,13)-sclr(:,12))/Ttotal,sum(sclr(:,15)-sclr(:,14))/Ttotal,Ttotal*100*24));
%make STL
if do3D
 [tri,xy] = export_stl(tri,xy,rho,0.5,[outnm '/rho.stl']);
 close all; trimesh(tri,xy(:,1),xy(:,2),xy(:,3),'edgecolor','k');
end;

function xnew = optC(xval, dfdx, dgdx, g, mvlimit)
%see E. Andreassen, et. al., Efficient topology optimization in MATLAB using 88 lines of code,
%    Structural and Multidisciplinary Optimization 43 (1), 1--16 (2011)
%%% INPUT %%%
% xval     : design variables
% dfdx     : objective sensitivity
% dgdx     : constraints sensitivity
% g        : constraint value (g<0)
% mvlimit  : movelimit, limits maximum change in individual elements of xval

%%% OUTPUT %%%
%xnew      : new set of design variables
i1 = 0; i2 = 1e6; nn = numel(dfdx);
while i2 - i1 > 1e-4
  lmid = (i1+i2)/2.;
  xnew = max(max(min(min(xval+mvlimit,xval.*sqrt(abs(-dfdx)./dgdx/lmid)),ones(nn,1)),xval-mvlimit),zeros(nn,1));
  if sum(xnew.*dgdx)/sum(dgdx) - g > 0.
    i1 = lmid;
  else
    i2 = lmid;
  end;
end;
	
function [tri,xy,bndmesh,geomfunc,options] = mesh_pizza(r1,N,theta)
%makes a mesh which is the difference between a disk with radius 1 
%  and a disk with radius r1. Only the first theta radians of the disks are used
%%% INPUT %%%
%r1      : inner radius
%N       : desired inverse edge length 
%theta   : the angle of the slice (theta > pi might give unexcpected results)

%%% OUTPUT %%%
%tri         : Element list, N x 3, where N is the number of elements
%xy          : Node coordinates, M x 2, where M is the number of nodes
%bndmesh.edg : boundare edge list, P x 2, where p is the number of boundary edges.
%              note that issorted(bndmesh.edg,2) == true
%bndmesh.IDs : IDs of boundary edges, P x 1, vertices belonging to edges with different IDs
%              are fixed (corners)
%geomfunc    : signed distance function describing the geometry
%options     : struct for the mesh adaptation. In this case values suitable for 2D are returned
options = gen_options();
options.verbose = false;
options.consRM = false;
[xy,tri,bndmesh,options,geomfunc] = mesh_rect(N,options);
bndmesh = bndmesh_polygon(tri,xy,bndmesh,options);
r = r1+xy(:,1)*(1.-r1); v = xy(:,2)*theta;
xy = [cos(v).*r,sin(v).*r];
%geomfunc
cr1      = @(x_) geom_circ(x_);
cr2      = @(x_) geom_circ(x_, r1);
g1       = @(x_) geom_rect(x_, [1.,1.], [0.,-1]);
g2       = @(x_) geom_rect(x_);
g2r      = @(x_) geom_rotate(x_, g2, theta);
c1       = @(x_) geom_diff(cr1(x_), g1(x_));
c2       = @(x_) geom_diff(c1(x_), g2r(x_));
geomfunc = @(x_) geom_diff(c2(x_),cr2(x_));

%destroy structure
Nmetric = metric_uniform(xy, [1/N 1/N]);
[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
options.consRM = true;

function [tri,xy,bndmesh,geomfunc,options] = mesh_cake(r1,r2,z1,rw,N,theta)
%The function returns an extruded version of mesh_pizza, except there is also a boundary ring on 
%on the bottom where bndmesh.IDs == 7
%%% INPUT %%%
%r1          : inner radius
%r2          : outer radius
%z1          : height
%rw          : inner radius of strip with bndmesh.IDs == 7 (outer radius of strip is 1)
%theta       : the angle of the extruded slice (theta > pi might give unexcpected results)

%%% OUTPUT %%%
%tri         : Element list, N x 4, where N is the number of elements
%xy          : Node coordinates, M x 3, where M is the number of nodes
%bndmesh.fac : boundare face list, P x 3, where p is the number of boundary faces.
%              note that issorted(bndmesh.fac,2) == true
%bndmesh.IDs : IDs of boundary faces, P x 1, edges belonging to faces with different IDs
%              are are to be respected (i.e. cannot be swapped)
%options     : struct for the mesh adaptation. In this case values suitable for 3D are returned

options = gen_options();
options.verbose = false;
options.prag_crs = true;
options.consRM = false;
options.smpRFN = 1;
options.consRFN = 2;
options.swap3D = false;
[xy,tri,bndmesh,options,geomfunc] = mesh_box([3,N,N],options,[]);
I1 = and(0.25 < xy(:,1), xy(:,1) < 0.50);
I2 = and(0.50 < xy(:,1), xy(:,1) < 0.75);
bndmesh.crnds = find(and(and(0.25 < xy(:,1), xy(:,1) < 0.75),xy(:,3) < 1e3*eps));
r = r1+xy(:,1)*(r2-r1); 
r(I1) = rw; r(I2) = 1;
v = xy(:,2)*theta;
z = xy(:,3)*z1;
xy = [r.*cos(v) r.*sin(v) z];
bndmesh = bndmesh_polyhedron(tri,xy,bndmesh,options);
facxy = (xy(bndmesh.fac(:,1),:) + xy(bndmesh.fac(:,2),:) + xy(bndmesh.fac(:,3),:))/3.; 
facr = sqrt(sum(facxy(:,1:2).^2,2)); xyr =  sqrt(sum(xy(:,1:2).^2,2));
I = and(facxy(:,3) < 1e3*eps,and(rw-options.geomtol < facr,facr < 1.+options.geomtol));
bndmesh.IDs(bndmesh.IDs==7) = max(bndmesh.IDs)+1;
bndmesh.IDs(I) = 7; 
Nmetric = metric_uniform(xy, [1/N 1/N 1/N]);
[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
%trimesh(bndmesh.fac,xy(:,1),xy(:,2),xy(:,3))

function ND = get_ND(tri,rho,dem)
rho_ = rho(tri);
if(size(tri,2) == 3)
ND = 4*sum((mean(rho_,2)-sum(rho_(:,[1 2 3 1 2 3]).*rho_(:,[1 2 3 2 3 1])/6,2)).*dem)/sum(dem);
else %3D
ND = 4*sum((mean(rho_,2)-sum(rho_(:,[1 2 3 4 1 2 3 4 1 3]).*rho_(:,[1 2 3 4 2 3 4 1 2 4]),2)/10).*dem)/sum(dem);
end;
	

function [tri,xy] = export_stl(tri,xy,c,cutv,flname)
%%% exports the volume satisfying cutv < c
%%% INPUT %%%
%tri      : element list, N x 4, where N is the number of elements
%xy       : node coordinates, M x 3, where M is the number of elements
%c        : scalar variable,  M x 1, to be cut
%cutv     : lower thresshold 
%flname   : string specifying the output pvd filename
asciiformat = false;
sclr = c';
[fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2,nd2fac,nd2edg] = bks_all3D(tri);
%new xy, edg2xy tri2edg: 12 34 24 13 14 23
cedg  = sum(c(edg) < cutv ,2) == 1;
edg1 = edg(cedg ,:);
v  = sclr(edg1);
edg2xy = zeros(size(edg,1),1); edg2xy(cedg) = size(xy,1)+1:size(xy,1)+sum(cedg);
xy = [xy; xy(edg1(:,1),:).*repmat((cutv - v(:,2))./(v(:,1)-v(:,2)),1,3) + xy(edg1(:,2),:).*repmat((cutv-v( :,1))./(v(:,2)-v(:,1)),1,3)];

%BOUNDARY
bndtri = tri(fac2tri(fac2tri(:,2)==0,1),:);
bndfac = fac(fac2tri(:,2)==0,:); 
bndtrifac = sort([bndtri bndfac]');
nfac = bndtrifac([bndtrifac(1,:)~=bndtrifac(2,:); and(bndtrifac(2:6,:)~=bndtrifac(1:5,:),bndtrifac(2:6,:)~=bndtrifac(3:7,:)); bndtrifac(6,:)~=bndtrifac(7,:)]);
bndtri = bndtri'; I = repmat(nfac',4,1) == bndtri; 
[C,R] = find(I); bndfac_ = reshape(bndtri(not(I)),3,size(bndtri,2));
I = or(C==1,C==3); bndfac_(1:2,I) = bndfac_([2 1],I); bndfac = bndfac_';
%full ones
I1 = cutv < sclr(bndfac);
bndtri = bndfac(sum(I1,2) == 3,:);
%tri ones
It = and(I1,repmat(sum(I1,2)==1,1,3));
[C,R] = find(It');
tri1 = bndfac_(It');
C1 = C+1; C1(C1==4) = 1; C2 = (C1+1); C2(C2==4) = 1;
edg1 = sort([tri1 bndfac(R+(C1-1)*size(bndfac,1))],2);
edg2 = sort([tri1 bndfac(R+(C2-1)*size(bndfac,1))],2);
bndtri = [bndtri; [tri1 edg2xy([bks_bndedg2edg(edg,nd2edg,edg1) bks_bndedg2edg(edg,nd2edg,edg2)])]];

%quad ones
Iq = and(I1,repmat(sum(I1,2)==2,1,3));
[C,R] = find(Iq'); C = reshape(C,2,size(C,1)/2)'; R = R(1:2:end);
I_ = and(C(:,1)==1, C(:,2)==3); C(I_,:) = C(I_,[2 1]);
Iq_ = and(not(I1),repmat(sum(I1,2)==2,1,3));
tri1 = bndfac_(Iq_');
qua1 = bndfac(R+(C(:,1)-1)*size(bndfac,1));
qua2 = bndfac(R+(C(:,2)-1)*size(bndfac,1));
edg1 = sort([tri1 qua1],2);
edg2 = sort([tri1 qua2],2);
qua = [qua1 qua2 edg2xy([bks_bndedg2edg(edg,nd2edg,edg1) bks_bndedg2edg(edg,nd2edg,edg2)])];
sqL1 = sum((xy(qua(:,1),:) - xy(qua(:,4),:).^2),2);
sqL2 = sum((xy(qua(:,2),:) - xy(qua(:,3),:).^2),2);
I1 = sqL1 < sqL2; nI1 = not(I1); %make shortest edge
bndtri = [bndtri; qua(I1,1:3); qua(I1,[2 4 3]); qua(nI1,[1 2 4]); qua(nI1,[1 4 3])];
     
%INTERIOR
tri2cedg_ = cedg(tri2edg)'; tri2edg_ = tri2edg';
Nedg = sum(tri2cedg_);
I1 = Nedg==3; tri2edg_ = tri2edg(I1,:)';
tri_ = reshape(edg2xy(tri2edg_(tri2cedg_(:,I1))),3,sum(I1))';
Ii = sum(cutv < sclr(tri(I1,:)),2) ~= 1;  %i.e. 3
tri_(Ii,1:2) = tri_(Ii,[2 1]);

I2 = Nedg==4; tri2edg_ = tri2edg(I2,:)';
[C,R] = find(cutv < sclr(tri(I2,:))'); C = reshape(C,2,sum(I2));
qua = reshape(edg2xy(tri2edg_(tri2cedg_(:,I2))),4,sum(I2))';
Ii = C(1,:)' ~= 1; %not(or(or(C(1,:)' == 1,C(1,:)' == 3),C(1,:)' == 5));

sqL1 = sum((xy(qua(:,1),:) - xy(qua(:,2),:)).^2,2);
sqL2 = sum((xy(qua(:,3),:) - xy(qua(:,4),:)).^2,2);
I1 = sqL1 < sqL2; %make shortest edge
I_ = and(I1,Ii); qua(I_,1:2) = qua(I_,[2 1]); I_ = and(not(I1),Ii); qua(I_,3:4) = qua(I_,[4 3]);
tri1 = [tri_; qua(I1,[1 2 3]); qua(I1,[2 1 4]); qua(not(I1),[2 3 4]); qua(not(I1),[1 4 3])];
tri = [bndtri(:,[2 1 3]); tri1];

%clean unused nodes
I = false(size(xy,1),1); I(tri(:)) = true;
old2new = zeros(size(xy,1),1); old2new(I) = 1:sum(I);
xy = xy(I,:);
tri = old2new(tri);

if false %0 < size(tri,1) && size(tri,1) < 1e6 
[edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg] = bks_all(tri);
tri1 = tri(edg2tri(:,1),:); tri2 = tri(edg2tri(:,2),:);
[C1,R] = find(repmat(edg(:,1)',3,1) == tri1');
[C2,R] = find(repmat(edg(:,2)',3,1) == tri1');
[D1,R] = find(repmat(edg(:,1)',3,1) == tri2');
[D2,R] = find(repmat(edg(:,2)',3,1) == tri2');
clckwsC = or(or(and(C1==1,C2==2),and(C1==2,C2==3)),and(C1==3,C2==1));
clckwsD = or(or(and(D1==1,D2==2),and(D1==2,D2==3)),and(D1==3,D2==1));
all(clckwsC == not(clckwsD))
sum(not(clckwsC == not(clckwsD)))

if sum(edg2tri(:,2)==0) ~=0
  error('Did not write non-manifold STL');
end;
else
  warning('Did not check validity - file too big');
end;


xy_ = xy'; tri_ = tri';
n = cross(xy(tri(:,2),:)-xy(tri(:,1),:),xy(tri(:,3),:)-xy(tri(:,1),:),2);
trixy = reshape(xy_(:,tri_),9,size(tri,1));
allfloat = [n'; trixy];
if not(asciiformat)
ndn = 'l'; %LittleEndian ('b' for BigEndian)
%stlwrite(flname, tri, xy);
fid = fopen(flname,'w');
stlheader = 'STLFILE HEADER WITH 80 CHARACTERS                             STLFILE HEADER END';
fwrite(fid,sprintf(stlheader),'char',ndn);
fwrite(fid,size(tri,1),'uint32',ndn);
for jj=1:size(tri,1)
fwrite(fid,allfloat(:,jj),'float32',ndn);
fwrite(fid,0,'uint16',ndn);
end;
fclose(fid);
else
fid = fopen(flname,'w');
output1 = 'solid MATLAB_Model\n';
output3 = 'endsolid MATLAB_Model';
output2 = sprintf('  facet normal %0.6f %0.6f %0.6f\n    outer loop\n      vertex %0.6f %0.6f %0.6f\n      vertex %0.6f %0.6f %0.6f\n      vertex %0.6f %0.6f %0.6f\n    endloop\n  endfacet\n',allfloat(:));
output = [output1 output2 output3];
fwrite(fid,sprintf(output),'char');
fclose(fid);
end;

OpenSCADpolyhedron = true;
if OpenSCADpolyhedron 
 flname_ = [flname(1:end-3) 'scad'];
 out1 = 'points = [\n';
 out3 = '];\n\n faces = [\n';
 out5 = '];\n\n polyhedron(points, faces, convexity = 10);';
 out2 = sprintf('[%12.12e, %12.12e, %12.12e],\n',xy'); out2 = out2(1:end-2);
 out4 = sprintf('[%d, %d, %d],\n',tri'-1); out4 = out4(1:end-2);
 output = [out1 out2 out3 out4 out5];
 fid = fopen(flname_,'w');
 fwrite(fid,sprintf(output),'char');
 fclose(fid);
end;
function options = gen_options()
    options.qualM = 1; %mesh quality metric: vasilevski(1), fast orthogonal (2), fast orthogonal vasilevski (3), steiner ellipse (4), steiner ellipse vasilevski (5)
    options.qualP = 0; %combines angles with quality in metric space (~=0), perhaps even with strong weight (>>1).
    options.consRM = 1; %conservative coarsening only acts when local minimum mesh quality is improved (1)
    options.fastRM = 1; %coarsening is performed on edges shorter than options.Lup only (1). In case of node removal, one can choose to prefer removal of nodes with associated with many short edges (fastRM==2)
    options.volRM = 0; %remove nodes associated with elements, whose volumes are smaller than 1 in metric space
    options.RMedg = 0; %use edge collapse rather than node removal for coarsening
    options.qualRM = 0; %collapse shortest edge in metric space rather than make best elements when remeshing in coarsening by node removal (0), deprecated
    options.rmnonconvec = 2; %do not try to convexify non-convex sets during node-removal or edge swapping (3D) (0), allow one concave point (1) or an infinite amount of them (2)
    options.debug  = 0; %test inverted elements (1), wrong area (2 && area~=0), correct internal value of element quality list (2)
    options.area = 0; %total area/volume of mesh, only relevant for non-curved geometries (~=0)
    options.prag_adapt = 2; %only use smoothing as post-processeing (1) or smoothing only (3) or smooth every other inner iteration (4). We can also do only swapping/smoothing (5) or only refinement/coarsening (6). Finally, there is also the possibility of calling each operation untill a convergence criteria is met (0)
    options.outerit = 1; % number of outer iterations in test cases
    options.innerit = 10; %number of inner iterations
    options.consMV = 1; % accept smoothing regardless of whether in improves the minimum local element quality (0)
    options.log_m_add = 0; %add mesh metrics in log space (1)
    options.nosparse = true; %use sparse function (speed increase)
    options.ismatlab = false;
    options.minchg = 1.;
    options.greedyCLR = 2; %greedy colouring algorithm (1), largest ID first (0) or most neighboughs first (fall back to ID in conflicts), greedy colouring with cromatic number reduction (2) and the same with preference to most neigh boughs (3, extremely slow)
    options.verbose = 1; %print mesh statistics (1) and info on operations (2)
    options.geomtol = 1e-12; %tolerance for detecting co-linear lines/faces, when no boundary mesh is supplied
    options.Llow = 1/sqrt(2); %lower threshold for edges in metric space
    options.Lup = sqrt(2); %upper  threshold for edges in metric space
    options.advRFN = 0; %considers four options for refinement of elements with 3 split edges in 2D
    options.smpRFN = 0; %only refine one edge in each element, perhaps using colouring (2)
    %not fully implemented:
    options.minA = 1e-8; %minimum area/volume of elements allowed
    options.min_edgL = 1e-4; %minimum edge length allowed
    options.spltc = 2; %split edges in middle (1) or in middle of metric space (2)
    
    options.swap3D = 1; %do face to edge swapping in 3D
    options.RMnd3D = 1; %prevent splitting of spheres instead of killing them after the fact (1,2) and ditch spheres based on failed edge generation only (2,3)
    options.consRFN = 0; %only refine elements, when it improves all the local mesh quality (1) or also when it only improves the worst local mesh quality (2)
    options.minqual = 1e-9; %minimum quality when options.consRM == 0, options.consRFN == 0 or options.consMV == 0
    options.mntn_bks = true; %maintain active set for nodes (edges (+fac, inverse))
    options.fastRFN = true; %do not refine elements smaller than options.Llow
    options.MVspd = 1.; %smoothing speed, 1 corresponds to laplacian smoothing
    options.MVit = 5; %number of smoothing iterations to use for post-processing
    options.prag_crs = false; %coarsen by collapsing to the shortest (possible) edge
    options.OPTnorm = inf; %norm to optimize in local operations
    options.MVinterp = false; %perform interpolation in smoothing
    options.fem_solve_it = false; %1000;   %maximum number of itrations for solver. 0 triggers direct solver
    options.fem_solve_tol = [];  %relative solver tolerance for iterative solver
    %options.TS = ; if set this specifies the ID of the boundary to be frozen

function [edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg,nd2tri] = bks_all(tri)
edga = [tri(:,2) tri(:,1); tri(:,3) tri(:,2); tri(:,1) tri(:,3)];
edga2tri = [1:size(tri,1) 1:size(tri,1) 1:size(tri,1)]';

[edg,I] = sort(edga,2); %biggest in last column
[edg,Is] = sortrows(edg);
tris = edga2tri(Is);
I = I(Is,1)==1; %true for non flipped
d =[true; or(edg(1:end-1,1)~=edg(2:end,1),...
             edg(1:end-1,2)~=edg(2:end,2))]; Nd = sum(d);
%d =[true; any(edg(1:end-1,:)~=edg(2:end,:),2)]; Nd = sum(d); %slower alternative
ad = not(d); Nad = sum(ad); 
d2 = [ad(2:end); false]; d2 = d2(d);
edg = edg(d,:);

% gen edg2tri
tris2 = zeros(Nd,1); tris2(d2) = tris(ad);
edg2tri = [tris(d) tris2];

%flip edg2tri according to edg permutation
Iflp=and(d2,I(d)); edg2tri(Iflp,:) = edg2tri(Iflp,[2 1]);

%gen edga2edg
edga2edg = zeros(size(edga,1),1);
edga2edg(Is(d))  = 1:size(edg,1);
edga2edg(Is(ad)) = find(d2);

% gen tri2edga tri2edg
tri2edga = reshape([1:size(tri,1)*3]',size(tri,1),3);
tri2edg = edga2edg(tri2edga);

if nargout <= 3
	return;
end;

% gen nd2edg
nd2edg = inv_table(edg);

if nargout == 8
	return;
end;
% gen nd2tri
nd2tri = inv_table(tri);


function [fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2,nd2fac,nd2edg,nd2tri] = bks_all3D(tri);
%gen faces
faca = [tri(:,1) tri(:,2) tri(:,4); ...
        tri(:,2) tri(:,1) tri(:,3); ...
        tri(:,3) tri(:,4) tri(:,2); ...
        tri(:,4) tri(:,3) tri(:,1)];
faca2tri = [1:size(tri,1) 1:size(tri,1) 1:size(tri,1) 1:size(tri,1)]';

[fac,I] = sort(faca,2); %biggest in last column, smaller in first
I = even_permut(I);
[fac,Is] = sortrows(fac);
tris = faca2tri(Is);
I = I(Is);
d =[true; or(or(fac(1:end-1,1)~=fac(2:end,1), ...
	       fac(1:end-1,2)~=fac(2:end,2)), ...
	       fac(1:end-1,3)~=fac(2:end,3))]; Nd = sum(d);
%d =[true; any(fac(1:end-1,:)~=fac(2:end,:),2)];  Nd = sum(d); %slower alternative
ad = not(d); Nad = sum(ad);
d2 = [ad(2:end); false]; d2 = d2(d);
fac = fac(d,:);

% gen fac2tri
tris2 = zeros(Nd,1); tris2(d2) = tris(ad); 
fac2tri = [tris(d) tris2];

%gen fac2tri2, the element to which the permutation of the face belongs 
fac2tri2 = faca2tri(Is(d));
%fac2tri2 = fac2tri(:,1); fac2tri2(i) = fac2tri(i,2);
%flip fac2tri according to fac permutation
i = and(d2,I(d)); fac2tri(i,:) = fac2tri(i,[2 1]); 

%gen faca2fac
faca2fac = zeros(size(faca,1),1);
faca2fac(Is(d))  = 1:size(fac,1);
faca2fac(Is(ad)) = find(d2);
% gen tri2faca, tri2fac
tri2faca = reshape([1:size(tri,1)*4]',size(tri,1),4);
tri2fac = faca2fac(tri2faca);

if nargout <= 3
	return;
end;

%gen edges
edga = [faca(:,2) faca(:,1); ...
        faca(:,3) faca(:,2); ...
        faca(:,1) faca(:,3)];
edga2faca = [1:size(faca,1) 1:size(faca,1) 1:size(faca,1)]';
edga2tri  = [faca2tri; faca2tri; faca2tri];
faca2edga = reshape(1:size(edga,1),size(faca,1),3);
fac2edga = faca2edga(Is(d),:);

edg = sort(edga,2); %biggest in last column
[edg,Ise] = sortrows(edg);
facs = edga2faca(Ise);
tris = edga2tri(Ise);
de =[true; or(edg(1:end-1,1)~=edg(2:end,1),...
              edg(1:end-1,2)~=edg(2:end,2))]; Nde = sum(de);
%de = [true; any(edg(1:end-1,:)~=edg(2:end,:),2)]; Nde = sum(de); %slower alternative
ade = not(de); Nad = sum(ade);
de2 = cumsum(de); %de2 = [ade(2:end); false]; de2 = de2(de);

edg = edg(de,:);
%gen edg2faca
edgs = cumsum(de); %edgs = zeros(size(edga,1),1); edgs(de) = 1; edgs = cumsum(edgs);
edg2faca = rpval2M(edgs,facs); edg2faca = rpval2M_clean(edg2faca);
%gen edg2fac
edg2faca_ = edg2faca; edg2faca_(edg2faca_==0) = numel(faca2fac)+1;
faca2fac_ = [faca2fac; 0];
edg2fac = faca2fac_(edg2faca_); edg2fac = rpval2M_clean(edg2fac);
%gen edg2tri
edg2tri = rpval2M(edgs,tris); edg2tri = rpval2M_clean(edg2tri);

%gen edga2edg
edga2edg = zeros(size(edga,1),1);
edga2edg(Ise(de))  = 1:size(edg,1);
edga2edg(Ise(ade)) = de2(ade); %find(de2);
% gen tri2edga, tri2edg
tri2edga = reshape([1:size(tri,1)*12]',size(tri,1),12);
tri2edg = reshape(edga2edg(tri2edga(repmat([true false true false true true false false true false true false],size(tri,1),1))),size(tri,1),6);
%tri2edg = sort(edga2edg(tri2edga)');
%tri2edg = reshape(tri2edg([true(1,size(tri,1); tri2edg(1:end-1,:)~= tri2edg(2:end,:)]),6,size(tri,1))';

%gen faca2edg, fac2edg (fac2edga does not make sense)
faca2edg = edga2edg(faca2edga);
fac2edg = faca2edg(Is(d),:);

if nargout == 21
	return;
end;

% gen nd2fac
nd2fac = inv_table(fac);
%fact = fac;
%[nds,I] = sort(fact(:));
%facs = [1:size(fac,1) 1:size(fac,1) 1:size(fac,1)]'; facs = facs(I);
%nd2fac =  rpval2M(nds,facs);

if nargout == 22
	return;
end;
% gen nd2edg
nd2edg = inv_table(edg);
%edgt = edg;
%[nds,I] = sort(edgt(:));
%edgs = [1:size(edg,1) 1:size(edg,1)]'; edgs = edgs(I);
%nd2edg =  rpval2M(nds,edgs);

if nargout == 23
	return;
end;
% gen nd2tri
nd2tri = inv_table(tri);
%trit = tri;
%[nds,I] = sort(trit(:));
%tris = [1:size(tri,1) 1:size(tri,1) 1:size(tri,1)  1:size(tri,1)]'; tris = tris(I);
%nd2tri = rpval2M(nds,tris);

function is_even = even_permut(I)
is_even = false(size(I,1),1);
is_even(all(I == repmat([1,2,3],size(I,1),1),2)) = true;
is_even(all(I == repmat([3,1,2],size(I,1),1),2)) = true;
is_even(all(I == repmat([2,3,1],size(I,1),1),2)) = true;
%[1 3 2]: 0
%[1 2 3]: 1
%[2 1 3]: 0
%[2 3 1]: 1
%[3 2 1]: 0
%[3 1 2]: 1


function [tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options)
%%% INPUT/OUTPUT %%%
%tri          : N x (d+1)    contains the vertex numbers of the elements
%xy           : M X d        contains the vertex coordinates
%Nmetric      : M x (3 or 6) containts the metric field with units of inverse length (not squared). It can also
%                            contain nodal fields to be transferred to the new mesh (these fields are stored in 
%                            the columns 4 and higher or 7 and higher (2D or 3D))
%bndmesh.edg  : P x 2        contains the vertex numbers of the boundary edges in case of 2D
%bndmesh.fac  : P x 3        contains the vertex numbers of the boundary faces in case of 3D
%bndmesh.IDs  : P x 1        contains the IDs of the boundary edges/faces
%bndmesh.triID: N x 1        optional field containing the region IDs of all elements
%bndmesh.xyold: M x d        optional field containing the coordinates of the input mesh. Setting this triggerst 
%                            exact interpolation of Nmetric
%bndmesh.elv  : N x Q        optional field containing Q element wise fields to be interpolatd to the new mesh
%geomfunc     :              empty for a polygon/polehedron domain. Otherwise it is a signed distance function
%                            that returns a double J x (d+1) double where the first column is the value and the 
%                            next d is the gradient
%triQ         : N x  1       element qualities
%
%%% FUNCTIONALITY %%%
%Computes boundary mesh, if it is not input
%options.prag_adapt==0     triggers the use of newadapt, otherwise 
%options.prag_adapt == 4   triggers refinement only, while options.prag_adapt == 3 triggerst smoothing only
%setting options.forpres   triggers plotting of the inner interations in files with names given by the 
%                          string in options.forpres


if (size(xy,2) == 2 && size(Nmetric,2) == 3) || (size(xy,2) == 3 && size(Nmetric,2) == 6)
triQ = elem_qual(tri, xy, Nmetric, options);
elseif size(xy,2) == 2
triQ = elem_qual(tri, xy, Nmetric(:,1:3), options);
else
triQ = elem_qual(tri, xy, Nmetric(:,1:6), options);	
end;

minTriQ = min(triQ);
if minTriQ < options.minqual
	warning(sprintf('Low quality elements in input(%0.0e)',minTriQ));
	options.minqual = minTriQ;
	if minTriQ < 0
		error('inverted element in input mesh');
	end;
end;

if isfield(bndmesh,'triID')
  bndmesh = bndmesh_triID(bndmesh,tri);
end;

if isfield(bndmesh,'edg') || isfield(bndmesh,'fac')
	if isfield(bndmesh,'crnds')
	bndmesh.crnds = [bndmesh.crnds; geom_crnds(bndmesh,1)];
	else
	bndmesh.crnds = geom_crnds(bndmesh,1);
	end; 
	bndmesh.crnds = unique(bndmesh.crnds);
end;

if size(xy,2) == 2 && not(isfield(bndmesh,'edg'))
	bndmesh = bndmesh_polygon(tri,xy,bndmesh,options);
elseif size(xy,2) == 3 && not(isfield(bndmesh,'fac'))
	bndmesh = bndmesh_polyhedron(tri,xy,bndmesh,options);
end;

if isfield(bndmesh,'xyold')
  nd2tri = inv_table(tri);
  bndmesh.trinew = nd2tri(:,1);
  bndmesh.triold = tri;
  bndmesh.ngh = bks_tri2edg2tri(tri);
  Nmetric_ = Nmetric;
end;

if options.mntn_bks
	bks = bks_init(tri);
else
	bks = [];
end;
if options.debug
	sanity_check(tri,xy,triQ,Nmetric,options);
end;
if options.prag_adapt == 4 
     for i=1:options.innerit
      [xy,tri,bndmesh,Nmetric,triQ,bks,ndone3] = adapt_add_nd(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);   
     end;
elseif options.prag_adapt == 3 %only smooth
    [xy,Nmetric,triQ,bks,ndone4] = adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options,options.innerit);
elseif options.prag_adapt
    [tri,xy,bndmesh,Nmetric,triQ,bks]  = pragadapt(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
else
[tri,xy,bndmesh,Nmetric,triQ,bks]  = newadapt(tri,xy,Nmetric,bndmesh,triQ,bks,bks,geomfunc,options);
end;
if isfield(bndmesh,'xyold');
	[Nmetric_,triO,s,Igood] = elem_interp(bndmesh.triold,bndmesh.xyold,bndmesh.trinew,xy,options,Nmetric_);
	Nmetric(Igood,:) = Nmetric_(Igood,:);
end;
if isfield(bndmesh,'elv')
 xy_ = squeeze(mean(reshape(xy(tri(:),:),size(tri,1),size(tri,2),size(xy,2)),2));
 bndmesh.trinew = bndmesh.trinew(tri(:,1));
 bndmesh.elv = bndmesh.elv(elem_find(bndmesh,true(size(tri,1),1),xy_,options),:);
end;
if options.debug
	sanity_check(tri,xy,triQ,Nmetric,options);
end;

function [tri,xy,bndmesh,Nmetric,triQ,bks]  = pragadapt(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options,again); 
	 save_plot(tri,xy,0,options,'initial mesh');
    for i=1:options.innerit
        if options.verbose == 1
	       disp_qual(triQ,tri,xy,i);
	end;
        %COARSENING
        if options.prag_adapt ~=5
        if options.RMedg
        [xy,tri,bndmesh,Nmetric,triQ,bks,ndone1] = adapt_rm_edg(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
        else
      	[tri,triQ,bks,bndmesh,ndone1,xy,Nmetric] = adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options);
        end;
        save_plot(tri,xy,(i-1)*4+2,options,sprintf('i=%0.0f, after coarsening',i));
        else
        ndone1 = 0.;
        end;
        %SWAPPING
        if options.prag_adapt ~=6
        [tri,triQ,bks,bndmesh,ndone2] = adapt_flipedg(tri,xy,Nmetric,bndmesh,triQ,bks,options); 
        save_plot(tri,xy,(i-1)*4+1,options,sprintf('i=%0.0f, after swapping',i));
        else
        ndone2=0.;
        end;
        %REFINEMENT
        if options.prag_adapt ~=5
        [xy,tri,bndmesh,Nmetric,triQ,bks,ndone3] = adapt_add_nd(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
        save_plot(tri,xy,(i-1)*4+3,options,sprintf('i=%0.0f, after refinement',i));
        else
        	ndone3 = 0.;
        end;
        %SMOOTHING
        if options.prag_adapt ~=6 && (options.prag_adapt == 2 || (options.prag_adapt==4 && mod(i,2) == 0))
            [xy,Nmetric,bndmesh,triQ,bks,ndone4] = adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options); 
                    save_plot(tri,xy,(i-1)*4+4,options,sprintf('i=%0.0f, after smoothing',i));
        else
            ndone4 = 0;
        end;
        if options.verbose == 2
  	     disp(sprintf('%2.0f%% coarsen, %2.0f%% flip, %2.0f%% refine, %2.0f%% move',ndone1*100,ndone2*100,ndone3*100,ndone4*100));
        end;
        if ndone1 == 0 && ndone2 == 0 && ndone3 == 0 && ndone4 == 0 && options.verbose ~= -1
            disp('Adaptation has converged'); break;
        end;
    end;
    if options.prag_adapt ~=6
    [xy,Nmetric,bndmesh,triQ,bks,ndone4] = adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options,options.MVit);
            save_plot(tri,xy,4*i+1,options,'after final smoothing');
    end;
    	if nargin==9
    		disp('Recalculation qualities as angles in degrees');
    		options.qualP = eps;
    		triQ = elem_angle(tri,xy,options);
    	end;
    if 0 < options.verbose
    	disp_qual(triQ,tri,xy);
    end;
    if options.qualP < 0 && nargin == 8
     options.qualP = -options.qualP;
     options.spltc = 1;
     options.smpRFN = 2;
     options.consRFN = 2;
     options.fastRFN = 0;
     options.consRM = 1;
     options.fastRM = 0;
     options.consMV = 1;
     disp(sprintf('Fixing angles larger than %0.1f degress', options.qualP*180/pi));
     triQ = elem_angle(tri,xy,options);
     disp_qual(triQ,tri,xy);
     for i=1:options.MVit
       [tri,triQ,bks,bndmesh,ndone1,xy,Nmetric] = adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options);
       [tri,triQ,bks,bndmesh,ndone2] = adapt_flipedg(tri,xy,Nmetric,bndmesh,triQ,bks,options); 
       [xy,tri,bndmesh,Nmetric,triQ,bks,ndone3] = adapt_add_nd(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
       [xy,Nmetric,bndmesh,triQ,bks,ndone4] = adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options); 
       disp_qual(triQ,tri,xy,i)
     end;
    end;

function [tri,xy,bndmesh,Nmetric,triQ,bks]  = newadapt(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options)
    options.consRM = true;
    for i=1:options.innerit
               if options.verbose
	       disp_qual(triQ,tri,xy,i);
	end;
        while true
        	   if options.RMedg
        	   	[xy,tri,bndmesh,Nmetric,triQ,bks,ndone1] = adapt_rm_edg(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
        	   else
	            [tri,triQ,bks,bndmesh,ndone1,xy,Nmetric] = adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options);
	   end;
            if ndone1 < 1e-2
                break
            end;
        end;
        [xy,tri,bndmesh,Nmetric,triQ,ndone3] = adapt_add_nd(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
        while true %for j=1:3
            [xy,Nmetric,bndmesh,triQ,bks,ndone4] = adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options); 
            if not(size(xy,2) == 3 && options.swap3D == 0)
            [tri,triQ,bks,ndone2] = adapt_flipedg(tri,xy,Nmetric,bndmesh,triQ,bks,options); 
            else
            ndone2 = 0;
            end;
            if  size(xy,2) == 3
	    [tri,triQ,bks,bndmesh,ndone2_] = adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options);
	    ndone2 = ndone2+ndone2_;
 	    end;
            if ndone2 < 1e-2 && ndone4 < 1e-2
                break
            end;
        end;
          if options.verbose
    	disp_qual(triQ,tri,xy);
     end;
        if ndone3  == 0 %&& ndone1 == 0 && ndone2 == 0
            disp('Adaptation has converged'); break;
        end;
    end;

function disp_qual(triQ,tri,xy,ii)
if nargin == 3
disp(sprintf('%6.0f nodes, %5.0f elements, %0.1f%% average quality (min=%2.1f%%)',size(xy,1),size(tri,1),mean(triQ)*100,min(triQ)*100));
else
disp(sprintf('%3.0f: %6.0f nodes, %5.0f elements, %0.1f%% average quality (min=%2.1f%%)',ii,size(xy,1),size(tri,1),mean(triQ)*100,min(triQ)*100));
end; %0.8 is maximum mean for uniform structured 2D, it is 0.7 in 3D (Vasilevski)

function save_plot(tri,xy,i4,options,intit)
if isfield(options,'forpres')
	if not(exist(options.forpres,'dir'))
		mkdir(options.forpres);
	end;
	if size(xy,2) == 2
	i4=i4-1;
	figure(); trimesh(tri,xy(:,1),xy(:,2),'color','k','marker','none'); xlim([-0.01 1.01]); ylim([-0.01 1.01]); axis('off'); axis('equal'); box('off');  title(intit,'fontsize',32);
	print([options.forpres sprintf('/%0.0f.png',i4)],'-dpng');
	close all;
	else %3D
        if numel(strfind(version,'R2010a'))~=0
            [fac,fac2tri,tri2fac] = bks_all3D(tri);
            trimesh(fac(fac2tri(:,2)==0,:),xy(:,1),xy(:,2),xy(:,3),'edgecolor','k','marker','none'); axis('off'); box('off');
            title(intit,'fontsize',32); print('-dpng',sprintf('./%s/%0.0f.png',options.forpres,i4));
        else
            save([options.forpres sprintf('/%0.0f.mat',i4)],'tri','xy','intit');
        end;
	end;
end;
function [tri,triQ,bks,bndmesh,ndone,xy,Nmetric] = adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options)
%if options.debug
	%sanity_check(tri,xy,triQ,Nmetric,options);
%end;
% GENERATE BOOKS
ndone = 0;
nd2tri = inv_table(tri);
if size(xy,2) == 3		 
        [fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2] = bks_all3D(tri);
	nd2bndfac = inv_table(bndmesh.fac);
	nd2bndfac = [nd2bndfac; zeros(size(xy,1)-size(nd2bndfac,1),size(nd2bndfac,2))];
	%nd2edg = inv_table(edg);
	bndnds_ = false(size(xy,1),1); bndnds_(fac(fac2tri(:,2)==0,:)) = true;%bndnds_(bndmesh.fac(:)) = true;
	[ngh3,ngh4] = bks_nd2tri2ndO(tri,nd2tri,bndnds_);
	ngh = bks_nd2tri2nd2(tri,nd2tri); %for colouring
	%bndnds = zeros(size(bndnds_)); bndnds(bndnds_) = find(bndnds_);
	[tmp,edg1,edg2,edg2fac,edg2ID]  = geom_crnds(bndmesh);
%	bndedg2edg = bks_bndedg2edg(edg,nd2edg,edg1);
	bndnds2_ = false(size(xy,1),1); bndnds2_(edg1) = true;
	if isfield(bndmesh,'triID') && size(edg2ID,2) > 2
		bndnds3_ = false(size(xy,1),1); bndnds3_(edg2(edg2ID(:,3)~=0,:)) = true; %belonging to three or more faces
		bndndsI_ = false(size(xy,1),1); bndndsI_(bndmesh.fac(:)) = true;	
		ngh3 = ngh3(not(bndnds3_(bndnds_)),:);
		%bndnds2_(and(bndndsI_,not(bndnds_))) = false; %killing the internal ones
		bndnds_(bndnds3_) = false;
		%bndndsI_ = and(bndndsI_,or(bndnds3_,not(bndnds_)));
		bndndsI_(bndnds_) = false; %belonging to interior face		
	end;
	if size(Nmetric,2) ~= 6
	Tmetric = Nmetric(:,7:end); Nmetric = Nmetric(:,1:6);
	end;
else
	bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.edg(:)) = true;
	if isfield(bndmesh,'triID')
		[edg,edg2tri,tri2edg] = bks_all(tri);
		bndnds2_ = bndnds_; bndnds2_(edg(edg2tri(:,2)==0,:)) = false;
		ngh = bks_nd2tri2ndO(tri,nd2tri,and(bndnds_,not(bndnds2_)));
		nd2bndedg = inv_table(bndmesh.edg);
		allIDs = bndmesh.triID(nd2tri(:,1));
	else
		ngh = bks_nd2tri2ndO(tri,nd2tri,bndnds_);
		bndnds2_ = false(size(bndnds_));
	end;
	if size(Nmetric,2) ~= 3
	Tmetric = Nmetric(:,4:end); Nmetric = Nmetric(:,1:3);
	end;
end;
%Do not remove corner nodes
nd2clr = bks_clr(ngh,options); nd2clr(bndmesh.crnds) = 0;
if isfield(options,'TS')
	if size(xy,2) == 3
		nd2clr(bndmesh.fac(bndmesh.IDs==options.TS,:)) = 0;
	else
		nd2clr(bndmesh.edg(bndmesh.IDs==options.TS,:)) = 0;
	end;
end;

if options.fastRM
    if options.volRM
      rbadnds = find_badvol(tri,xy,Nmetric,ngh,nd2tri,options);
    else
    if options.fastRM == 2
    if size(xy,2) == 2 && not(isfield(bndmesh,'triID'))
     edg = bks_all(tri);
    end;
    badedg = find_badedg(xy,Nmetric,edg,options); %FIND BAD EDGES
    rbadnds = badedg2badnds(badedg,edg,size(xy));
    elseif options.fastRM == 1
    rbadnds = find_badnds(xy,Nmetric,ngh,options);
    else
    rbadnds = (1:size(xy,1))';
    end;
    end;
    %update colours to reflect bad nodes
    rbadnds_ = true(size(ngh,1),1); rbadnds_(rbadnds) = false; 
    %rbadnds_(and(not(bndnds_),sum(ngh~=0,2)<=size(tri,2)^2-3)) = false; 
    nd2clr(rbadnds_) = 0;
end;

if options.consRM
  triQtb = inf(size(nd2tri)); I = nd2tri~=0;
  triQtb(I) = triQ(nd2tri(I)); triQtb = min(triQtb,[],2);  
else
  triQtb = [];
end;

newtri = [];
newQ = [];
newIDt = [];
delnds = false(size(nd2tri,1),1); 
ndelnds = false(size(nd2tri,1),1);
newfacs = []; newIDs = []; newtriN2 = []; %for 3D coarsening only
clrs = max(nd2clr);
for i=1:clrs
    badnds_ = nd2clr==i;
    if options.mntn_bks
    	badnds_ = and(badnds_,bks.rmnd);
    end;
    badnds = find(badnds_);
    if numel(badnds) == 0 %nothing to do here
        continue;
    end;
    if size(xy,2) == 2
    	circles = ngh(badnds,:); 
    	edg_ = [];
    else
    	spheres = ngh4(:,:,badnds_);
    	spheres = spheres(:,:,1:max(sum(spheres(1,:,:)~=0,3)));
    	Ifac = and(bndnds_,badnds_); 
    	if any(Ifac)
    	badnds_fac = find(Ifac);
    	circles = ngh3(badnds_(bndnds_),:);
%    	if numel(geomfunc) ~= 0
%    	nvec = geomfunc(xy(Ifac,:)); nvec = -nvec(:,2:4);
%    	else
    	nvec = cross(xy(circles(:,2),:)-xy(badnds_fac,:),xy(circles(:,1),:)-xy(badnds_fac,:),2); %points out
    	badbndfac = sort([badnds_fac circles(:,1:2)],2);
    	circID = bndmesh.IDs(bks_bndedg2edg(bndmesh.fac,nd2bndfac,badbndfac));
    	end;
    end;
    if size(xy,2) == 2 || any(Ifac)
    	circles = circles(:,1:max(sum(circles~=0,2)));
    end;
    if size(xy,2) == 2
     if isfield(bndmesh,'triID')
      if any(and(badnds_,bndnds2_))
    	 [circles,badnds,badIDs] = split_2D_circles(circles,badnds,allIDs,bndnds2_,bndmesh,nd2bndedg,nd2tri,tri);
      else
    	badIDs = allIDs(badnds);
      end;
     else
    	badIDs = badnds;
     end;
     [newtri_,qualityN,badnds,newtriN_,badbadnds_,newIDt_] = fill_circles(circles,Nmetric,xy,badnds,badIDs,triQtb,options,edg_);
    else %3D
     if any(Ifac)
      if any(bndnds_(badnds))
      [spheres,newfacs_,newIDs_,newtriN2_,badnds,badbadnds_] = close_spheres(circles,circID,nvec,spheres,badnds,xy,Nmetric,badnds_fac,bndnds2_,bndmesh,nd2bndfac,options);  
      newfacs  = [newfacs ; newfacs_];
      newIDs   = [newIDs  ; newIDs_];
      newtriN2 = [newtriN2; newtriN2_];
      ndelnds(badbadnds_) = true;
      badnds_(badbadnds_) = false;
      end;
     end; %any fac
     if numel(badnds) == 0 %nothing to do here
        continue;
     end;
     %mesh spheres
     if isfield(bndmesh,'triID')
    	badnds_fac2 = find(and(badnds_,bndndsI_));
    	if numel(badnds_fac2) ~= 0
    	[spheres,newfacs_,newIDs_,newtriN2_,badnds,badbadnds_] = split_spheres(spheres,tri,badnds,xy,Nmetric,badnds_fac2,bndndsI_,bndnds2_,bndmesh,nd2bndfac,nd2tri(badnds_fac2,:)',options);
   	newfacs  = [newfacs ; newfacs_];
        newIDs   = [newIDs  ; newIDs_];
    	newtriN2 = [newtriN2; newtriN2_];
    	ndelnds(badbadnds_) = true;
    	end;
    	badIDs = get_badIDs3D(spheres(:,1,:),tri,badnds,bndmesh,nd2tri(badnds,:)'); 
     else
    	badIDs = badnds;
     end;
     [newtri_,qualityN,badnds,newtriN_,badbadnds_,newIDt_] = fill_spheres(spheres,Nmetric,xy,badnds,badIDs,triQtb,options);
    end; %triangulation
    if numel(badnds) == 0 %we could not improve the situation
        continue;
    end;    
    
    if options.debug == 2 && options.area ~= 0
     %sanity_check3(badnds,newtri_,newtriN_)
     badnds = unique(badnds);
     nd2tri_ = nd2tri(badnds,:);
     [tmp,volo_] = elem_inv(tri(nd2tri_(nd2tri_~=0),:),xy);
     volo = zeros(size(nd2tri_)); volo(nd2tri_~=0) = volo_;
     [newtriN_,I] = sort(newtriN_); newtri_ = newtri_(I,:); qualityN = qualityN(I); newIDt_ = newIDt_(I);
     [tmp,voln_] = elem_inv(newtri_,xy);
     glb2lcl = zeros(max(newtriN_),1);
     if numel(newtriN_) ~= 1
     glb2lcl(newtriN_([newtriN_(1:end-1)~=newtriN_(2:end); true])) = 1:numel(badnds);
     nd2tri_ = inv_table(glb2lcl(newtriN_))';
     else
    	glb2lcl(newtriN_) = 1; nd2tri_ = 1;
     end;
     voln = zeros(size(nd2tri_)); voln(nd2tri_~=0) = voln_;
     if any(abs(sum(voln,1)'-sum(volo,2)) > max(options.minA,1e-12))
    	save for_debug3D.mat; error('domain volume reduced in coarsening, set options.area=0 for curved geometries');
     end;
    end;
    
    %update tri, coordinate and Nmetric table
    ndelnds(badbadnds_) = true; %for active set
    delnds(badnds) = true;
    % fix tri table (new node numbers)
    newtri = [newtri; newtri_];
    newQ = [newQ; qualityN];       
    newIDt = [newIDt; newIDt_];
       
    if i < clrs %update without running gen_books
        %uncolor neighbours to badnds
        clrlssnds = ngh(badnds,:);
        nd2clr(clrlssnds(clrlssnds~=0)) = 0; 
    end;
end;
if any(delnds)
ndone = nnz(delnds);
deltri_ = nd2tri(delnds,:);
I = true(size(tri,1),1); I(deltri_(deltri_~=0)) = false;
tri = [tri(I,:); newtri];
triQ = [triQ(I,:); newQ];
if isfield(bndmesh,'triID')
bndmesh.triID = [bndmesh.triID(I); newIDt];
end;

I = true(size(xy,1),1); I(delnds) = false; xy = xy(I,:); Nmetric = Nmetric(I,:);
if exist('Tmetric','var')
Nmetric = [Nmetric Tmetric(I,:)];
end;
if isfield(bndmesh,'xyold')
      bndmesh.trinew = bndmesh.trinew(I);
end;
old2new = mvndtri(delnds); tri = old2new(tri); 
if options.mntn_bks
bks.rmnd(ndelnds) = false;
bks.mvnd(newtri) = true;
bks.rmnd(newtri) = true;
bks.mvnd = bks.mvnd(I);
bks.rmnd = bks.rmnd(I);
end;
if size(xy,2) == 3
nI = not(I); nI = nI(newtriN2); 
delbndfac = nd2bndfac(delnds,:); 
I = true(size(bndmesh.fac,1),1); I(delbndfac(delbndfac~=0)) = false;
bndmesh.fac = [bndmesh.fac(I,:); sort(newfacs(nI,:),2)];
bndmesh.IDs = [bndmesh.IDs(I,:); newIDs(nI)];
bndmesh.fac = old2new(bndmesh.fac);
else
bndmesh = update_bndmesh(bndmesh,delnds);
bndmesh.edg = old2new(bndmesh.edg);
end;
bndmesh.crnds = old2new(bndmesh.crnds);

ndone = ndone/size(ngh,1);

if options.debug
   sanity_check(tri,xy,triQ,Nmetric,options);
end;
elseif exist('Tmetric')
Nmetric = [Nmetric Tmetric];
end;

function [spheres,newfacs,newIDs,newtriN,badnds,badbadnds] = close_spheres(circles,circID,nvec,spheres,badnds,xy,Nmetric,badnds_fac,bndnds2_,bndmesh,nd2bndfac,options)
%fix bad edge nodes
IfacEDG = bndnds2_(badnds_fac);
if any(IfacEDG)
 badnds_edg = badnds_fac(IfacEDG);
[circles2,circID2,nvec2] = split_3D_circles(badnds_edg,circles(bndnds2_(badnds_fac),:),bndmesh,xy,nd2bndfac);
 I1 = 1:sum(IfacEDG); I2 = sum(IfacEDG)+1:size(circles2,1);
 circles(IfacEDG,:) = circles2(I1,:);
 circles = [circles; circles2(I2,:)];
circID(IfacEDG) = circID2(I1);
circID = [circID; circID2(I2)];
nvec(IfacEDG,:) = nvec2(I1,:);
nvec = [nvec; nvec2(I2,:)];

badnds_fac_real = [badnds_fac; badnds_edg];
else
 badnds_fac_real = badnds_fac;
end; %any edge

%%mesh boundaries (close spheres)
nvec = nvec./repmat(sqrt(sum(nvec.^2,2)),1,3);
[newfacs,qualityN,badnds_,newtriN,badbadnds,newRs,newedg,newedgN] = fill_circles(circles,Nmetric,xy,badnds_fac_real,(1:numel(circID))',[],options,nvec);
newIDs = circID(newRs);
if numel(badbadnds) ~= 0
 Ikeep = true(max(badnds),1); Ikeep(badbadnds) = false;
 spheres = spheres(:,:,Ikeep(badnds));
 badnds = badnds(Ikeep(badnds));
else
 Ikeep = true(max(badnds_fac_real),1);
end;

%CHECK INVERSION AND AREA?
if options.debug == 2 && numel(newfacs) ~= 0
 nvec = nvec(newRs,:);
 Iedg = [IfacEDG; true(sum(IfacEDG),1)];
 I = Ikeep(badnds_fac_real);
 sanity_check3(xy,circles(I,:),badnds_fac_real(I),newtriN,nvec,newfacs,options);
end; %debug


if numel(newfacs) ~= 0
 [newtriN,I] = sort(newtriN); newIDs = newIDs(I); newfacs = newfacs(I,:);
if numel(newtriN) > 1
 d = [newtriN(1:end-1)~=newtriN(2:end); true];
else
 d = true;
end;
if isfield(options,'concave3D') && numel(newedgN) ~= 0 %we cannot allow self-intersecting spheres
    if any(bndnds2_(badnds))
    I_ = [and(Ikeep(badnds_fac),bndnds2_(badnds_fac)); false(size(badnds_edg))]; NN = sum(circles~=0,2);
    newedg = [newedg; circles(I_,1) circles(find(I_)+(NN(I_)-1)*size(circles,1))];
    newedgN = [newedgN; badnds_fac(I_)];
    end;
    badbadnds_ = kill_self_intersecting(spheres, xy, badnds, newedg, newedgN);
    if numel(badbadnds_) ~= 0
        warning('COARSENING CONTRAINED DUE TO INTERSECTING BOUNDARY CAVITY');
        Ikeep = true(max(badnds),1); Ikeep(badbadnds_) = false;
        spheres = spheres(:,:,Ikeep(badnds));
        badnds = badnds(Ikeep(badnds));
        newfacs = newfacs(Ikeep(newtriN),:);
        d = d(Ikeep(newtriN));
        newIDs  = newIDs(Ikeep(newtriN));
        newtriN = newtriN(Ikeep(newtriN));
        badbadnds = [badbadnds; badbadnds_];
    end;
end;

nN = diff([0; find(d)]); nMax = max(nN);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN',nMax,1); I = C <= CnN; 
spheres2 = zeros(3,nMax,numel(badnds));
glb2lcl = zeros(max(badnds),1); glb2lcl(badnds) = 1:numel(badnds);
spheres2(:,reshape(C(I),size(newfacs,1),1)+(glb2lcl(newtriN)-1)*size(spheres2,2)) = newfacs';
spheres = [spheres spheres2];
end;


function tri = mvndtri(badnds) %reduces node numbers
tri = zeros(size(badnds)); tri(not(badnds)) = 1:nnz(not(badnds));

function badedg = find_badedg(xy,Nmetric,edg,options)
edgnowl = elem_qual(edg,xy,Nmetric,options);
badedg = find(options.Llow > edgnowl);

function rbadnds = badedg2badnds(badedg,edg,nds)
rbadnds = [];
% BAD EDGES TO BAD NODES
if numel(badedg) == 0
    return;
end;
for i = 7+23*(nds(2)==3):-1:1
    badnds = edg(badedg,:); badnds = sort(badnds(:)); badnds_ = badnds(find([badnds(1:end-1)~=badnds(2:end); true]));
    nrbadnds = badnds_(find(diff([0; find([badnds(1:end-1)~=badnds(2:end); true])]) == i)); %we take the nodes connected to 5 bad edges first, then 4 ... 
    %I_ = false(nds,1); I_(nrbadnds) = true; nrbadnds = find(I_); %nrbadnds = nrbadnds(find(not(ismember(nrbadnds,crnds)))); 
    rbadnds = [rbadnds; nrbadnds];
    I_ = false(nds(1),1); I_(rbadnds) = true; edg_ = edg(badedg,:); I = reshape(I_(edg_),size(edg_)); %I = ismember(edg(badedg,:),nrbadnds); 
    badedg = badedg(not(or(I(:,1),I(:,2))));
    if numel(badedg) == 0
        break;
    end;
end;
    
function rbadnds = find_badnds(xy,Nmetric,ngh,options)
n1 = repmat((1:size(xy,1))',1,size(ngh,2));
edg = [n1(ngh~=0) ngh(ngh~=0)];
edgnowl = elem_qual(edg,xy,Nmetric,options);
rbadnds = edg(options.Llow > edgnowl,1);
rbadnds_ = false(size(xy,1),1); rbadnds_(rbadnds) = true; %rbadnds_(crnds) = false;
rbadnds = find(rbadnds_);

function badnds_ = find_badvol(tri,xy,Nmetric,ngh,nd2tri,options)
options_ = options; options_.qualM = 10;
volM = elem_qual(tri,xy,Nmetric,options_);
volM_ = zeros(size(nd2tri)); nd2tri_ = nd2tri(nd2tri~=0);
volM_(nd2tri~=0) = volM(nd2tri_);
volM = sum(volM_,2)./sum(nd2tri~=0,2); %mean
if options.volRM == 2
 badnds_ = volM < 1.;
 return;
end;
badnds = volM < 1.;
nds = size(xy,1);
badnds_ = false(nds,1);
badndsf = zeros(nds,1); badndsf(badnds) = find(badnds);
nghI = ngh ~= 0;

while any(badnds)
 volM(not(badnds)) = inf;
 cmp = inf(size(ngh)); cmp(nghI) = volM(ngh(nghI));
 cmp2 = zeros(size(ngh)); cmp2(nghI) = badndsf(ngh(nghI));
 I = all(or(repmat(volM,1,size(ngh,2)) < cmp, and(and(repmat(volM,1,size(ngh,2))~=inf,repmat(volM,1,size(ngh,2)) == cmp),repmat(badndsf,1,size(ngh,2))>cmp2)),2); %sum(I)
 badnds_(I) = true;
 nRM = ngh(I,:); nRM = nRM(nRM(:)~=0);
 badnds(nRM(:)) = false; badnds(I) = false; badndsf(I) = 0;
end;


function bndmesh = update_bndmesh(bndmesh,badnds_)
%badnds_ = false(max([bndmesh.edg(:); badnds']),1); badnds_(badnds) = true;
affedg = badnds_(bndmesh.edg);
affedg_ = or(affedg(:,1),affedg(:,2));
if any(affedg_)
    newedg = bndmesh.edg(not(affedg_),:);
    newIDs = bndmesh.IDs(not(affedg_));
    tmp = bndmesh.edg'; affedg_badnd = tmp(affedg'); affIDs = bndmesh.IDs(affedg_);
    [affedg_badnd,I] = sort(affedg_badnd);
    oldedg = bndmesh.edg(affedg_,:); oldedg = oldedg(I,:); affIDs = affIDs(I);
    I = repmat(affedg_badnd,1,2) ~= oldedg;
    oldedg = oldedg'; 
    newedg = [newedg; sort(reshape(oldedg(I'),2,sum(affedg_)/2)',2)];
    newIDs = [newIDs; affIDs(2:2:end)];
    bndmesh.edg = newedg;
    bndmesh.IDs = newIDs;
end;

function [ocircles,circID,nvec] = split_3D_circles(badnds,circles,bndmesh,xy,nd2bndfac)
nbad = numel(badnds);
NN = sum(circles~=0,2);
badnds_ = repmat(badnds,1,size(circles,2));
Ic = circles'~=0;
[Cc,R] = find(Ic);
thmp2 = [2:size(circles,2)+1]'; thmp1 = [0:size(circles,2)-1]'; 
Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
Cr = thmp2(Cc); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
fac = zeros(numel(R),3);
fac(:,1) = circles(R+(Cc-1)*nbad);
fac(:,2) = circles(R+(Cr-1)*nbad);
fac(:,3) = badnds(R);
fac = sort(fac,2);
IDs = zeros(size(circles));
IDs(R+(Cc-1)*nbad) = bndmesh.IDs(bks_bndedg2edg(bndmesh.fac,nd2bndfac,fac));
d = and([false(nbad,1) IDs(:,1:end-1)~=IDs(:,2:end)],circles~=0); I = sum(d,2) == 2;
NN1 = zeros(nbad,1);
[shft,R_] = find(d(not(I),:)');
NN1(not(I)) = shft-1;
if any(I)
[shft,R_] = find(d(I,:)');
NN1(I) = shft(2:2:end)-shft(1:2:end);
I = repmat(I',size(circles,2),1); Ic_ = Ic(:,I(1,:)); I = I(Ic);
RI = R(I); CcI = Cc(I); 
shft1 = repmat(shft(1:2:end)'-1,size(circles,2),1); shft1 = shft1(Ic_);
CsI = CcI-shft1; Ip = CsI<1; CsI(Ip) = NN(RI(Ip))+CsI(Ip);
circles(RI+(CsI-1)*nbad) = circles(RI+(CcI-1)*nbad);
IDs(RI+(CsI-1)*nbad) = IDs(RI+(CcI-1)*nbad);
end;
circID = [IDs(:,1); IDs((1:nbad)'+NN1*nbad)];
%circID = [IDs((1:nbad)'+(NN1-1)*nbad); IDs(:,1)];
%nvec =  cross([xy(circles(:,2),:)-xy(badnds,:); ...
              %xy(circles((1:nbad)'+(NN1+1)*nbad),:)-xy(badnds,:)], ...
             %[xy(circles(:,1),:)-xy(badnds,:); ... 
              %xy(circles((1:nbad)'+NN1*nbad),:)-xy(badnds,:)],2);
%nvec =  cross([xy(circles(:,2),:)-xy(badnds,:); ...
              %xy(circles((1:nbad)'+NN1*nbad),:)-xy(badnds,:)], ...
             %[xy(circles(:,1),:)-xy(badnds,:); ... 
              %xy(circles((1:nbad)'+(NN1-1)*nbad),:)-xy(badnds,:)],2);
nvec =  cross([xy(circles(:,2),:)-xy(badnds,:); ...
              xy(circles(:,1),:)-xy(badnds,:)], ...
             [xy(circles(:,1),:)-xy(badnds,:); ... 
              xy(circles((1:nbad)'+(NN-1)*nbad),:)-xy(badnds,:)],2);
I2 = NN1(R)<Cc; I1 = Cc<=1+NN1(R); 
ocircles = zeros(nbad*2,size(circles,2));
ocircles(R(I1)+(Cc(I1)-1)*nbad*2) = circles(R(I1)+(Cc(I1)-1)*nbad);
ocircles(R(I2)+(Cl(I2)-1)*nbad*2+nbad) = circles(R(I2)+(Cc(I2)-1)*nbad);
ocircles(nbad+1:end,end) = circles(:,1);
%ocircles(nbad+1:end,:) = fix_circles(ocircles(nbad+1:end,:));


function [tri,triQ,bks,bndmesh,ndone] = adapt_flipedg(tri,xy,Nmetric,bndmesh,triQ,bks,options)
% GENERATE BOOKS
ndone = 0;
if size(xy,2) == 2 || options.swap3D
if size(xy,2) == 3
[fac,fac2tri,tri2fac] = bks_all3D(tri);
edg = fac;
nd2edg = inv_table(fac);
edg2tri = fac2tri;
tri2edg = tri2fac;
nbndedg = bks_bndedg2edg(fac,nd2edg,bndmesh.fac);
if size(Nmetric,2) ~= 6
Nmetric = Nmetric(:,1:6);
end;
else
	[edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg] = bks_all(tri);
	nbndedg = bks_bndedg2edg(edg,nd2edg,bndmesh.edg);
	if size(Nmetric,2) ~= 3
	Nmetric = Nmetric(:,1:3);
	end;
end;
%calculate ngh
edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg);
%boundary edges

%color
edg2clr = bks_clr(edg_ngh,options); edg2clr(nbndedg) = 0;
edgalt = bks_edgalt(edg,tri,edg2tri,nbndedg);
clrs = max(edg2clr); %size(edgclr,1);
edgs = size(edgalt,1);
%calculate boundary edge numbers (could be avoided by updating edg instead of regenerating it)

newtri = [];
newQ = [];
newIDs = [];
flipedg = false(size(edg,1),1);

for i=1:clrs
    % FIND EDGES TO FLIP    
    edgi = edg2clr==i;
    if not(any(edgi)) %nothing to do here(should almost never trigger)
    	continue;
    end;
    [edgi,newtri_,Nquality] = find_flipedg(xy,edg,edg2tri,edgalt,edgi,Nmetric,triQ,options);
    if any(edgi) 
            flipedg(edgi) = true;
            newtri = [newtri; newtri_];
            newQ = [newQ; Nquality];
            if isfield(bndmesh,'triID')
            	newIDs = [newIDs; repmat(bndmesh.triID(edg2tri(edgi,1)),size(edg,2),1)];
            end;
            if i<clrs %we are uncolouring the affected edges
                clrlssedg = edg_ngh(edgi,:);
                edg2clr(clrlssedg(clrlssedg ~= 0)) = 0;
            end;
            if options.debug
            	 sanity_check(tri,xy,triQ,Nmetric,options);
            end;
    end;
end;
if any(flipedg)
Nedgi = nnz(flipedg);
ndone = Nedgi;
deltri = reshape(edg2tri(flipedg,:),2*Nedgi,1);
tri(deltri,:) = newtri(1:2*Nedgi,:);
triQ(deltri) = newQ(1:2*Nedgi);
if isfield(bndmesh,'triID')
	bndmesh.triID(deltri) = newIDs(1:2*Nedgi);
end;
if size(edg,2) == 3
tri = [tri; newtri(2*Nedgi+1:end,:)];
triQ = [triQ; newQ(2*Nedgi+1:end)];
if isfield(bndmesh,'triID')
	bndmesh.triID = [bndmesh.triID; newIDs(2*Nedgi+1:end)];
end;
end;
if options.mntn_bks
bks.mvnd(newtri) = true;
bks.rmnd(newtri) = true;
%newedg = edgalt(edgi,:); delfac = find(edgi);
%if size(edg,2) == 3
%newfac = [repmat(newedg,3,1) reshape(edg(edgi,:),Nedgi*3,1)];
%end;
end;
ndone = ndone / size(edg,1);
if options.debug
    sanity_check(tri,xy,triQ,Nmetric,options);
end;
end;
end; %swap3D
if size(xy,2) == 3
	[tri,triQ,bks,bndmesh,ndone_] = adapt_flipedg3D(tri,xy,Nmetric,bndmesh,triQ,bks,options);
	ndone = ndone + ndone_;
end;

function [edgi,newtri,Nquality] = find_flipedg(xy,edg,edg2tri,edgalt,edgi,Nmetric,triQ,options);
dim = size(edg,2);
edgs = nnz(edgi);
if dim == 2	
ntri = [[edgalt(edgi,[2 1]) edg(edgi,1)]; [edgalt(edgi,:) edg(edgi,2)]];
else %3D
ntri = [[edgalt(edgi,:) edg(edgi,[2 1])]; [edgalt(edgi,:) edg(edgi,[3 2])];  [edgalt(edgi,:) edg(edgi,[1 3])]];
end;

qualityO = triQ(edg2tri(edgi,:)); 
qualityOmin = min(reshape(qualityO,edgs,2),[],2);
qualityN = elem_qual(ntri,xy,Nmetric,options,1); 
qualityNmin = min(reshape(qualityN,edgs,dim),[],2);

I = qualityNmin*options.minchg > qualityOmin;
edgi(edgi) = I;
I = repmat(I,dim,1);
newtri = ntri(I,:);
Nquality = qualityN(I,:);


function [tri,triQ,bks,bndmesh,ndone] = adapt_flipedg3D(tri,xy,Nmetric,bndmesh,triQ,bks,options)
if options.debug
	sanity_check(tri,xy,triQ,Nmetric,options);
end;
% GENERATE BOOKS
ndone = 0;
[fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2,nd2fac,nd2edg,nd2tri] = bks_all3D(tri);
nd2tri = edg2tri;
[tmp,edg1,edg2,edg2fac,edg2ID] = geom_crnds(bndmesh); 
bndedg_ = false(size(edg,1),1); 
bndedg22edg = bks_bndedg2edg(edg,nd2edg,edg2);
bndedg2edg = bks_bndedg2edg(edg,nd2edg,edg1);
bndedg_(bndedg22edg) = true;
if not(isfield(bndmesh,'triID'))
	ngh3 = bks_nd2tri2ndO(tri,nd2tri,bndedg_,edg);
else
	truebndedg = false(size(edg,1),1);
	truebndedg(fac2edg(fac2tri(:,2)==0,:)) = true;
	ngh3 = bks_nd2tri2ndO(tri,nd2tri,truebndedg,edg);
end;
ngh = bks_edg2tri2edg(edg2tri,tri2edg);
%3D: Do not flip edges on the boundary of the boundary

nd2clr = bks_clr(ngh,options); nd2clr(bndedg2edg) = 0;
if isfield(options,'TS')
	bndfac2fac = bks_bndedg2edg(fac,nd2fac,bndmesh.fac);
	nd2clr(fac2edg(bndfac2fac(bndmesh.IDs==options.TS),:)) = 0;
end;
%if isfield(bndmesh,'triID')  %filthy hack
	%nd2clr(and(not(truebndedg),bndedg_)) = 0;
%end;

triQtb = inf(size(nd2tri)); I = nd2tri~=0;
triQtb(I) = triQ(nd2tri(I)); triQtb = min(triQtb,[],2);
  
newtri = [];
newIDs = [];
newQ = [];
delnds = false(size(nd2tri,1),1); 
ndelnds = false(size(nd2tri,1),1);
clrs = max(nd2clr);
for i=1:clrs
    badnds_ = nd2clr==i; 
    badnds  = find(badnds_);
    if numel(badnds) == 0 %nothing to do here
        continue;
    end;
    circles = ngh3(badnds,:);
    circles = circles(:,1:max(sum(circles~=0,2)));
    if isfield(bndmesh,'triID')
    	badnds2_ = and(not(truebndedg(badnds_)),bndedg_(badnds_));
    	if any(badnds2_)
    		[circles,badnds] = split_circles_3D(circles,tri,edg(badnds,:),badnds,badnds2_,bndmesh,nd2tri(badnds,:));
    	end;
    	badIDs = get_badIDs(circles,tri,edg(badnds,:),bndmesh,nd2tri(badnds,:)');
    else
    	badIDs = badnds;
    end;
    [newtri_,qualityN,badnds,newtriN_,badbadnds_,newIDs_] = fill_circles(circles,Nmetric,xy,badnds,badIDs,triQtb,options,edg(badnds,:));
    
    if numel(badnds) == 0 %we could not improve the situation
        continue;
    end;    
    %update tri, coordinate and Nmetric table
    ndelnds(badbadnds_) = true;
    delnds(badnds) = true;
    % fix tri table (new node numbers)
    newIDs = [newIDs; newIDs_];
    newtri = [newtri; newtri_];
    newQ = [newQ; qualityN];       
       
    if i < clrs %update without running gen_books
        %uncolor neighbours to badnds
        clrlssnds = ngh(badnds,:);
        nd2clr(clrlssnds(clrlssnds~=0)) = 0; 
    end;
end;
if any(delnds)
ndone = nnz(delnds);
deltri_ = nd2tri(delnds,:);
I = true(size(tri,1),1); I(deltri_(deltri_~=0)) = false;
tri = [tri(I,:); newtri];
triQ = [triQ(I,:); newQ];
if isfield(bndmesh,'triID');
bndmesh.triID = [bndmesh.triID(I); newIDs];
end;

bndmesh = update_bndmesh_3D(bndmesh,delnds,edg2fac,edg2,edg2ID, bndedg2edg,bndedg22edg);

if options.mntn_bks
bks.mvnd(newtri) = true;
bks.rmnd(newtri) = true;
end;
ndone = ndone/size(ngh,1);

if options.debug
   sanity_check(tri,xy,triQ,Nmetric,options);
end;
end;


function bndmesh =  update_bndmesh_3D(bndmesh,deledg_,edg2fac,edg2,edg2ID,bndedg2edg,bndedg22edg)
%deledg_ = false(max(max(deledg),max(bndedg22edg)),1); deledg_(deledg) = true;
flipedg = deledg_(bndedg22edg);
if not(any(flipedg))
	return;
end;
DMcnv = zeros(max(edg2(:)),1); DMcnv(edg2(:)) = 1; DMcnv(DMcnv~=0) = 1:nnz(DMcnv);
I_ = false(size(deledg_,1),1); I_(bndedg2edg) = true;
edgalt = bks_edgalt(DMcnv(edg2),DMcnv(bndmesh.fac),edg2fac,find(I_(bndedg22edg)));
DMcnvINV = find(DMcnv);
edgalt(edgalt(:,1)~=0,:) = DMcnvINV(edgalt(edgalt(:,1)~=0,:));
flipfac = edg2fac(flipedg,[1 2]);
newfac = [[edg2(flipedg,1) edgalt(flipedg,:)]; [edg2(flipedg,2) edgalt(flipedg,:)]];
newIDs = repmat(edg2ID(flipedg,1),2,1);
bndmesh.fac(flipfac(:),:) = sort(newfac,2);
bndmesh.IDs(flipfac(:),:) = newIDs;


function [circles,badnds] = split_circles_3D(circles,tri,edg,badnds,badnds2_,bndmesh,nd2tri_);
badnds2 = badnds(badnds2_);
nd2tri_ = nd2tri_(badnds2_,:)';
nbad = numel(badnds2);
edg = edg(badnds2_,:);
circles_ = circles(badnds2_,:);
NN = sum(circles_~=0,2);
badnds_ = repmat(badnds2,1,size(circles,2));
Ic = circles_'~=0;
[Cc,R] = find(Ic);
thmp2 = [2:size(circles_,2)+1]'; thmp1 = [0:size(circles_,2)-1]'; 
Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
Cr = thmp2(Cc); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
triC = zeros(numel(R),4);
triC(:,1) = circles_(R+(Cc-1)*nbad); triC(:,3) = edg(R,1);
triC(:,2) = circles_(R+(Cr-1)*nbad); triC(:,4) = edg(R,2);
triC = sort(triC,2);
tri_ = zeros([4 size(nd2tri_,1) nbad]);
tri_(:,nd2tri_~=0) = tri(nd2tri_(nd2tri_~=0),:)';
tri_ = sort(tri_);
%tri1 = reshape(squeeze(tri_(1,:,R)),size(nd2tri_,1),numel(R)); %reshape active for nbad==1
tri1 = squeeze(tri_(1,:,R)); tri2 = squeeze(tri_(2,:,R)); 
tri3 = squeeze(tri_(3,:,R)); tri4 = squeeze(tri_(4,:,R));
[C,R_] = find(and(and(tri1 == repmat(triC(:,1)',size(nd2tri_,1),1),tri2 == repmat(triC(:,2)',size(nd2tri_,1),1)),and(tri3 == repmat(triC(:,3)',size(nd2tri_,1),1),tri4 == repmat(triC(:,4)',size(nd2tri_,1),1))));
triF = nd2tri_(C+(R(R_)-1)*size(nd2tri_,1));
IDs = zeros(size(circles_)); 
IDs(R+(Cc-1)*nbad) =  bndmesh.triID(triF);
d = and([false(nbad,1) IDs(:,1:end-1)~=IDs(:,2:end)],circles_~=0); I = sum(d,2) == 2;
NN1 = zeros(nbad,1);
[shft,R_] = find(d(not(I),:)');
NN1(not(I)) = shft-1;
if any(I)
[shft,R_] = find(d(I,:)');
NN1(I) = shft(2:2:end)-shft(1:2:end);
I = repmat(I',size(circles_,2),1); Ic_ = Ic(:,I(1,:)); I = I(Ic);
RI = R(I); CcI = Cc(I); 
shft1 = repmat(shft(1:2:end)'-1,size(circles_,2),1); shft1 = shft1(Ic_);
CsI = CcI-shft1; Ip = CsI<1; CsI(Ip) = NN(RI(Ip))+CsI(Ip);
circles_(RI+(CsI-1)*nbad) = circles_(RI+(CcI-1)*nbad);
IDs(RI+(CsI-1)*nbad) = IDs(RI+(CcI-1)*nbad);
end;
I2 = NN1(R)<Cc; I1 = Cc<=1+NN1(R); 
ocircles = zeros(nbad*2,size(circles_,2));
ocircles(R(I1)+(Cc(I1)-1)*nbad*2) = circles_(R(I1)+(Cc(I1)-1)*nbad);
ocircles(R(I2)+(Cl(I2)-1)*nbad*2+nbad) = circles_(R(I2)+(Cc(I2)-1)*nbad);
ocircles(nbad+1:end,end) = circles_(:,1);
ocircles = circles_fix(ocircles);
badnds = [badnds; badnds2];
circles(badnds2_,:) = ocircles(1:nbad,:);
circles = [circles; ocircles(nbad+1:end,:)];

function badIDs = get_badIDs(circles,tri,badedg,bndmesh,nd2tri_)
nbad = size(badedg,1);
triC = sort([badedg circles(:,[1 2])],2);
tri_ = zeros([4 size(nd2tri_,1) nbad]);
tri_(:,nd2tri_~=0) = tri(nd2tri_(nd2tri_~=0),:)';
tri_ = sort(tri_);
tri1 = reshape(squeeze(tri_(1,:,:)),size(nd2tri_,1),nbad); %reshape active for nbad==1
tri2 = reshape(squeeze(tri_(2,:,:)),size(nd2tri_,1),nbad);
tri3 = reshape(squeeze(tri_(3,:,:)),size(nd2tri_,1),nbad);
tri4 = reshape(squeeze(tri_(4,:,:)),size(nd2tri_,1),nbad);
[C,R] = find(and(and(tri1 == repmat(triC(:,1)',size(nd2tri_,1),1),tri2 == repmat(triC(:,2)',size(nd2tri_,1),1)),and(tri3 == repmat(triC(:,3)',size(nd2tri_,1),1),tri4 == repmat(triC(:,4)',size(nd2tri_,1),1))));
triF = nd2tri_(C+(R-1)*size(nd2tri_,1));
badIDs = bndmesh.triID(triF);	

function circles = circles_fix(circles)
ocircles = circles';
[Cg,Rg] = find(ocircles~=0);
[C ,R ] = find(ocircles==0);
move = zeros(size(circles)); move(R+(C-1)*size(circles,1)) = -1; move = cumsum(move')';
Cg = Cg + reshape(move(Rg+(Cg-1)*size(circles,1)),size(Cg));
circles = zeros(size(circles));
circles(Rg+(Cg-1)*size(circles,1)) = ocircles(ocircles~=0);


function [xy,tri,bndmesh,Nmetric,triQ,bks,ndone] = adapt_add_nd(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options)
%GENERATE BOOKS
if size(xy,2) == 3
	[xy,tri,bndmesh,Nmetric,triQ,bks,ndone] = adapt_add_nd3D(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
	return;
end;

ndone = 0;
[edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg] = bks_all(tri);
bndedg2edg = bks_bndedg2edg(edg,nd2edg,bndmesh.edg);
[edgnowL,c,wght1] = elem_qual(edg,xy,Nmetric(:,1:3),options);
splitedg = options.Lup<edgnowL;
if not(options.fastRFN)
	splitedg = true(size(splitedg));
end;

if not(any(splitedg)) 
    return
end;
if 0 < options.qualP
[triangle,badedg,badangle] = elem_angle(tri,xy,options);
badedg2edg = bks_bndedg2edg(edg,nd2edg,badedg);
splitedg = false(size(splitedg));
splitedg(badedg2edg) = true;
allangle = zeros(size(splitedg)); allangle(badedg2edg) = badangle;
edgnowL = allangle;
end;
if options.smpRFN == 1%split a single edge per element
   splitedg = splitedg1(splitedg,edgnowL,tri2edg,edg2tri);
end;
if isfield(options,'TS')
	splitedg(bndedg2edg(bndmesh.IDs==options.TS)) = false;
end;

if options.smpRFN == 2 %using colouring
[newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,xyf,ndone] = calc_new_wrap_wrap(tri,xy,Nmetric,triQ,splitedg,c,wght1,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options);
else
[newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,ndone] = calc_new_wrap(tri,xy,Nmetric,triQ,splitedg,c,wght1,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options);
xyf = find(splitedg);
end;
Nsplt = nnz(splitedg);
if Nsplt==0
	return;
end;

newnds = size(xy,1)+(1:sum(splitedg))';


if options.mntn_bks
bks.mvnd(tri(deltric~=0,:)) = true;
bks.rmnd(tri(deltric~=0,:)) = true;
bks.mvnd =  [bks.mvnd; repmat(true,Nsplt,1)];
bks.rmnd =  [bks.rmnd; repmat(true,Nsplt,1)];
end;

xy = nxy;
Nmetric = nNmetric;
%newtri = elem_fix(newtri,xy);
tri = [tri(deltric==0,:); newtri];
triQ = [triQ(deltric==0); q3]; 
if isfield(bndmesh,'triID')
	bndmesh.triID = [bndmesh.triID(deltric==0); bndmesh.triID(oldtri)];
end;
if isfield(bndmesh,'xyold')
   bndmesh.trinew = [bndmesh.trinew; elem_find(bndmesh,edg(splitedg,1),xy(newnds,:),options)];
end;


% UPDATE BOUNDARY MESH all(bndmesh.edg == sort(bndmesh.edg,2)) == true
%calculate boundary edge numbers (could be avoided by updating edg instead of regenerating it)
edgs = size(edg,1);
% delete split edges and add new ones
splitedgbnd = splitedg(bndedg2edg);
if any(splitedgbnd)
    %newnds = size(xy,1)-sum(splitedg)+(1:sum(splitedg))';
    edg2newnd = zeros(edgs,1); edg2newnd(xyf) = newnds;
    newnds = edg2newnd(bndedg2edg(splitedgbnd));
    newedg = bndmesh.edg(not(splitedgbnd),:);
    newIDs = bndmesh.IDs(not(splitedgbnd ));
    newedg = [newedg; [edg(bndedg2edg(splitedgbnd),1) newnds; edg(bndedg2edg(splitedgbnd),2) newnds]];
    newIDs = [newIDs; repmat(bndmesh.IDs(splitedg(bndedg2edg)),2,1)];
    bndmesh.edg = newedg;
    bndmesh.IDs = newIDs;
end;


if options.debug
    sanity_check(tri,xy,triQ,Nmetric(:,1:3),options);
end;

function [newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,xyf,ndone] = calc_new_wrap_wrap(tri,xy,Nmetric,triQ,splitedg,c,wght1,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options)

edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg);
edg2clr = bks_clr(edg_ngh,options); clrs = max(edg2clr);
ndone = 0; newedg = zeros(0,2); oldtri = []; nNmetric = zeros(0,size(Nmetric,2)); nxy = zeros(0,2); q3 = []; newtri = zeros(0,3); splitedgO = false(size(edg,1),1); deltric = zeros(size(tri,1),1); newnds = []; xyf = []; Nsplt = 0;
for i=1:clrs
	splitedg_ = and(splitedg,edg2clr==i);
	[newtri_,q3_,nxy_,nNmetric_,oldtri_,deltric_,newedg_,splitedg_,ndone_] = calc_new_wrap(tri,xy,Nmetric,triQ,splitedg_,c,wght1,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options);
	splitedgO(splitedg_) = true;
	fxrr = [(1:size(xy,1)) Nsplt+size(xy,1)+(1:nnz(splitedg_))]';
	newtri_ = fxrr(newtri_); newedg_ = fxrr(newedg_);
	Nsplt = Nsplt+nnz(splitedg_);
	ndone = ndone + ndone_;
	newtri = [newtri; newtri_];
	q3 = [q3; q3_];
	nxy = [nxy; nxy_(size(xy,1)+1:end,:)];
	xyf = [xyf; find(splitedg_)];
	nNmetric = [nNmetric; nNmetric_(size(xy,1)+1:end,:)];
	newedg = [newedg; newedg_];
	oldtri = [oldtri; oldtri_];
	
	deltric(deltric_~=0) = 1;
	if i ~= clrs
	edg2clr(tri2edg(oldtri_,:)) = 0;
	end;
end; %for
nxy = [xy; nxy];
nNmetric = [Nmetric; nNmetric];
splitedg = splitedgO; 

function [newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,ndone] = calc_new_wrap(tri,xy,Nmetric,triQ,splitedg,c,wght1,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options)
while true
[newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,ndone] = calc_new(tri,xy,Nmetric,splitedg,c,wght1,edg,edga,tri2edg,tri2edga,bndedg2edg,geomfunc,options);
badE = [];
goodE = true(size(q3,1),1);
if options.consRFN == 1
	goodE = q3*options.minchg >= triQ(oldtri);
elseif options.minA ~= 0 || numel(geomfunc)~=0 || options.debug
	goodE = q3 >= options.minqual;
end;
if options.consRFN == 2 && options.smpRFN
badE = calc_badT(oldtri,q3,triQ,splitedg,tri2edg,edg2tri); 
else
badE = oldtri(not(goodE));
end;
if numel(badE) ~= 0 && not(options.consRFN) && options.verbose ~= -1
	warning(sprintf('RECOMPUTING REFINEMENT due to inverted element (%1.1e,%1.1e), you might want to help us by inserting fixed nodes on your curved concave boundary.',mean(nxy(newtri(find(q3<options.minqual,1),:),:))));
end;
	
if numel(badE) ~= 0
	badedg = tri2edg(badE,:);
	if options.smpRFN
	splitedgo = splitedg;
	splitedg(badedg(:)) = false;
	splitedg_ = splitedg(splitedgo);
	nNmetric = nNmetric([true(size(xy,1),1); splitedg_],:);
	nxy = nxy([true(size(xy,1),1); splitedg_],:);
	newtri = newtri(repmat(splitedg_,2,1),:);
	oldtri = oldtri(repmat(splitedg_,2,1),:);
	q3 = q3(repmat(splitedg_,2,1));
	newedg = newedg(repmat(splitedg_,3,1),:);
	deltric = sum(splitedg(tri2edg),2);
	ndone = sum(splitedg)/size(edg,1);
	else
	splitedg(badedg(:)) = false;
	end;
else
	break;
end;
end;

function badT = calc_badT(oldtri,newq,triQ,splitedg,tri2edg,edg2tri)
edg2tri = edg2tri(splitedg,:);
edgQ = inf(size(edg2tri)); edgQ(edg2tri~=0) = triQ(edg2tri(edg2tri~=0));
edgQ_ = zeros(size(splitedg)); edgQ_(splitedg) = min(edgQ,[],2);
newq_ = min(reshape(newq,size(newq,1)/2,2),[],2);
oldtri = oldtri(1:size(oldtri,1)/2);
tri2edg = tri2edg(oldtri,:)'; edg = tri2edg(splitedg(tri2edg));
Ibad = newq_ < edgQ_(edg);
badT = oldtri(Ibad);


function [newtri,q3,xy,Nmetric,oldtri,deltric,newedg,splitedg,ndone] = calc_new(tri,xy,Nmetric,splitedg,c,wght1,edg,edga,tri2edg,tri2edga,bndedg2edg,geomfunc,options)

newcoords = c(splitedg,:);

%fix newcoords in case of curved geomtry
if numel(geomfunc)~=0
	bndedg = false(size(edg,1),1);
	bndedg(bndedg2edg) = true;
	I = bndedg(splitedg);
	dist = geomfunc(newcoords(I,:));
	newcoords(I,:) = newcoords(I,:) - dist(:,2:3).*repmat(dist(:,1),1,2);
	% we will not allow interior nodes outside the geometry:
	nI = not(I);
	dist = geomfunc(newcoords(nI,:));
	goodI = true(size(I)); goodI(nI) = dist(:,1) > options.geomtol;
	if any(not(goodI))
		if options.verbose ~= -1
		warning(sprintf('REFINEMENT CONSTRAINED DUE TO ATTEMPTED CREATION OF %0.0f INTERIOR NODE(S) OUTSIDE GEOMETRY',sum(not(goodI))));
		end;
		splitedg(splitedg) = goodI;
		newcoords = newcoords(goodI,:);
	end;
end;
ndone = sum(splitedg)/size(edg,1);

%% disp(sprintf('### MESH REFINEMENT: %0.1f%% more nodes',size(newcoords,1)));

I_ = find(splitedg); 
I2 = zeros(size(edg,1),1); I2(splitedg) = 1:numel(I_);
deltric = sum(splitedg(tri2edg),2);
deltri = I2(tri2edg);
deltri_ = deltri+size(xy,1);
%deltrit = deltri'; deltrit(deltrit==0) = numel(I_)+1; I__ = [I_;0];

%deltri contains the splitedg numbers for each element
%deltri_ contains the new node numbers for each element

%UPDATE Nmetric and calculate new coordinates
ndsO = size(xy,1);
xy = [xy; newcoords];
if size(Nmetric,2) ~= 3
Tmetric = Nmetric(:,4:end); Nmetric = Nmetric(:,1:3); 
wght1 = repmat(wght1(splitedg),1,size(Tmetric,2));
Tmetric = [Tmetric; Tmetric(edg(splitedg,1),:).*wght1 + Tmetric(edg(splitedg,2),:).*(1.-wght1)];
end;
Nmetric = [Nmetric; metric_avg(edg(splitedg,:),Nmetric,options)];
	 	
% SPLIT ELEMENTS WITH 3 LONG EDGES
I = deltric==3;
if any(I)
    N3 = nnz(I);
    if not(options.advRFN)
    	newtri = [deltri_(I,:); ...
		[tri(I,2) deltri_(I,[2 1])]; ...
		[tri(I,3) deltri_(I,[3 2])]; ...
	        [tri(I,1) deltri_(I,[1 3])]];
        q3 = elem_qual(newtri,xy,Nmetric,options);
    else   	
    nodes = [deltri_(I,:) tri(I,:)]; 
    chc = [1 2 3 5 2 1 6 3 2 4 1 3; ...
	   5 2 3 5 3 1 6 3 2 4 1 3; ...
    	   6 3 1 5 2 1 6 1 2 4 1 3; ...
    	   4 1 2 5 2 1 6 3 2 4 2 3];
    newtri = zeros(4*N3,3); q3 = zeros(size(newtri,1),1); qmin = zeros(N3,1);
    for i=1:4
	 newtri_ = [nodes(:,chc(i,[1 2 3])); nodes(:,chc(i,[4 5 6])); nodes(:,chc(i,[7 8 9])); nodes(:,chc(i,[10 11 12]))];
	 q3_ = elem_qual(newtri_,xy,Nmetric,options);
	 qmin_ = min(reshape(q3_,N3,4),[],2);
	 better = qmin < qmin_;
	 qmin(better) = qmin_(better);
	 better4 = repmat(better,4,1);
	 newtri(better4,:) = newtri_(better4,:);
	 q3(better4) = q3_(better4);
    end;%fori
    end;
    oldtri = repmat(find(I),4,1); 
    newedg = [[tri(I,2) deltri_(I,1)]; [tri(I,2) deltri_(I,2)]; ...
              [tri(I,3) deltri_(I,2)]; [tri(I,3) deltri_(I,3)]; ...
              [tri(I,1) deltri_(I,1)]; [tri(I,1) deltri_(I,3)]; ...
              newtri([false(N3,1); true(3*N3,1)],[2 3])];
else
    newtri = zeros(0,3); q3 = []; oldtri = []; newedg = zeros(0,2);
end;


% SPLIT ELEMENTS WITH 1 LONG EDGE
I = deltric==1;
if any(I)
    [C1,R_] = find(deltri(I,:)'~=0);
    C1a = C1+1; C1a(C1a==4) = 1; C1b = C1-1; C1b(C1b==0) = 3; R = find(I);
    oppnd = tri(R+size(tri,1)*(C1b-1));
    ooppnd = deltri_(R+size(tri,1)*(C1-1));
    newtri2 = [[oppnd tri(R+size(tri,1)*(C1-1)) ooppnd]; ...
               [oppnd ooppnd tri(R+size(tri,1)*(C1a-1))]];
    newtri = [newtri; newtri2];     
    q3 = [q3; elem_qual(newtri2,xy,Nmetric,options)];
    oldtri = [oldtri; repmat(R,2,1)];
    newedg = [newedg; [ooppnd oppnd]; newtri2(:,[2 3])];
end; 

% SPLITE ELEMENTS WITH 2 LONG EDGES
I = deltric==2;
if any(I)    
    [C2,R] = find(deltri(I,:)'==0);
R = find(I);
    C2b = C2-1; C2b(C2b==0) = 3;
    shrnd = tri(R+(C2b-1)*size(tri,1));
    edga2 = tri2edga(R+(C2-1)*size(tri,1));   
    deltri2_ = deltri_(I,:)';
    nshrnd = reshape(deltri2_(deltri2_~=ndsO),2,size(deltri2_,2))';
    tswitch = C2 ~= 2;
    nshrnd(tswitch,:) = nshrnd(tswitch,[2 1]);
    edgO = [edga(edga2,:) nshrnd];
    nI2 = size(edgO,1);
    tri1 = [edgO(:,[2 1 4]); edgO(:,[2 4 3])];
    tri2 = [edgO(:,[3 1 4]); edgO(:,[1 3 2])]; 
    %% DETERMINE HOW TO SPLIT
    quality1 = elem_qual(tri1,xy,Nmetric,options); 
    quality1m = min(reshape(quality1,nI2,2),[],2);
    quality2 = elem_qual(tri2,xy,Nmetric,options); 
    quality2m = min(reshape(quality2,nI2,2),[],2);
    option1_ = quality2m<quality1m; 
    option1 = repmat(option1_,2,1);
    option2 = not(option1); 
    newtri3 = [[shrnd nshrnd]; tri1(option1,:); tri2(option2,:)];
    newtri = [newtri; newtri3];
    q3 = [q3; elem_qual([shrnd nshrnd],xy,Nmetric,options); quality1(option1); quality2(option2)];
    R = find(I);
    oldtri = [oldtri; R; repmat(R(option1_),2,1); repmat(R(not(option1_)),2,1)];
    newedg = [newedg; nshrnd; newtri3([false(2*nI2,1); true(nI2,1)],[1 2])];
end;
if exist('Tmetric','var')
Nmetric = [Nmetric Tmetric];
end;


function splitedg_ = splitedg1(splitedg,edgnowL,tri2edg,edg2tri)
edgs = size(splitedg,1);
splitedg_ = false(edgs,1);
splitedgf = zeros(edgs,1); splitedgf(splitedg) = find(splitedg);
edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg);
nghI = edg_ngh~=0;
while any(splitedg)
edgnowL(not(splitedg)) = 0;
cmp = zeros(size(edg_ngh)); cmp(nghI) = edgnowL(edg_ngh(nghI));
cmp2 = zeros(size(edg_ngh)); cmp2(nghI) = splitedgf(edg_ngh(nghI));
I = all(or(repmat(edgnowL,1,size(edg_ngh,2)) > cmp,and(and(repmat(edgnowL,1,size(edg_ngh,2))~=0,repmat(edgnowL,1,size(edg_ngh,2)) == cmp),repmat(splitedgf,1,size(edg_ngh,2))>cmp2)),2); %sum(I)
splitedg_(I) = true;
nsplit = edg_ngh(I,:); nsplit = nsplit(nsplit(:)~=0);
splitedg(nsplit(:)) = false; splitedg(I) = false; splitedgf(I) = 0; 
end;
function [xy,tri,bndmesh,Nmetric,triQ,bks,ndone] = adapt_add_nd3D(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options)
%GENERATE BOOKS
[fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2,nd2fac] = bks_all3D(tri);
bndfac2fac = bks_bndedg2edg(fac,nd2fac,bndmesh.fac);
ndone = 0;
[edgnowL,c,wght1] = elem_qual(edg,xy,Nmetric(:,1:6),options);
splitedg = options.Lup<edgnowL;
if not(options.fastRFN)
	splitedg = true(size(splitedg));
end;
if not(any(splitedg)) 
    return;
end;
if 0 < options.qualP
nd2edg = inv_table(edg);
[triangle,badedg,badangle] = elem_angle(tri,xy,options);
%badedg2edg = bks_bndedg2edg(edg,nd2edg,badedg);
splitedg = false(size(splitedg));
splitedg(badedg2edg) = true;
allangle = zeros(size(splitedg)); allangle(badedg2edg) = badangle;
edgnowL = allangle;
end;
if numel(geomfunc) ~= 0
    nd2edg = inv_table(edg);
    [crnds,edg1,edg2] = geom_crnds(bndmesh); bndedg2edg = bks_bndedg2edg(edg,nd2edg,edg1);
else
    bndedg2edg = [];
end;

if options.smpRFN == 1 %split a single edge per element
   splitedg = splitedg1e(splitedg,edgnowL,tri2edg,edg2tri);
end;
if isfield(options,'TS')
	splitedg(fac2edg(bndfac2fac(bndmesh.IDs==options.TS),:)) = false;
end;

if size(Nmetric,2) ~= 6
Tmetric = Nmetric(:,7:end);
Nmetric = Nmetric(:,1:6);
end;

if options.smpRFN == 2
edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg);
edg2clr = bks_clr(edg_ngh,options); clrs = max(edg2clr);
ndone = 0; newfacs = zeros(0,3); newIDs = []; oldtri = []; nNmetric = zeros(0,size(Nmetric,2)); nxy = zeros(0,3); newq = []; newtri = zeros(0,4); splitedgO = false(size(edg,1),1); deltric = zeros(size(tri,1),1); keepbndfacs=true(size(bndmesh.IDs)); Nsplt = 0;
for i=1:clrs
	splitedg_ = and(splitedg,edg2clr==i);
        [newtri_,newq_,nxy_,nNmetric_,oldtri_,deltric_,keepbndfacs_,newfacs_,newIDs_,splitedg_,ndone_] = clc_new_wrap(tri,xy,Nmetric,triQ,splitedg_,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options);
	splitedgO(splitedg_) = true;
	fxrr = [(1:size(xy,1)) Nsplt+size(xy,1)+(1:nnz(splitedg_))]';
	newtri_ = fxrr(newtri_); newfacs_ = fxrr(newfacs_);
	Nsplt = Nsplt+nnz(splitedg_);
	ndone = ndone + ndone_;
	newtri = [newtri; newtri_];
	newq = [newq; newq_];
	nxy = [nxy; nxy_(size(xy,1)+1:end,:)];
	nNmetric = [nNmetric; nNmetric_(size(xy,1)+1:end,:)];
	newfacs = [newfacs; newfacs_];
	newIDs = [newIDs; newIDs_];
	oldtri = [oldtri; oldtri_];
	keepbndfacs(not(keepbndfacs_)) = false;
	deltric(deltric_~=0) = 1;
	if i ~= clrs
	edg2clr(tri2edg(oldtri_,:)) = 0;
	end;
end; %for
nxy = [xy; nxy];
nNmetric = [Nmetric; nNmetric];
splitedg = splitedgO; 
else
[newtri,newq,nxy,nNmetric,oldtri,deltric,keepbndfacs,newfacs,newIDs,splitedg,ndone] = clc_new_wrap(tri,xy,Nmetric,triQ,splitedg,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options);
end;
if not(any(splitedg))
	if exist('Tmetric')
	Nmetric = [Nmetric Tmetric];
	end;
	return;
end;
xy = nxy;
Nmetric = nNmetric;
bndmesh.fac = [bndmesh.fac(keepbndfacs,:); sort(newfacs,2)];
bndmesh.IDs = [bndmesh.IDs(keepbndfacs); newIDs];
if isfield(bndmesh,'triID')
	bndmesh.triID = [bndmesh.triID(deltric==0); bndmesh.triID(oldtri)];
end;
if isfield(bndmesh,'xyold')
   bndmesh.trinew = [bndmesh.trinew; elem_find(bndmesh,edg(splitedg,1),xy(size(nd2fac,1)+1:end,:),options)];
end;
if exist('Tmetric');
wght1 = repmat(wght1(splitedg),1,size(Tmetric,2));
Tmetric_ = Tmetric(edg(splitedg,1),:).*wght1 + Tmetric(edg(splitedg,2),:).*(1.-wght1);
Nmetric = [Nmetric [Tmetric; Tmetric_]];	
end;

%trio=tri;

if options.mntn_bks && any(splitedg)
bks.mvnd(tri(deltric~=0,:)) = true;
bks.rmnd(tri(deltric~=0,:)) = true;
Nsplt = size(xy,1)-size(nd2fac,1);
bks.mvnd =  [bks.mvnd; true(Nsplt,1)];
bks.rmnd =  [bks.rmnd; true(Nsplt,1)];
end;

tri = [tri(deltric==0,:); newtri];
triQ = [triQ(deltric==0); newq]; 

if options.debug
sanity_check(tri,xy,triQ,Nmetric(:,1:6),options);
if numel(geomfunc) ~= 0 && size(xy,2) == 3
bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.fac) = true;
dist = geomfunc(xy(bndnds_,:));
if max(abs(dist(:,1))) > options.geomtol
	error('boundary node far from boundary');
end;
end;
end;

function [newtri,newq,nxy,nNmetric,oldtri,deltric,keepbndfacs,newfacs,newIDs,splitedg,ndone] = clc_new_wrap(tri,xy,Nmetric,triQ,splitedg,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options)


while true
[newtri,newq,nxy,nNmetric,oldtri,deltric,keepbndfacs,newfacs,newIDs,splitedg,ndone] = clc_new(tri,xy,Nmetric,splitedg,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options);
goodE = true(size(newq,1),1); 
if options.consRFN == 1
	goodE = newq*options.minchg >= triQ(oldtri);
elseif options.minA ~= 0 || numel(geomfunc)~=0 || options.debug
	goodE = newq >= options.minqual;
end; 
if options.consRFN == 2 && options.smpRFN
badE = calc_badE(oldtri,newq,triQ,splitedg,tri2edg,edg2tri); 
%[badE,goodE,splitedg,goodnd] = calc_badE(oldtri,newq,triQ,splitedg,tri2edg,edg2tri); 
%deltric(oldtri(not(goodE))) = 0; oldtri = oldtri(goodE); newq = newq(goodE);  old2new = zeros(size(goodnd)); old2new(goodnd) = 1:sum(goodnd); old2new =  [(1:size(xy,1))'; old2new+size(xy,1)]; newtri = old2new(newtri(goodE,:)); goodnd = [true(size(xy,1),1); goodnd]; nxy = nxy(goodnd,:); nNmetric = nNmetric(goodnd,:); goodfacs = all(goodnd(newfacs),2); newfacs = old2new(newfacs(goodfacs,:)); newIDs = newIDs(goodfacs);  keepbndfacs = true(size(bndfac2fac)); keepbndfacs(any(splitedg(fac2edg(bndfac2fac,:)),2)) = false; ndone = sum(splitedg)/size(edg,1); break;
else
badE = oldtri(not(goodE));
end;

if numel(badE) ~= 0 && not(options.consRFN) && options.verbose ~= -1
	warning(sprintf('RECOMPUTING REFINEMENT due to inverted element (%1.1e,%1.1e,%1.1e), you might want to help us by inserting fixed nodes on your curved boundary.',mean(nxy(newtri(find(newq<options.minqual,1),:),:))));
end;

if numel(badE) ~= 0
	badedg = tri2edg(badE,:);
	if options.smpRFN
		deledg = false(size(splitedg)); deledg(badedg) = true; deledg(not(splitedg)) = false;
		edg2tri_ = edg2tri(deledg,:); ndeltri = true(size(tri,1),1);
		ndeltri(edg2tri_(edg2tri_~=0)) = false;
		newtri = newtri(ndeltri(oldtri),:);
		newq   =   newq(ndeltri(oldtri));
		oldtri = oldtri(ndeltri(oldtri));
		deltric(not(ndeltri)) = 0;
		I = [true(size(xy,1),1); not(deledg(splitedg))];
		nNmetric = nNmetric(I,:);
		old2new = zeros(size(nxy,1),1);
		nxy = nxy(I,:);
		old2new(I) = 1:size(nxy,1);
		newIDs  =  newIDs(I(newfacs(:,1)));
		newfacs = newfacs(I(newfacs(:,1)),:);
		keepbndfacs(any(deledg(fac2edg(bndfac2fac,:)),2)) = true;
		splitedg(deledg) = false;
		newtri(:,1)  = old2new(newtri(:,1));
		newfacs(:,1) = old2new(newfacs(:,1));
		ndone = sum(splitedg)/size(edg2fac,1);
		break;
	else
		splitedg(badedg(:)) = false;
	end;
else
	 break; 
end;
end; %while


function [newtri,newq,xy,Nmetric,oldtri,deltric,keepbndfacs,newfacs,newIDs,splitedg,ndone] = clc_new(tri,xy,Nmetric,splitedg,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options)

newcoords = c(splitedg,:);
%fix newcoords in case of curved geomtry
if numel(geomfunc)~=0 
	bndedg  = false(size(edg,1),1);
 	bndedg2 = false(size(edg,1),1);
	bndedg(fac2edg(bndfac2fac,:)) = true; 
	I = bndedg(splitedg);
	bndedg(bndedg2edg) = false;
	bndedg2(bndedg2edg) = true;
	I1 = bndedg(splitedg);
	I2 = bndedg2(splitedg);
	dist = zeros(sum(I),4);
	dist(I1(I),:) = geomfunc{1}(newcoords(I1,:));
	if any(I2) %numel(bndedg2edg) ~= 0
	dist(I2(I),:) = geomfunc{2}(newcoords(I2,:));
	end;
	newcoords(I,:) = newcoords(I,:) - dist(:,2:4).*repmat(dist(:,1),1,3);	
	% we will not allow interior nodes outside the geometry:
	nI = not(I);
	dist = geomfunc{1}(newcoords(nI,:));
	goodI = true(size(I)); goodI(nI) = dist(:,1) > -1e3*eps;
	if any(not(goodI))
		if options.verbose ~= -1
		warning(sprintf('REFINEMENT CONSTRAINED DUE TO ATTEMPTED CREATION OF %0.0f INTERIOR NODE(S) OUTSIDE GEOMETRY',sum(not(goodI))));
		end;
		splitedg(splitedg) = goodI;
		newcoords = newcoords(goodI,:);
	end;
end;

ndone = sum(splitedg)/size(edg2fac,1);

I_ = find(splitedg);
I2 = zeros(size(edg2fac,1),1); I2(splitedg) = 1:numel(I_);
deltric = sum(splitedg(tri2edg),2);
deltri = I2(tri2edg);
deltri_ = deltri+size(xy,1);
%deltri contains the splitedg numbers for each element
%deltri_ contains the new node numbers for each element

%UPDATE Nmetric
Nmetric = [Nmetric; metric_avg(edg(splitedg,:),Nmetric,options)]; 
xy = [xy; newcoords];
[xy,Nmetric,keepbndfacs,newfacs,newIDs,delfacc,fac2,fac2q] = mesh_faces(splitedg,I2,edg2fac,fac2edg,fac2edga,edga,bndfac2fac,xy,Nmetric,bndmesh,options); %,tri,fac2tri2
J_ = find(delfacc==2);
delfacc2 = zeros(size(fac,1),1); delfacc2(J_) = 1:numel(J_);

% SPLIT ELEMENTS WITH 1 LONG EDGE
%[tmp,vol] = elem_inv(tri,xy);
I = deltric==1;
if any(I)
    [newtri,newq] = mesh1long(deltri_(I,:),xy,edga,tri2edga(I,:),Nmetric,options);
    oldtri = repmat(find(I),2,1);
    %[tmp1,vol1] = elem_inv(newtri,xy); disp([1 sum(vol1) sum(vol(I))]);
else
    newtri = []; newq = []; oldtri = [];
end;

% SPLIT ELEMENTS WITH 2 LONG EDGES
I = deltric==2;
if any(I)
    [newtri2,q2,oldtri2] = mesh2long(tri(I,:),deltri_(I,:),xy,edga,tri2edga(I,:),fac2,fac2q,delfacc2,tri2fac(I,:),fac2tri2,find(I),Nmetric,options);
    %[tmp2,vol2] = elem_inv(newtri2,xy);  disp([2 sum(vol2) sum(vol(I))]);
    newtri = [newtri; newtri2]; newq = [newq; q2]; oldtri = [oldtri; oldtri2]; 
end; 

% SPLIT ELEMENTS WITH 3 LONG EDGES
I = deltric==3;
if any(I)
    [newtri3,q3,xy,Nmetric,oldtri3] = mesh3long(deltri_(I,:),xy,edga,tri2edga(I,:),fac2,fac2q,delfacc2,tri2fac(I,:),fac2tri2,find(I),Nmetric,options);
    %[tmp3,vol3] = elem_inv(newtri3,xy);  disp([3 sum(vol3) sum(vol(I))]);
    newtri = [newtri; newtri3]; newq = [newq; q3]; oldtri = [oldtri; oldtri3];
end; 

%% SPLIT ELEMENTS WITH 4 LONG EDGES
I = deltric==4; 
if any(I)
    [newtri4,q4,xy,Nmetric,oldtri4] = mesh4long(deltri_(I,:),xy,edga,tri2edga(I,:),fac2,fac2q,delfacc2,tri2fac(I,:),fac2tri2,find(I),Nmetric,options);
    %[tmp4,vol4] = elem_inv(newtri4,xy); disp([4 sum(vol4) sum(vol(I))]);
    newtri = [newtri; newtri4]; newq = [newq; q4]; oldtri = [oldtri; oldtri4];
end;

% SPLIT ELEMENTS WITH 5 LONG EDGES
I = deltric==5;
if any(I)
    [newtri5,q5,oldtri5] = mesh5long(deltri_(I,:),xy,edga,tri2edga(I,:),fac2,fac2q,delfacc2,tri2fac(I,:),fac2tri2,find(I),Nmetric,options);
    %[tmp5,vol5] = elem_inv(newtri5,xy); disp([5 sum(vol5) sum(vol(I))]);
    newtri = [newtri; newtri5]; newq = [newq; q5]; oldtri = [oldtri; oldtri5];
end;

% SPLIT ELEMENTS WITH 6 LONG EDGES
I = deltric==6;
if any(I)
    [newtri6,q6] = mesh6long(tri(I,:),deltri_(I,:),xy,Nmetric,options);
  %  [tmp6,vol6] = elem_inv(newtri6,xy); disp([6 sum(vol6) sum(vol(I))]);
    newtri = [newtri; newtri6]; newq = [newq; q6]; 
    oldtri = [oldtri; repmat(find(I),8,1)];
end; 
if exist('Tmetric','var')
Nmetric = [Nmetric Tmetric];
end;


function splitedg_ = splitedg1e(splitedg,edgnowL,tri2edg,edg2tri)
edgs = size(splitedg,1);
splitedg_ = false(edgs,1);
splitedgf = zeros(edgs,1); splitedgf(splitedg) = find(splitedg);
edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg); nghI = edg_ngh~=0;
edgnowL(not(splitedg)) = 0;
while any(splitedg)
edg_ngh_ = edg_ngh(splitedg,:);
nI = sum(splitedg);
otherL = zeros(nI,size(edg_ngh,2)); otherL(nghI(splitedg,:)) = edgnowL(edg_ngh_(nghI(splitedg,:)));
otherI = zeros(nI,size(edg_ngh,2)); otherI(nghI(splitedg,:)) = splitedgf(edg_ngh_(nghI(splitedg,:)));
selfL = repmat(edgnowL(splitedg),1,size(edg_ngh,2));
selfI = repmat(splitedgf(splitedg),1,size(edg_ngh,2));
I = false(edgs,1);
I(splitedg) = all(or(selfL > otherL,and(and(selfL~=0, selfL == otherL),selfI>otherI)),2); 
splitedg_(I) = true;
nsplit = edg_ngh(I,:); nsplit = nsplit(nsplit(:)~=0);
splitedg(nsplit(:)) = false; edgnowL(nsplit(:)) = 0;
splitedg(I) = false; splitedgf(I) = 0; 
end;

function [badE,goodE,splitedg_,goodnd] = calc_badE(oldtri,newq,triQ,splitedg,tri2edg,edg2tri)
edg2tri = edg2tri(splitedg,:);
edgQ = inf(size(edg2tri)); edgQ(edg2tri~=0) = triQ(edg2tri(edg2tri~=0));
edgQ_ = zeros(size(splitedg)); edgQ_(splitedg) = min(edgQ,[],2);
newq_ = min(reshape(newq,size(newq,1)/2,2),[],2);
oldtri = oldtri(1:size(oldtri,1)/2);
tri2edg = tri2edg(oldtri,:)'; edg = tri2edg(splitedg(tri2edg));
Ibad = newq_ < edgQ_(edg);
badE = oldtri(Ibad);
if nargout == 1
return;
end;
splitedg_ = splitedg;
splitedg_(edg(Ibad)) = false;
goodnd = splitedg_(splitedg);
goodE = repmat(splitedg_(edg),2,1);




function [abest,bbest] = surfqual(fac,indA,indB,xy,Nmetric,options)
N = size(fac,1); 
qA = elem_qual([fac(:,indA(1,:));fac(:,indA(2,:))],xy,Nmetric,options);
qB = elem_qual([fac(:,indB(1,:));fac(:,indB(2,:))],xy,Nmetric,options);
qAm = min(reshape(qA,N,2),[],2);
qBm = min(reshape(qB,N,2),[],2);
abest = qAm > qBm; bbest= not(abest);



function [xy,Nmetric,keepbndfacs,newfacs,newIDs,delfacc,fac2,fac2q] = mesh_faces(splitedg,I2,edg2fac,fac2edg,fac2edga,edga,bndfac2fac,xy,Nmetric,bndmesh,options) %,tri,fac2tri2
splitfac = false(size(fac2edg,1),1);
splitfac_ = edg2fac(splitedg,:);
splitfac(splitfac_(splitfac_~=0)) = true;
J_ = find(splitfac);
delfacc = sum(splitedg(fac2edg),2); delfacc_ = delfacc(splitfac);
delfac = I2(fac2edg(splitfac,:));
delfac_ = delfac+max(edga(:)); 

N1 = sum(delfacc==1); N2 = sum(delfacc==2); N3 = sum(delfacc==3);
keepbndfacs = true(size(bndfac2fac)); 
bndfac2fac_ = zeros(size(delfacc)); bndfac2fac_(bndfac2fac) = 1:numel(bndfac2fac);
% One split edge per face
if N1~=0

fac2edga1 = fac2edga(delfacc==1,:); %fac1_ = fac(delfacc==1,:);
delfac1_ = delfac_(delfacc_==1,:);
It = delfac(delfacc_==1,:)'~=0; [C,R] = find(It);
fac1 = zeros(N1,4);
fac1(:,1) = delfac1_(R+(C-1)*N1);
fac1(:,[2 3]) = edga(fac2edga1(R+(C-1)*N1),:);
Cp = C+1; Cp(Cp==4) = 1;
fac1(:,4) = edga(fac2edga1(R+(Cp-1)*N1),1);

%bndmesh
It = delfacc(bndfac2fac) == 1; keepbndfacs(It) = false;
bndfac2fac1_ = bndfac2fac_(delfacc==1); bndfac2fac1 = bndfac2fac1_~=0;
newfacs = [fac1(bndfac2fac1,[1 2 4]); fac1(bndfac2fac1,[1 3 4])];
newIDs = repmat(bndmesh.IDs(bndfac2fac1_(bndfac2fac1)),2,1);
else
newfacs = []; newIDs = [];
end;

% Two split edges per face
if N2~=0
delfac2_ = delfac_(delfacc_==2,:);
fac2edga2 = fac2edga(delfacc==2,:);
It = delfac(delfacc_==2,:)'~=0; [C,R] = find(It); [C1,R1] = find(not(It)); C = reshape(C,2,N2);
%cc2 = C1;
%delfac_t2 = delfac_(delfacc_==2,:)';
%fac4 = [edga(fac2edga2(R1+(C(1,:)'-1)*N2),:) edga(fac2edga2(R1+(C(2,:)'-1)*N2),:) C1];
%[all(fac4(or(C1==1,C1==3),1) == fac4(or(C1==1,C1==3),4)) all(fac4(C1==2,2) == fac4(C1==2,3))]
fac2 = zeros(N2,5); fac2(:,[3 4]) = [delfac2_(R1+(C(1,:)'-1)*N2) delfac2_(R1+(C(2,:)'-1)*N2)]; % reshape(delfac_t2(It),2,N2)'; 
It2 = C1==2; I13 = not(It2);
fac2(It2,[3 4]) = fac2(It2,[4 3]);
if any(I13)
fac2(I13,[1 2 5]) = [edga(fac2edga2(R1(I13)+(C(1,I13)'-1)*N2),2) edga(fac2edga2(R1(I13)+(C(2,I13)'-1)*N2),:)]; end;
if any(It2)
fac2(It2,[2 1 5]) = [edga(fac2edga2(R1(It2)+(C(1,It2)'-1)*N2),1) edga(fac2edga2(R1(It2)+(C(2,It2)'-1)*N2),[2 1])]; end;
qualA = elem_qual([fac2(:,[1 2 3]); ...
		 fac2(:,[4 2 3])],xy,Nmetric,options);
qualB = elem_qual([fac2(:,[2 1 4]); ...
		 fac2(:,[3 1 4])],xy,Nmetric,options);
qualAm = min(reshape(qualA,N2,2),[],2);
qualBm = min(reshape(qualB,N2,2),[],2);
fac2q = qualAm > qualBm; fac2nq = not(fac2q);

%[fac2q,fac2nq] = surfqual(fac2,[[1 2 3];[4 2 3]],[[2 1 4]; [3 1 4]],xy,Nmetric,options); 
%tri2 = tri(fac2tri2(delfacc==2),:)';
%It = or(or(tri2 == repmat(fac2(:,1)',4,1),tri2 == repmat(fac2(:,2)',4,1)),tri2 == repmat(fac2(:,5)',4,1));
%loneND = tri2(not(It));
%[tmp,vol] = elem_inv([fac2(:,[1 2 5]) loneND],xy); [tmp cc2]

%bndmesh
It = delfacc(bndfac2fac) == 2; keepbndfacs(It) = false;
bndfac2fac2_ = bndfac2fac_(delfacc==2); bndfac2fac2 = bndfac2fac2_~=0;
ItA_ = and(bndfac2fac2,fac2q ); ItA = ItA_(bndfac2fac2);
ItB_ = and(bndfac2fac2,fac2nq); ItB = ItB_(bndfac2fac2);
newfacs2 = zeros(sum(bndfac2fac2),9);
newfacs2(ItA,[1 4 7]) = fac2(ItA_,[1 2 3]);
newfacs2(ItA,[2 5 8]) = fac2(ItA_,[4 2 3]);
newfacs2(ItB,[1 4 7]) = fac2(ItB_,[2 1 4]);
newfacs2(ItB,[2 5 8]) = fac2(ItB_,[3 1 4]);
newfacs2(:,[3 6 9]) = fac2(bndfac2fac2,[3 4 5]);
newfacs = [newfacs; reshape(newfacs2,sum(It)*3,3)];
newIDs = [newIDs; repmat(bndmesh.IDs(bndfac2fac2_(bndfac2fac2)),3,1)];
else
fac2 = []; fac2q = []; cc2 = [];
end;
if N3~=0
delfac3_ = delfac_(delfacc_==3,:);
fac2edga3 = fac2edga(delfacc==3,:);
fac3 = zeros(N3,6);
fac3(:,[1 2 3]) = delfac3_;
fac3(:,[4 5]) = edga(fac2edga3(:,3),[2 1]);
fac3(:,6) = edga(fac2edga3(:,1),1);
%21=32 opposite of edga1
%31=12 opposite of edga2
%11=22 opposite of edga3

%bndmesh
It = delfacc(bndfac2fac) == 3; keepbndfacs(It) = false;
bndfac2fac3_ = bndfac2fac_(delfacc==3); bndfac2fac3 = bndfac2fac3_~=0;
newfacs = [newfacs; fac3(bndfac2fac3,[1 2 3]); fac3(bndfac2fac3,[2 3 4]); fac3(bndfac2fac3,[3 1 5]); fac3(bndfac2fac3,[1 2 6])];
newIDs = [newIDs; repmat(bndmesh.IDs(bndfac2fac3_(bndfac2fac3)),4,1)];
end;
%bndmesh.fac = [bndmesh.fac(keepbndfacs,:); sort(newfacs,2)];
%bndmesh.IDs = [bndmesh.IDs(keepbndfacs); newIDs];



function [newtri1,q1] = mesh1long(deltri1_,xy,edga,tri2edga1,Nmetric,options)
%deltri1_ = deltri_(I,:);
%tri2edga1 = tri2edga(I,:); 
tzr = max(edga(:)); N1= size(tri2edga1,1);
Ca1  = [1 3 5 6 9 11]';
deltri1_n = deltri1_';
I1 = deltri1_n~=tzr;
[C,R] = find(I1);
spltnd = deltri1_n(I1);
edgaC1 = tri2edga1(R+(Ca1(C)-1)*N1);
thmsk = [2 1 4 3 6 5]';
edgaC2 = tri2edga1(R+(Ca1(thmsk(C))-1)*N1); 
nodes = [edga(edgaC1,:) edga(edgaC2,:)];
flp = or(C==1,C==2); nodes(flp,[3 4]) = nodes(flp,[4 3]);
newtri1 = [[spltnd nodes(:,1) nodes(:,[3 4])]; ...
	  [spltnd nodes(:,2) nodes(:,[4 3])]];
q1 = elem_qual(newtri1,xy,Nmetric,options,1);

%[tmp1,vol1] = elem_inv(newtri1,xy);
%[tmp,vol] = elem_inv(tri(I,:),xy);
%[sum(vol) sum(vol1)]



function [xy,Nmetric,bndmesh,triQ,bks,ndone] = adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options,repeat)
ndone = 0;
nds = size(xy,1);
if nargin == 8
	repeat = 1;
end;
if size(xy,2) == 3
	[tmp,edg1,edg2] = geom_crnds(bndmesh); bndmesh.edg = edg1;
	bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.edg(:)) = true;
	bndnds2_ = false(size(xy,1),1); bndnds2_(bndmesh.fac(:)) = true;
	if size(Nmetric,2) ~= 6
	Tmetric = Nmetric(:,7:end);
	Nmetric = Nmetric(:,1:6);
	end;
else
	bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.edg(:)) = true;
	if size(Nmetric,2) ~= 3
	Tmetric = Nmetric(:,4:end);
	Nmetric = Nmetric(:,1:3);
	end;
end;

nd2tri = inv_table(tri);
ngh = bks_nd2tri2nd2(tri,nd2tri);

nd2clr = bks_clr(ngh,options); nd2clr(bndmesh.crnds) = 0;
if isfield(options,'TS')
	if size(xy,2) == 3
		nd2clr(bndmesh.fac(bndmesh.IDs==options.TS,:)) = 0;
	else
		nd2clr(bndmesh.edg(bndmesh.IDs==options.TS,:)) = 0;
	end;
end;
clrs = max(nd2clr);


%BOUNDARY NODES ARE AVERAGED WITH NEIGHBORING BOUNDARY NODES ONLY
if size(xy,2) == 3
bndngh = gen_bndngh3(bndmesh,bndnds2_); 
ngh(bndnds2_,1:size(bndngh,2)) = bndngh;
ngh(bndnds2_,size(bndngh,2)+1:end) = 0;
end;

if any(bndnds_)
bndngh = gen_bndngh2(bndmesh,bndnds_); ngh(bndnds_,[1 2]) = bndngh;
ngh(bndnds_,3:end) = 0;
end;
%bndngh = gen_bndngh2(bndmesh,nds);  ngh(bndnds_(1:end-1),[1 2]) = bndngh(bndnds_(1:end-1),:);

if size(xy,2) == 3
	bndnds1_ = bndnds_; bndnds_ = bndnds2_;
	bndmesh = rmfield(bndmesh,'edg');
end;

NN = sum(ngh~=0,2); 
ngh_ = ngh;
%start moving
for jj=1:repeat
for i=1:clrs
    inds = nd2clr==i;
    if options.mntn_bks
    	inds = and(inds,bks.mvnd);
    end;
    if not(any(inds))
    	continue;
    end;
    posnewi = calc_posnew(xy,ngh(inds,:),NN(inds),Nmetric,inds,options);
    if options.MVspd ~= 1.
    	posnewi = posnewi*options.MVspd+xy(inds,:)*(1-options.MVspd);
    end;
    posnew = xy; posnew(inds,:) = posnewi;
    % FIX CURVED BOUNDARY 
    indsB = and(inds,bndnds_);
    if any(indsB) && numel(geomfunc) ~= 0
    	  if size(xy,2) == 2
            dist = geomfunc(posnew(indsB,:));
            else
            indsB1 = and(indsB,not(bndnds1_));
            indsB2 = and(indsB,bndnds1_);
	  dist = zeros(sum(indsB),4);
	  dist(indsB1(indsB),:) = geomfunc{1}(posnew(indsB1,:));
	  if any(indsB2) %any(bndnds1_)
	  dist(indsB2(indsB),:) = geomfunc{2}(posnew(indsB2,:));
	  end;
            end;
            posnew(indsB,:) = posnew(indsB,:) - dist(:,2:size(xy,2)+1) .* repmat(dist(:,1),1,size(xy,2));	
    end;
    if numel(geomfunc) ~= 0
    	%we will not allow interior nodes outside the geometry:
    	indsnB = and(inds,not(bndnds_));
    	if size(xy,2) ~= 2
    	  dist = geomfunc{1}(posnew(indsnB,:));
    	else
    	  dist = geomfunc(posnew(indsnB,:));
    	end;
    	good = dist(:,1) > options.geomtol;
    	if any(not(good))
          if options.verbose ~= -1    		
    	warning(sprintf('SMOOTHING CONSTRAINED DUE TO ATTEMPTED CREATION OF %0.0f INTERIOR NODE(S) OUTSIDE GEOMETRY',sum(not(good))));
    	end;
    	inds(indsnB) = good;
    	%moveback = false(size(xy,1),1); moveback(indsnB) = not(good);
    	%posnew(moveback,:) = xy(moveback,:);
    	end;
    end;

    tomove = inds;
    
    qualold  = gen_maxedgE(xy   ,ngh(tomove,:),tri,nd2tri(tomove,:),Nmetric,bndnds_,triQ,options,1); %old element qualities

    qualnew = gen_maxedgE(posnew,ngh(tomove,:),tri,nd2tri(tomove,:),Nmetric,bndnds_,triQ,options,0); %new element qualities
    if options.consMV
	worstimproved = min(qualnew,[],2)*options.minchg > min(qualold,[],2);
    else
    	worstimproved = min(qualnew,[],2) > options.minqual; %accepting all (valid) smoothing
    end;
    %update coordinates and triQ
    tomove(tomove) = worstimproved;
    if sum(tomove) ~= 0
        ndone = ndone + sum(tomove);
     if options.MVinterp
       nd2tri_ = nd2tri(tomove,:);
       if options.MVinterp == 2 && not(exist('Tmetric'))
       	Nmetric(tomove,:) = elem_interp(tri,xy,nd2tri_,posnew(tomove,:),options,Nmetric);
       elseif options.MVinterp == 2
       	Nmetric_ = elem_interp(tri,xy,nd2tri_,posnew(tomove,:),options,[Nmetric Tmetric]);
       	Nmetric(tomove,:) = Nmetric_(:,1:3*size(xy,2)-3);
       	Tmetric(tomove,:) = Nmetric_(:,3*size(xy,2)-2:end);
       elseif exist('Tmetric')
         Tmetric(tomove,:) = elem_interp(tri,xy,nd2tri_,posnew(tomove,:),options,Tmetric);
       end;
       if(exist('Tmetric'))
         Tmetric(tomove,:) = elem_interp(tri,xy,nd2tri_,posnew(tomove,:),options,Tmetric);
       end;
     end;
     if isfield(bndmesh,'xyold')
      bndmesh.trinew(tomove) = elem_find(bndmesh,tomove,posnew(tomove,:),options);
     end;
    end;
    
    xy(tomove,:) = posnew(tomove,:);
    if options.mntn_bks
    bks.mvnd(and(inds,not(tomove))) = false;
    bks.rmnd(tomove) = true; %mvnd is already true
    chnd = ngh(tomove,:); chnd = chnd(chnd~=0);
    bks.rmnd(chnd) = true;
    bks.mvnd(chnd) = true;
    end;
    
    qualnew = qualnew(worstimproved,:);
    afftri  = nd2tri(tomove,:);
    I = afftri ~= 0;
    triQ(afftri(I)) = qualnew(I);
    if options.debug
	sanity_check(tri,xy,triQ,Nmetric,options);
	if numel(geomfunc) ~= 0 && size(xy,2) == 3
    bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.fac) = true;
    dist = geomfunc(xy(bndnds_,:));
    if max(abs(dist(:,1))) > options.geomtol
	error('boundary node far from boundary');
    end;
    end;
    end;%debug
end;%colors
end;%repeat
if exist('Tmetric')
Nmetric = [Nmetric Tmetric];
end;
ndone = ndone / size(xy,1);

function [qual,afftri] = gen_maxedgE(xye,ngh,tri,nd2tri,Nmetric,bndnds_,triQ,options,triQuse)
%qual is matrix (nds,1), that gives the qualities for elements with the nodes 
    I = nd2tri~=0; [R,C] = find(I); nds = sum(I(:));
    afftri_ = nd2tri(R+(C-1)*size(nd2tri,1));
    if triQuse
	    quality = triQ(afftri_);    	
    else
	    ntri = [reshape(tri(afftri_,:),nds,size(tri,2))];  %reshape, when nds==1?
	    quality = elem_qual(ntri,xye,Nmetric,options,1);
    end;
qual   = Inf(size(nd2tri));         qual(R+(C-1)*size(qual,1)) = quality;
if nargout == 2
    afftri = zeros(size(nd2tri)); afftri(R+(C-1)*size(qual,1)) = afftri_;
end;


function posnewi = calc_posnew(xy_,ngh_,NN_,Nmetric,inds,options)
dngh = size(ngh_);
I = ngh_~=0; 
nghI = ngh_(I);
[R,C] = find(I);
self = repmat(reshape(find(inds)',dngh(1),1),1,dngh(2));
edg_ = [reshape(self(I),sum(I(:)),1) reshape(ngh_(I),sum(I(:)),1)];
L = elem_qual(edg_,xy_,Nmetric,options);
Lo = zeros(dngh); Lo(I) = L;
x = zeros(size(ngh_)); x(I) = xy_(nghI,1);
y = zeros(size(ngh_)); y(I) = xy_(nghI,2);
posnewi = [sum(Lo.*x,2) sum(Lo.*y,2)]./repmat(sum(Lo,2),1,2);
if size(xy_,2) == 3
z = zeros(size(ngh_)); z(I) = xy_(nghI,3);
posnewi = [posnewi sum(Lo.*z,2)./sum(Lo,2)];
end;



function ngh = gen_bndngh2(bndmesh,bndnds_)
nbnd = sum(bndnds_);
bndnds = zeros(size(bndnds_,1),1); bndnds(bndnds_) = 1:nbnd; 
bndndsI = [find(bndnds_); 0];
nd2edg = inv_table(bndnds(bndmesh.edg));
nd2edg = nd2edg(:,[1 2]); %we will not move corners anyway
edg = [bndnds(bndmesh.edg); [0 NaN]]; nd2edg(nd2edg==0) = size(edg,1); %lone nodes only have a single edge
self = repmat((1:nbnd),2,1);
edg1 = edg(nd2edg(:,1),:)';
edg2 = edg(nd2edg(:,2),:)';
ngh = [edg1(edg1 ~=  self) edg2(and(edg2 ~= self, not(isnan(edg2))))]; 
ngh(ngh==0) = nbnd+1;
ngh = bndndsI(ngh);

function ngh = gen_bndngh3(bndmesh,bndnds_);
nbnd = sum(bndnds_);
bndnds = zeros(size(bndnds_,1),1); bndnds(bndnds_) = 1:nbnd; 
bndndsI = [find(bndnds_); 0];
nd2fac = inv_table(bndnds(bndmesh.fac));
ngh = bks_nd2tri2nd2(bndnds(bndmesh.fac),nd2fac);
ngh(ngh==0) = nbnd+1;
ngh = bndndsI(ngh);
function [eigL,eigR] = analyt_eig(mm,myeps)
if nargin  == 1
	myeps = eps*1e3;
end;
if size(mm,2) == 3 %2d
	eigL = zeros(size(mm,1),2);
	eigL(:,1) = 0.5*(mm(:,1)+mm(:,3)-sqrt((mm(:,1)-mm(:,3)).^2+4*mm(:,2).^2));
	eigL(:,2) = 0.5*(mm(:,1)+mm(:,3)+sqrt((mm(:,1)-mm(:,3)).^2+4*mm(:,2).^2));
	v1 = [ones(size(mm,1),1) zeros(size(mm,1),1)];
	I1 = abs(0.5*sqrt((mm(:,1)-mm(:,3)).^2+4*mm(:,2).^2))<myeps;
	I2 = abs(mm(:,2)) < myeps; 
	nI = and(not(I1),not(I2)); 
	eigL(I2,[1 2]) = mm(I2,[1 3]);
	v1(nI,:) = [-mm(nI,2) mm(nI,1)-eigL(nI,1)]; 
	v1n = v1./repmat(sqrt(sum(v1.^2,2)),1,2);
	eigR = [v1n(:,1) v1n(:,2) -v1n(:,2) v1n(:,1)]; 
else %3D
      H11 = mm(:,1);
      H12 = mm(:,2);
      H22 = mm(:,3);
      H13 = mm(:,4);
      H23 = mm(:,5);
      H33 = mm(:,6);
      p1 = H12.^2 + H13.^2 + H23.^2;
      zeroC = zeros(size(p1,1),1);
      onesC = ones(size(p1,1),1);
      eig1 = H11; eig2 = H22; eig3 = H33;
      v1 = [onesC, zeroC, zeroC];
      v2 = [zeroC, onesC, zeroC];
      v3 = [zeroC, zeroC, onesC];
      % A is not diagonal.                       
      nI = abs(p1) > myeps;
      if nnz(nI) ~= 0
      p1 = p1(nI);
      H11 = H11(nI); H12 = H12(nI); H22 = H22(nI);
      H13 = H13(nI); H23 = H23(nI); H33 = H33(nI);
      q = (H11+H22+H33)/3.;
      %H11 = H11./q; H12 = H12./q; H22 = H22./q; H13 = H13./q; H23 = H23./q; H33 = H33./q; p1 = p1./q.^2; q = ones(size(q));
      p2 = (H11-q).^2 + (H22-q).^2 + (H33-q).^2 + 2.*p1;
      p = sqrt(p2 / 6.);
      Id = [onesC,zeroC,onesC,zeroC,zeroC,onesC]; %identity matrix
      HH = [H11,H12,H22,H13,H23,H33];

      B = repmat((1./p),1,6) .* (HH-repmat(q,1,6).*Id(nI,:)); 
      %detB = B11*B22*B33+2*(B12*B23*B13)-B13*B22*B13-B12*B12*B33-B11*B23*B23
      detB = B(:,1).*B(:,3).*B(:,6)+2*(B(:,2).*B(:,5).*B(:,4))-B(:,4).*B(:,3).*B(:,4)-B(:,2).*B(:,2).*B(:,6)-B(:,1).*B(:,5).*B(:,5);
      
      %calc r
      r = detB / 2. ;
      rsmall = r<=-1.+myeps;
      rbig   = r>= 1.-myeps;
      rgood = and(rsmall==false,rbig==false);
      phi = zeros(size(H11,1),1);
      phi(rsmall) = pi / 3.;
      phi(rbig)   = 0.;
      phi(rgood)  = acos(r(rgood)) / 3.;
      
      eig1(nI) = q + 2.*p.*cos(phi);
      eig3(nI) = q + 2.*p.*cos(phi + (2.*pi/3.));
      eig2(nI) = 3.*q - eig1(nI) - eig3(nI);
      %eig1 = eig1.*q; eig2 = eig2.*q; eig3 = eig3.*q;
      v1(nI,1) = H22.*H33 - H23.^2 + eig1(nI).*(eig1(nI)-H33-H22);
      v1(nI,2) = H12.*(eig1(nI)-H33)+H13.*H23;
      v1(nI,3) = H13.*(eig1(nI)-H22)+H12.*H23;
      v2(nI,1) = H12.*(eig2(nI)-H33)+H23.*H13;
      v2(nI,2) = H11.*H33 - H13.^2 + eig2(nI).*(eig2(nI)-H11-H33);
      v2(nI,3) = H23.*(eig2(nI)-H11)+H12.*H13;
      v3(nI,1) = H13.*(eig3(nI)-H22)+H23.*H12;
      v3(nI,2) = H23.*(eig3(nI)-H11)+H13.*H12;
      v3(nI,3) = H11.*H22 - H12.^2 + eig3(nI).*(eig3(nI)-H11-H22);
      L1 = sqrt(sum((v1(nI,:).^2),2));
      L2 = sqrt(sum((v2(nI,:).^2),2));
      L3 = sqrt(sum((v3(nI,:).^2),2));
      v1(nI,:) = v1(nI,:)./repmat(L1,1,3);
      v2(nI,:) = v2(nI,:)./repmat(L2,1,3);
      v3(nI,:) = v3(nI,:)./repmat(L3,1,3);
      nI(nI) = all([L1 L2 L3]~=0,2); %,abs(1-abs(r(rgood))) > ;
      end;
      eigL = [eig1,eig2,eig3];
      eigR = [v1(:,1),v1(:,2),v1(:,3),...
              v2(:,1),v2(:,2),v2(:,3),...
              v3(:,1),v3(:,2),v3(:,3)];
      %[sum(not(nI)) sum(not(all([L1 L2 L3]~=0,2))) sum(not(sum([L1 L2 L3]==0,2)<2)) sum(not(sum(abs(mm(nI,:) -analyt_prod(analyt_fulleig(eigL(nI,:)),eigR(nI,:))),2) < myeps)) sum(not(abs(mm(:,2).^2 + mm(:,4).^2 + mm(:,5).^2) > myeps))]
      nI(nI) = sum(abs(mm(nI,:) -analyt_prod(analyt_fulleig(eigL(nI,:)),eigR(nI,:))),2) < myeps;
      for ii = find(not(nI))' %ii = 1:size(eigL,1)
      	%reshape(mm(ii,[1 2 4 2 3 5 4 5 6]),3,3)
	[v,L] = eig(reshape(mm(ii,[1 2 4 2 3 5 4 5 6]),3,3));
	eigL(ii,:) = diag(L);
	eigR(ii,:) = v(:);
      end;
end;


function out = analyt_fulleig(eigL)
zeron = zeros(size(eigL,1),1);
if size(eigL,2) == 2
	out = [eigL(:,1) zeron eigL(:,2)];
else %3D
	out = [eigL(:,1) zeron eigL(:,2) zeron zeron eigL(:,3)];
end;   
function out = analyt_prod(eigL_,eigR)
if size(eigR,2) == 4 || size(eigR,2) == 9 %rotation
	if size(eigL_,2)==3
	inds = [1 2;
	        2 3];
	indA = [1 2;
	        3 4];
	else %3D
	inds = [1 2 4;
	        2 3 5;
	        4 5 6];
	indA = [1 2 3;
	        4 5 6;
	        7 8 9];
	end;
	indB = indA';
	out = zeros(size(eigL_));
	for ii=1:size(inds,1)
	for jj=1:size(inds,1)
	for mm=1:size(inds,1)
	for nn=1:size(inds,1)
		if ii<nn
			continue
		end;
		out(:,inds(ii,nn)) = out(:,inds(ii,nn)) + eigR(:,indB(ii,jj)).*eigL_(:,inds(jj,mm)).*eigR(:,indA(mm,nn));
	end;
	end;
	end;
	end;
else %size(eigR,2) == 3 || size(eigR,2) == 6, coord transform
	if size(eigR,2)==3
	inds = [1 2;
	        2 3];
	else %3D
	inds = [1 2 4;
	        2 3 5;
	        4 5 6];
	end;
	out = zeros(size(eigL_));
	for ii=1:size(inds,1)
	for jj=1:size(inds,1)
		out(:,ii) = out(:,ii) + eigL_(:,jj).*eigR(:,inds(ii,jj));
	end;
	end;
end;


function bndmesh = bndmesh_polygon(tri,xy,bndmesh,options)
edga = [tri(:,2) tri(:,1); tri(:,3) tri(:,2); tri(:,1) tri(:,3)];
[edg,I] = sort(edga,2); %biggest in last column
[edg,Is] = sortrows(edg);
I = I(Is,1)==1; %true for non flipped
d =[or(edg(1,1)~=edg(2,1),edg(1,2)~=edg(2,2)); and(or(edg(2:end-1,1)~=edg(1:end-2,1), edg(2:end-1,2)~=edg(1:end-2,2)),...
    or(edg(2:end-1,1)~=edg(3:end,1),edg(2:end-1,2)~=edg(3:end,2))); or(edg(end,1)~=edg(end-1,1),edg(end,2)~=edg(end-1,2))]; 
edg = edg(d,:);
edg(I(d),:) = edg(I(d),[2 1]);
t = xy(edg(:,1),:)-xy(edg(:,2),:); t = t./repmat(sqrt(sum(t.^2,2)),1,2);
edg(I(d),:) = edg(I(d),[2 1]);
dist = t(:,2).*xy(edg(:,1),1) - t(:,1).*xy(edg(:,1),2);
[tmp, I] = sortrows([t dist]); edg = edg(I,:);
d = [true; any(abs(tmp(1:end-1,:)-tmp(2:end,:)) > options.geomtol,2)];
IDs = cumsum(d);

if options.verbose && options.verbose ~= -1
 disp(sprintf('Found %i co-linear edges', max(IDs)));
end;
bndmesh.edg = edg;
bndmesh.IDs = IDs;
if isfield(bndmesh,'crnds')
bndmesh.crnds = [bndmesh.crnds; geom_crnds(bndmesh)];
else
bndmesh.crnds = geom_crnds(bndmesh);
end;


function bndmesh = bndmesh_polyhedron(tri,xy,bndmesh,options)
faca = [tri(:,1) tri(:,2) tri(:,4); ...
        tri(:,2) tri(:,1) tri(:,3); ...
        tri(:,3) tri(:,4) tri(:,2); ...
        tri(:,4) tri(:,3) tri(:,1)];

[fac,I] = sort(faca,2); %biggest in last column, smaller in first
is_even = false(size(I,1),1);
is_even(all(I == repmat([1,2,3],size(I,1),1),2)) = true;
is_even(all(I == repmat([3,1,2],size(I,1),1),2)) = true;
is_even(all(I == repmat([2,3,1],size(I,1),1),2)) = true;
[fac,Is] = sortrows(fac);
is_even = is_even(Is);
d =[or(or(fac(1,1)~=fac(2,1), ...
	fac(1,2)~=fac(2,2)), ...
          fac(1,3)~=fac(2,3));
and(or(or(fac(2:end-1,1)~=fac(1:end-2,1), ...
          fac(2:end-1,2)~=fac(1:end-2,2)), ...
	fac(2:end-1,3)~=fac(1:end-2,3)), ...
    or(or(fac(2:end-1,1)~=fac(3:end,1), ...
          fac(2:end-1,2)~=fac(3:end,2)), ...
          fac(2:end-1,3)~=fac(3:end,3)));
    or(or(fac(end-1,1)~=fac(end,1), ...
          fac(end-1,2)~=fac(end,2)), ...
          fac(end-1,3)~=fac(end,3))]; 
fac = fac(d,:);
fac(is_even(d),[1 2]) = fac(is_even(d),[2 1]);
n = cross(xy(fac(:,2),:)-xy(fac(:,1),:),xy(fac(:,3),:)-xy(fac(:,1),:),2); n = n./repmat(sqrt(sum(n.^2,2)),1,3);
fac(is_even(d),[1 2]) = fac(is_even(d),[2 1]);
dist = sum(n.*xy(fac(:,1),:),2);
[tmp, I] = sortrows([n dist]); fac = fac(I,:);
d = [true; any(abs(tmp(1:end-1,:)-tmp(2:end,:)) > options.geomtol,2)];
IDs = cumsum(d);

if options.verbose && options.verbose ~= -1
 disp(sprintf('Found %i co-linear faces', max(IDs)));
end;
bndmesh.fac = fac;
bndmesh.IDs = IDs;
if isfield(bndmesh,'crnds')
bndmesh.crnds = [bndmesh.crnds; geom_crnds(bndmesh,1)];
else
bndmesh.crnds = geom_crnds(bndmesh,1);
end;

function [edg2clr,clr2edg] = bks_clr(edg_ngh,options) 
%%% ALGORITHMS
%options.greedyCLR == 3  greedy with Largest degree first post-processing
%options.greedyCLR == 2  greedy with Jones Plassmann post-processing
%options.greedyCLR == 1  greedy (colouring)
%options.greedyCLR == 0  Jones Plassmann
%options.greedyCLR == -1 Largest degree first 

if options.debug == 2 %check input
for i=1:size(edg_ngh,1)
if not(all(any(edg_ngh(edg_ngh(i,edg_ngh(i,:)~=0),:) == i,2)))
	error(sprintf('input to coloring algorithm is flawed (%0.0f)',i));
end; end; end;

edgs = size(edg_ngh,1);
if not(options.greedyCLR>0) && mean(sum(edg_ngh~=0,2)==size(edg_ngh,2)) < 0.5
	edg2clr = clr_JP(edg_ngh,options);
else
mxtrs = 70;
nghI = edg_ngh~=0;
trs = 0; clrs = size(edg_ngh,2)+(mean(sum(edg_ngh~=0,2)==size(edg_ngh,2)) > 0.5);
edg2clr = ceil(rand(edgs,1)*clrs); %1+mod(1:edgs,clrs);
while 1
    if trs > mxtrs
        trs = 0;
        clrs = clrs + 1;
        warning(sprintf('### EDGE GROUPING: GROUP SIZE INCREASED(%0.0d) ###',clrs-1));
    end;
    if trs == 0
        badedg = (1:edgs)';
    end;
    badedgclr = edg_ngh(badedg,:); badedgclr(nghI(badedg,:)) = edg2clr(badedgclr(nghI(badedg,:)));
    cmpclr = edg2clr(repmat(badedg,1,size(edg_ngh,2)));
    I = sum(cmpclr == badedgclr,2)>0; badedg = badedg(I);
    if numel(badedg) == 0
        break;
    end;
    % SUGGEST NEW COLORS TO BAD EDGES
    trs = trs + 1;
    I = find(cmpclr(I,:) == badedgclr(I,:));
    badedg_ = 1:edgs+1; badedg_(badedg) = 0; 
    tknclrs = edg_ngh(badedg,:); tknclrs(nghI(badedg,:)) = badedg_(tknclrs(nghI(badedg,:)));
    
    %do not consider colours of badedgs as taken
    tknclrs(tknclrs~=0) = edg2clr(tknclrs(tknclrs~=0));
    if max(sum(tknclrs ~= 0,2)) == clrs
        trs = mxtrs + 1;
        continue;
    end;
    Nbad = numel(badedg);
    llclrs = repmat((1:clrs),Nbad,1); 
    [R,tmp,C] = find(tknclrs); 
    llclrs(R+(C-1)*Nbad) = 0;
    vlclrs = sum(llclrs~=0,2); llclrs = sort(llclrs,2);
    C = clrs+1-ceil(rand(Nbad,1).*vlclrs)'; nwclrs = llclrs((1:Nbad)+(C-1)*Nbad);
    edg2clr(badedg) = nwclrs;
end;

if options.greedyCLR > 1
%REDUCE CROMATIC NUMBER
if options.greedyCLR == 3
Nn = sum(edg_ngh~=0,2);
end;

while true
done = 1;
for i=1:clrs
	inds = edg2clr==i; %inds = edg2clr_(1:end-1)==i;
	finds = find(inds); 
	Ni = numel(finds);
	if Ni == 0
		continue;
	end;
	tknclrs = edg_ngh(inds,:); tknclrs(nghI(inds,:)) = edg2clr(tknclrs(nghI(inds,:))); 
	tknclrs = [edg2clr(finds) tknclrs];
	llclrs =  repmat((1:clrs),Ni,1); 
	[R,C,Cc] =  find(tknclrs);
	llclrs(R+(Cc-1)*Ni) = clrs+1;
	minclr = min(llclrs,[],2);
	if options.greedyCLR == 3
	%COLOUR MOST NEIGHBOUGHS FIRST
	llclrs =  repmat((1:clrs),Ni,1); 
	llclrs(R+(Cc-1)*Ni) = 0;
	maxclr = max(llclrs,[],2);
	nNngh = edg_ngh(inds,:); nNngh(nghI(inds,:)) = Nn(nNngh(nghI(inds,:)));
	nNngh = reshape(nNngh,Ni,size(edg_ngh,2));  %reshaping when Ni==1
	I2 = repmat(Nn(finds),1,size(nNngh,2)) < nNngh;
	bgngh = zeros(size(nNngh)); bgngh(I2) = tknclrs(I2);
	I = and(and(max(bgngh,[],2)<minclr,minclr < edg2clr(finds)),Nn(finds)~=0);
	Im = and(Nn(finds)~=0,not(I));
	edg2clr(finds(Im)) = maxclr(Im);
	Nn(finds(I)) = 0;
	else
	I = minclr < edg2clr(finds);
	end;
	edg2clr(finds(I)) = minclr(I);	
	if any(I)
		done = 0;
	end;
end; %for
if done == 1
	break;
end;
end;%while true
end; %reduce cromatic number
end; %greedyCLR


if nargout == 2
    clr2edg = env_table(edg2clr);
end;
if options.debug
    % SIMPLE CHECK
    edg2clr_ = NaN(size(edg_ngh)); edg2clr_(edg_ngh~=0) = edg2clr(edg_ngh(edg_ngh~=0));
    cmpmat = repmat(edg2clr(1:edgs),1,size(edg_ngh,2)); cmpmat = edg2clr_ == cmpmat;
    if any(cmpmat(:)) %0~=sum(sum(cmpmat))
        [R,C] = find(cmpmat); 
        error(sprintf('EGDE NUMBER %0.0f HAS WRONG COLOUR',R(1)))
    end;
end;


function edg2clr = clr_JP(edg_ngh,options)
edgs = size(edg_ngh,1);
nghI = edg_ngh~=0;
nN = sum(edg_ngh~=0,2);
clrs = size(edg_ngh,2)+1;
edg2clr = zeros(edgs,1);
badedg = (1:edgs)';
while numel(badedg) ~= 0
    badedgf_ = [zeros(edgs,1)];
    badedgf_(badedg) = badedg;
    %edges to color (they have the largest number):
    badedg_ngh = edg_ngh(badedg,:); badedg_ngh(nghI(badedg,:)) = badedgf_(badedg_ngh(nghI(badedg,:)));
    if options.greedyCLR == 0
	    colorI = badedg > max(badedg_ngh,[],2); 
    else
	    %edges to color have the most neighbours (largest number in case of conflicts):
	    badedgN = badedgf_; badedgN(badedg) = nN(badedg);
	    badedg_nghN = edg_ngh(badedg,:); badedg_nghN(nghI(badedg,:)) = badedgN(badedg_nghN(nghI(badedg,:)));
	    [maxN,C] = max(badedg_nghN,[],2); I = repmat(maxN,1,size(edg_ngh,2)) == badedg_nghN; 
	    maxNedg = zeros(size(badedg_ngh)); maxNedg(I) = badedg_ngh(I); maxNedg = max(maxNedg,[],2);
	    colorI = or(nN(badedg) > maxN,and(nN(badedg) == maxN,badedg>maxNedg));
    end;
        
    badedgI = badedg(colorI); 
    badedg_ = 1:edgs; badedg_(badedgI) = 0; 
    tknclrs = edg_ngh(badedgI,:); tknclrs(nghI(badedgI,:)) = badedg_(tknclrs(nghI(badedgI,:)));
    %do not consider colours of badedgs as taken
    tknclrs(tknclrs~=0) = edg2clr(tknclrs(tknclrs~=0));

    tknclrs = [zeros(numel(badedgI),clrs-size(tknclrs,2)) tknclrs];
    R = repmat((1:numel(badedgI))',1,size(tknclrs,2)); I = find(tknclrs); C = tknclrs(I);
    tknclrs = zeros(size(tknclrs)); tknclrs(R(I)+(C-1)*numel(badedgI)) = C;
    llclrs = repmat((1:clrs),numel(badedgI),1); llclrs(find(tknclrs)) = 0;
    vlclrs = sum(llclrs>0,2); llclrs = sort(llclrs,2);
    C = clrs-vlclrs'+1; nwclrs = llclrs((1:numel(badedgI))+(C-1)*numel(badedgI));
    edg2clr(badedgI) = nwclrs;
    clrs = max(clrs,max(nwclrs)+1);
    badedgf_(badedgI) = 0; badedg = badedgf_(find(badedgf_));
end;

function bndedg2edg = bks_bndedg2edg(edg,nd2edg,bndedg)
if size(edg,2) == 2
	nd2edg_ = nd2edg(bndedg(:,1),:); I = nd2edg_ ~= 0;
	edg2 = zeros(size(nd2edg_)); edg2(I) = edg(nd2edg_(I),2);
	II = edg2 == repmat(bndedg(:,2),1,size(nd2edg,2));
	nd2edg_ = nd2edg_'; bndedg2edg = nd2edg_(II');
else
	nd2fac_ = nd2edg(bndedg(:,1),:); I = nd2fac_ ~= 0;
	fac2 = zeros(size(nd2fac_)); fac2(I) = edg(nd2fac_(I),2);
	fac3 = zeros(size(nd2fac_)); fac3(I) = edg(nd2fac_(I),3);
	II = and(fac2 == repmat(bndedg(:,2),1,size(nd2edg,2)),fac3 == repmat(bndedg(:,3),1,size(nd2edg,2)));
	nd2fac_ = nd2fac_'; bndedg2edg = nd2fac_(II');
end;
function edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg,Iedg)
tris = size(tri2edg,1);
edgs = size(edg2tri,1);
if nargin == 2
	Iedg = (1:edgs)';
end;

switch size(tri2edg,2)
case 3
I = edg2tri~=0;
edg2tri_ = edg2tri(I);
N1 = zeros(edgs,3); N1(I(:,1),:) = tri2edg(edg2tri(I(:,1),1),:);
N2 = zeros(edgs,3); N2(I(:,2),:) = tri2edg(edg2tri(I(:,2),2),:);
edg_ngh = [N1 N2];
edg_ngh(edg_ngh == repmat(Iedg,1,6)) = 0;
edg_ngh = sort(edg_ngh,2); edg_ngh = edg_ngh(:,3:6);
case 4 %3D face, fac2tri2fac
I = edg2tri~=0;
edg2tri_ = edg2tri(I);
N1 = zeros(edgs,4); N1(I(:,1),:) = tri2edg(edg2tri(I(:,1),1),:);
N2 = zeros(edgs,4); N2(I(:,2),:) = tri2edg(edg2tri(I(:,2),2),:);
fac_ngh = [N1 N2];
%fac_ngh = [tri2edg_(edg2tri_(:,1),:) tri2edg_(edg2tri_(:,2),:)];
fac_ngh(fac_ngh == repmat(Iedg,1,8)) = 0;
edg_ngh = sort(fac_ngh,2); edg_ngh = edg_ngh(:,3:8);%3D swapping
case 6 %3D edg
edg2tri = edg2tri';
tri2edg = tri2edg';
I = edg2tri~=0;
edg2tri_ = edg2tri(I);
edg_ngh = zeros(6*size(edg2tri,1),edgs); edg_ngh(repmat(I,6,1)) = tri2edg(:,edg2tri_);
edg_ngh(edg_ngh == repmat(Iedg',size(edg_ngh,1),1)) = 0;
edg_ngh = sort(edg_ngh);
I = and(edg_ngh~=0,[edg_ngh(1:end-1,:)~=edg_ngh(2:end,:); true(1,edgs)]);
nN = sum(I);
%vals = edg_ngh(I);
%R = edg_(I);
%edg_ngh = rpval2M(vals,R);
edg_ngh(not(I)) = 0;
edg_ngh = sort(edg_ngh);
edg_ngh = edg_ngh(end-max(nN)+1:end,:)';
end;
%

function bks = bks_init(tri)
nds = max(tri(:));
bks.rmnd = true(nds,1);
bks.mvnd = true(nds,1);
function ngh = bks_nd2tri2nd2(tri,nd2tri,Ind)
%SIMPLE ALTERNATIVE WITHOUT ORDERING
if nargin == 2
	Ind = 1:size(nd2tri,1);
end;
nd2tri = nd2tri';
I = nd2tri~=0; nd2tri_ = nd2tri(I);
ngh1 = zeros(size(nd2tri)); ngh1(I) = tri(nd2tri_,1);
ngh2 = zeros(size(nd2tri)); ngh2(I) = tri(nd2tri_,2);
ngh3 = zeros(size(nd2tri)); ngh3(I) = tri(nd2tri_,3);
ngh = [ngh1; ngh2; ngh3];
if size(tri,2) == 4
 ngh4 = zeros(size(nd2tri)); ngh4(I) = tri(nd2tri_,4);
 ngh = [ngh; ngh4];
end;
ngh(ngh==repmat(Ind,size(ngh,1),1)) = 0;
ngh = sort(ngh);
ngh([false(1,size(ngh,2)); ngh(1:end-1,:)==ngh(2:end,:)]) = 0;
nN = sum(ngh~=0);
nds = ngh(ngh~=0);
ngh = zeros(max(nN),size(ngh,2));
ngh(repmat((1:max(nN))',1,size(ngh,2)) <= repmat(nN,max(nN),1)) = nds;
ngh = ngh';
function [ngh,ngh2] = bks_nd2tri2ndO(tri,nd2tri,Ibnd,edg)
% SLOWER ALTERNATIVE FOR ORDERING NGH
tris = size(tri,1);
nds = size(nd2tri,1);
bnds = sum(Ibnd);
Inds = not(Ibnd);
maxnghI = max(sum(nd2tri > 0,2));
maxnghb = max(sum(nd2tri > 0,2)+Ibnd); maxngh = max(maxnghb,maxnghI);
tri_ = [tri; repmat(0,1,size(tri,2))];
trin = nd2tri; trin(trin == 0) = tris+1;
trin = [trin repmat(tris+1,nds,maxnghb-maxnghI)]; %pad tris+1 when maxnghb>maxnghI
triN = sum(trin~=tris+1,2)+Ibnd;
trind = [tri_(trin',1) tri_(trin',2) tri_(trin',3)];
empty = (trin==tris+1)'; 
if size(tri,2) == 4
	trind = [trind tri_(trin',4)];
end;
if nargin==4
cnnds1 = reshape(repmat(edg(:,1)',maxngh*4,1),4,nds*maxngh)';
cnnds2 = reshape(repmat(edg(:,2)',maxngh*4,1),4,nds*maxngh)';
I1 = cnnds1 ~= trind; I1(empty(:),1) = 0;
I2 = cnnds2 ~= trind; I2(empty(:),2) = 0;
I = and(I1,I2);
trind = trind'; trind = reshape(trind(I')',2,size(trind,2))';
I = any([and(I1(:,2) == 0,I2(:,1) == 0) and(I1(:,3) == 0,I2(:,2) == 0) and(I1(:,4) == 0,I2(:,3) == 0) and(I1(:,1) == 0,I2(:,3) == 0) and(I1(:,2) == 0,I2(:,4) == 0) and(I1(:,4) == 0,I2(:,1) == 0)],2);
else %not swap
dimp = size(tri,2);
cnnds = reshape(repmat(1:nds,maxngh*dimp,1),dimp,nds*maxngh)';
I = cnnds ~= trind; I(empty(:),1) = 0;
trind = trind'; trind = reshape(trind(I')',dimp-1,size(trind,2))';
if size(tri,2) == 3
I = I(:,2) == 1;
else %3D
I = or(I(:,2) == 0, I(:,4) == 0);
end;
end;

trind_ = trind(I,1); trind(I,1) = trind(I,2); trind(I,2) = trind_;
if nargout == 2 %size(tri,2) == 4 && nargin~=4 %
ngh=trind; 
ngh2 = reshape(ngh',3,maxngh,nds);
%ngh = reshape(ngh',3,maxngh,nds); ngh2 = mk_tcrcls(ngh(:,:,Ibnd)); return;
[trind,nds,maxngh,triN] = mk_tcrcls(ngh2(:,:,Ibnd));
ngh = zeros(nds,maxngh);
ngh(:,[1 2]) = trind(1:maxngh:end,:);
else
ngh = zeros(nds,maxngh);
ngh(Inds,[1 2]) = trind(reshape([Inds'; false(maxngh-1,nds)],nds*maxngh,1),:);
%find starting boundary node
Ibnd_ = repmat(Ibnd',maxngh,1);
trind_ = reshape(trind(Ibnd_(:),:)',2*maxngh,bnds);
truss_ = repmat([true;false],maxngh,bnds); C = repmat(1:bnds,2*maxngh,1);
[trind_,R] = sort(trind_); truss=zeros(size(R)); truss(:) = truss_(R(:)+(C(:)-1)*2*maxngh);
I = and([trind_(1,:)~=trind_(2,:); and(trind_(2:end-1,:)~=trind_(3:end,:),trind_(2:end-1,:)~=trind_(1:end-2,:)); trind_(end,:)~=trind_(end-1,:)],truss);
thmsk = zeros(2*maxngh,1); thmsk(1:2:2*maxngh) = 1:maxngh;
I_ = false(maxngh,bnds); I_(thmsk(R(I))+(0:bnds-1)'*maxngh) = true;
trind_ = trind(Ibnd_(:),:);
ngh(Ibnd,[1 2]) = trind_(I_(:),:);
end;
R=(1:nds)';

for i=3:maxngh
    C = min(i-1,triN-1);
    ndsrch = repmat(ngh(R+(C-1)*nds)',maxngh,1);
    ngh(R,i) = trind(ndsrch(:) == trind(:,1),2).*(triN>=i);
end;

function [trind,nds,maxngh,Nn] = mk_tcrcls(circles)
%save for_debug3Da.mat;
nds = size(circles,3);
nds_ = size(circles,2);
nd1 = reshape(repmat(1:size(circles,3),3*size(circles,2),1),size(circles));
nd2 = reshape(repmat((1:3)',1,size(circles,2)*size(circles,3)),size(circles));
circles_ = reshape(circles,3*size(circles,2),size(circles,3));
edga2 = circles_(reshape([(2:3:3*nds_); (3:3:3*nds_); (1:3:3*nds_)],3*nds_,1),:);
[edga,I2] = sort([circles_(:) edga2(:)],2);
[edga,I] = sortrows([reshape(nd1,3*nds*nds_,1) edga]);
%d = and(edga(:,2)~=0,any([edga(1,:) ~= edga(2,:); and(edga(2:end-1,:) ~= edga(3:end,:),edga(2:end-1,:) ~= edga(1:end-2,:)); edga(end,:) ~= edga(end-1,:)],2));
d = and(edga(:,2)~=0,[any(edga(1,:) ~= edga(2,:),2); and(any(edga(2:end-1,:) ~= edga(3:end,:),2),any(edga(2:end-1,:) ~= edga(1:end-2,:),2)); any(edga(end,:) ~= edga(end-1,:),2)]);
%save for_debug3D.mat; error('just stop3D');
edga = edga(d,:);
Iflp = I2(I(d),1) == 2;
edga(Iflp,:) = edga(Iflp,[1 3 2]);
Nn = diff([0; find([edga(1:end-1,1)~=edga(2:end,1); true])]); maxngh = max(Nn);
R = ones(size(edga,1),1); R(1+cumsum(Nn(1:end-1))) = 1+maxngh-Nn(1:end-1);
R = cumsum(R);
trind = zeros(nds*maxngh,2);
trind(R,:) = edga(:,[2 3]);


function tri2edg2tri = bks_tri2edg2tri(tri)
if size(tri,2) == 3
 [edg,edg2tri,tri2edg] = bks_all(tri);
 tri2tri = [edg2tri(tri2edg(:,1),:)'; edg2tri(tri2edg(:,2),:)'; edg2tri(tri2edg(:,3),:)'];
 tri2edg2tri = reshape(tri2tri(tri2tri~=repmat(1:size(tri,1),6,1)),3,size(tri,1))';
else
 [fac,fac2tri,tri2fac] = bks_all3D(tri);
 tri2tri = [fac2tri(tri2fac(:,1),:)'; fac2tri(tri2fac(:,2),:)'; fac2tri(tri2fac(:,3),:)'; fac2tri(tri2fac(:,4),:)'];
 tri2edg2tri = reshape(tri2tri(tri2tri~=repmat(1:size(tri,1),8,1)),4,size(tri,1))';
end;


function triO = elem_find(bndmesh,I,xy_,options)
%This function finds element numbers (triO) containing a set of coordinates (xy_) using a good 
%starting element number, bndmesh.trinew(I). If the initial guess is wrong, the guess is moved to
%one of the neighbough elements
%%% INPUT %%%
%bndmesh.triold  : element list, N x 3 or N x 4 (2D or 3D), where N is the number of elements
%bndmesh.xyold   : node coordinates, M x 2 or M x 3 (2D or 3D), where M is the number of nodes
%bndmesh.ngh     : element neighboughs, N x 3 or N x4 (2D or 3D)
%bndmesh.trinew  : old element numbers containing new node coordinates, P x 1, where P is the number of new nodes
%I               : list of new node number where bndmesh.trinew might be wrong
%xy_             : list of coordinates for said nodes
%options         : struct containing options.geomtol which is used as tolerance

%%% OUTPUT %%%
%triO            : old element numbers containing xy_ 
if size(xy_,2) == 2
	ind = [3 1 2]'; 
else
	ind = [4 1 2 3]';
end;
tri = bndmesh.trinew(I);
triO = zeros(size(xy_,1),1);
R = (1:size(xy_,1))';
while true
[tmp,triO_,s,Igood] = elem_interp(bndmesh.triold,bndmesh.xyold,tri,xy_,options); % sum(not(Igood))
nI = not(Igood); RC = triO_+(ind(s)-1)*size(bndmesh.triold,1);
II = bndmesh.ngh(RC(nI))==0;
if any(II) && options.verbose ~= -1
   If1 = find(nI); xye = xy_(If1(find(II,1)),:);
   warning(sprintf('extrapolating to %1.1e %1.1e %1.1e',xye));
end;	
Igood(nI) = II; nI(nI) = not(II); %extrapolation
triO(R(Igood)) = triO_(Igood);
if all(Igood)
	break
end;
%I(Igood) = false;
%nI = not(Igood); RC = triO_(nI)+(ind(s(nI))-1)*size(bndmesh.triold,1);
tri = bndmesh.ngh(RC(nI)); %tri(nI);
xy_ = xy_(nI,:);
R = R(nI);
end;

function [Nmetric,triO,s,Igood] = elem_interp(tri,xy,nd2tri,xy_,options,Nmetric)
%the function finds the elements (triO) containing the input coordinates (xy_)
%%% INPUT %%%
%tri        : element list, N x 3 or N x 4 (2D or 3D), where N is the number of elements
%xy         : node coordinates, M x 2 or M x 3 (2D or 3D), where M is the number of nodes
%nd2tri     : masked array of potential elements to check, Q x P, where Q is the number of elements 
%             to find and P is the maximum number of potential elements. Code is somewhat simplified for P==1.
%xy_        : coordinates inside the desired elements, Q x 2 or Q x 3 (2D or 3D)
%options    : struct containing the tolerance options.geomtol, which controls how far outside the coordinate can be
%Nmetric    : optional input to be interpolated, M x R, where R is the number of fields to be interpolated

%%% OUTPUT %%%
%Nmetric    : interpolated fields, Q x R, if input is not given [] is returned
%triO       : the best candidate elements containing xy_
%s          : the local edge number closest to xy_, Q x 1. Relevant for iterating towards the correct 
%             element (see elem_find)
%Igood      : boolean Q x 1. True for correct elements (i.e. no extrapolation).

dim = size(xy,2);
pts = size(xy_,1); dimpts = dim*pts;
if size(nd2tri,2) == 1
 nd2tri_ = nd2tri;
 R = (1:size(nd2tri,1))'; C = ones(size(R)); I = true(size(R)); Isumdim = dim*numel(R);
else
NN = sum(nd2tri~=0,2); maxNN = max(NN);
xy_ = repmat(xy_',[1,1,maxNN]); I = repmat(1:maxNN,pts,1) <= repmat(NN,1,maxNN); xy_ = xy_(:,I)';
[R,C] = find(nd2tri~=0);
I = nd2tri(:) ~= 0; nd2tri_ = reshape(nd2tri(I),sum(I),1); Isum = size(nd2tri_,1); Isumdim = Isum*dim;
end;
t1 = xy(tri(nd2tri_,2),:) - xy(tri(nd2tri_,1),:);
t2 = xy(tri(nd2tri_,3),:) - xy(tri(nd2tri_,1),:);
%C = repmat(1:pts*dim,dim,1);
if dim==2
 digs = [-1 0 1]; All = zeros(Isumdim,3);
 All(1:2:Isumdim,2) = t1(:,1);
 All(1:2:Isumdim,1) = t1(:,2);
 All(2:2:Isumdim,3) = t2(:,1);
 All(2:2:Isumdim,2) = t2(:,2);

 A = spdiags(All,digs,Isumdim,Isumdim);
 X = A\reshape((xy_-xy(tri(nd2tri_,1),:))',Isumdim,1);
 s1 = X(1:2:Isumdim); s2 = X(2:2:Isumdim); s3 = 1-s1-s2;
 spike = zeros(size(s1));
 I_ = and(s1 <  s2, s1 <= s3); spike(I_) = -s1(I_);
 I_ = and(s2 <= s1, s2 <  s3); spike(I_) = -s2(I_);
 I_ = and(s3 <= s2, s3 <  s1); spike(I_) = -s3(I_);
 spikes = inf(size(nd2tri)); spikes(I) = spike;
 [tmp,Cbest] = sort(spikes,2);
 I_ = Cbest(R)==C; [Rbest,I] = sort(R(I_));
 triO = nd2tri(Rbest+(Cbest(:,1)-1)*pts);
 s1 = s1(I_); s1 = s1(I);
 s2 = s2(I_); s2 = s2(I);
 s3 = s3(I_); s3 = s3(I);
 s = [s1 s2 s3];
 if  nargin == 6 && any(options.geomtol < tmp(:,1)) && options.verbose ~= -1
   xye = xy_(I_,:); xye = xye(I,:); xye = xye(find(options.geomtol < tmp(:,1),1),:); 
   warning(sprintf('extrapolating to (%1.1e,%1.1e)',xye));
 end;
 if nargin == 6
 Nmetric = Nmetric(tri(triO,1),:).*repmat(s3,1,size(Nmetric,2)) ...
         + Nmetric(tri(triO,2),:).*repmat(s1,1,size(Nmetric,2)) ...
         + Nmetric(tri(triO,3),:).*repmat(s2,1,size(Nmetric,2));
 else
 Nmetric = [];
 end;
else %3D
 t3 = xy(tri(nd2tri_,4),:) - xy(tri(nd2tri_,1),:);
 digs = -2:2;  All = zeros(Isumdim,5);
 All(1:3:Isumdim,3) = t1(:,1);
 All(1:3:Isumdim,2) = t1(:,2);
 All(1:3:Isumdim,1) = t1(:,3);
 All(2:3:Isumdim,4) = t2(:,1);
 All(2:3:Isumdim,3) = t2(:,2);
 All(2:3:Isumdim,2) = t2(:,3);
 All(3:3:Isumdim,5) = t3(:,1);
 All(3:3:Isumdim,4) = t3(:,2);
 All(3:3:Isumdim,3) = t3(:,3);
 A = spdiags(All,digs,Isumdim,Isumdim);
 X = A\reshape((xy_-xy(tri(nd2tri_,1),:))',Isumdim,1);
 s1 = X(1:3:Isumdim); s2 = X(2:3:Isumdim); s3 = X(3:3:Isumdim); s4 = 1-s1-s2-s3;
 spike = inf(size(s1));
 I_ = and(and(s1 <= s2, s1 <  s3), s1 <  s4); spike(I_) = -s1(I_);
 I_ = and(and(s2 <  s1, s2 <  s3), s2 <= s4); spike(I_) = -s2(I_);
 I_ = and(and(s3 <= s1, s3 <= s2), s3 <  s4); spike(I_) = -s3(I_);
 I_ = and(and(s4 <= s1, s4 <  s2), s4 <= s3); spike(I_) = -s4(I_);
 spikes = inf(size(nd2tri)); spikes(I) = spike;
 [tmp,Cbest] = sort(spikes,2);
 I_ = Cbest(R)==C; [Rbest,I] = sort(R(I_));
 triO = nd2tri(Rbest+(Cbest(:,1)-1)*pts);
 s1 = s1(I_); s1 = s1(I); s2 = s2(I_); s2 = s2(I);
 s3 = s3(I_); s3 = s3(I); s4 = s4(I_); s4 = s4(I);
 s = [s1 s2 s3 s4];
 if nargin == 6 && any(options.geomtol < tmp(:,1)) && options.verbose ~= -1
   xye = xy_(I_,:); xye = xye(I,:); xye = xye(find(options.geomtol < tmp(:,1),1),:); 
   warning(sprintf('extrapolating to (%1.1e,%1.1e,%1.1e)',xye));
 end;
 if nargin == 6
 Nmetric = Nmetric(tri(triO,1),:).*repmat(s4,1,size(Nmetric,2)) ...
         + Nmetric(tri(triO,2),:).*repmat(s1,1,size(Nmetric,2)) ...
         + Nmetric(tri(triO,3),:).*repmat(s2,1,size(Nmetric,2)) ...
         + Nmetric(tri(triO,4),:).*repmat(s3,1,size(Nmetric,2));
 else
 Nmetric = [];
 end;
end;
if nargout > 2
  [tmp,s] = min(s,[],2); Igood = -options.geomtol < tmp;
  %ind = [3 1 2]; ind3D = [4 1 2 3]; 
end;


function [quality,c,wght1] = elem_qual(tri,xy,Nmetric,options,nvec)
%options.qualM == 1 Vassilevski
%options.qualM == 2 Orth, LW(H), worst
%options.qualM == 3 Orth, LW(H), worst vassilevski
%options.qualM == 4 Orth, ellipse, worst
%options.qualM == 5 Orth, ellipse, worst vassilevski
%options.qualM == 6 Orth, LW(H), alt worst
%options.qualM == 7 Orth, ellipse, alt worst
%options.qualM == 8 Home brew condition (sliver functional), now it is shephard
%options.qualM == 9 Condition number
%options.qualM == 10 volume
%options.qualM == 11 frenchy
%options.qualM == 12 Pascal
%options.qualM == 13 Pascal sq (real)
%opions.qualP   > 0  angle functional



if nargin == 5 && size(tri,2) ~= 2 %we wont bother calculating quality of inverted elements
  quality = repmat(-1,size(tri,1),1);
  if size(nvec,2) == 3
    [Ibad,area] = elem_inv(tri,xy,nvec); 
  else
    [Ibad,area] = elem_inv(tri,xy); 
  end;
  if size(tri,2) == 4 || size(tri,2) == 3
  	Ibad = area < options.minA;
  end;
  if any(not(Ibad))
  quality(not(Ibad)) = elem_qual_slow(tri(not(Ibad),:),xy,Nmetric,options);
  end;
else
  if nargout == 3
	[quality,c,wght1] = elem_qual_slow(tri,xy,Nmetric,options);
  elseif nargout == 2
  	[quality,c] = elem_qual_slow(tri,xy,Nmetric,options);
  else
  	quality = elem_qual_slow(tri,xy,Nmetric,options);
  end;
end;
	
function [quality,c,wght1] = elem_qual_slow(tri,xy,Nmetric,options)	
if size(tri,2) == 2
	if nargout == 3
	[quality,c,wght1] = gen_rel_edgL(tri,xy,Nmetric,options);
	elseif nargout == 2
	[quality,c] = gen_rel_edgL(tri,xy,Nmetric,options);
	else
	quality = gen_rel_edgL(tri,xy,Nmetric,options);
	end;
else
        if options.qualP > 0
  	  quality = elem_angle(tri,xy,options);
  	  return;
         end;
if options.qualM == 1 || options.qualM > 7 
quality = vassilevski(tri,xy,Nmetric,options);
end;
if options.qualM == 2 || options.qualM == 3 || options.qualM == 6
    quality = worst_HL(tri,xy,Nmetric,options);
end;
if options.qualM == 4 || options.qualM == 5 || options.qualM == 7
     quality = steiner_ell_metric(tri,xy,Nmetric,options);
end;
end;


function quality = vassilevski(tri,xy,Nmetric,options)
if size(tri,2) == 3
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2] = tri_in_metric_space(tri,xy,Nmetric,options);
	L = sqrt([sum(v1.^2,2) sum(v2.^2,2) sum((v1-v2).^2,2)]);
else
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2, xy4I, v3] = tri_in_metric_space(tri,xy,Nmetric,options);
	L = sqrt([sum(v1.^2,2) sum(v2.^2,2) sum(v3.^2,2) sum((v1-v2).^2,2) sum((v1-v3).^2,2) sum((v2-v3).^2,2)]);
end;
if options.qualM == 10
quality = area; return;
end;
if size(tri,2) == 3
if options.qualM == 1
quality = area./mean(L,2).^2.*myf(mean(L,2)); 
elseif options.qualM == 12
quality = area./mean(L,2).^2; 
elseif options.qualM == 13
quality = sign(area).*(area./mean(L,2).^2).^2;
elseif options.qualM == 9 || options.qualM == 8 || options.qualM == 11 
quality = area./mean(L.^2,2);
end;
else %3D
	[areaA1,H1] = gen_crosspro(xy1I,xy2I,xy3I,xy4I);
	[areaA2,H2] = gen_crosspro(xy4I,xy2I,xy3I,xy1I);
	[areaA3,H3] = gen_crosspro(xy1I,xy4I,xy3I,xy2I);
        [areaA4,H4] = gen_crosspro(xy1I,xy2I,xy4I,xy3I);
        L2 = [areaA1 areaA2 areaA3 areaA4];
        if options.qualM == 1
        quality = area./mean(L,2).^3.*myf(mean(L,2));
        %quality = area./mean(L,2).^3.*myf(mean(L2,2));
        %quality = area./mean(L.^2,2).^(3/2);
        %quality = area./mean(L2,2).^(3/2).*myf(mean(L,2));
        %quality = area./mean(L2,2).^(3/2).*mean(myf(L),2);
        elseif options.qualM == 12
         quality = area./mean(L,2).^3;
        elseif options.qualM == 13
         quality = sign(area).*(area./mean(L,2).^3).^2;
        elseif options.qualM == 11
         quality = area./mean(L.^2,2).^(3/2); %frenchy
        elseif options.qualM == 8
       	%quality = area./mean(L,2).^3; 
       	%quality = area./mean(L,2)./mean(L2,2); 
       	quality = area.^2./mean(L.^2,2).^3; 
        else %options.qualM == 9
        %quality = area./mean(L,2).^3.*myf(mean(L2,2)); 
        quality = area./sqrt(mean(L.^2,2).*mean(L2.^2,2));
        end;
	%qual2Ds = elem_qual([tri(:,1:3); tri(:,2:4); tri(:,[1 3 4]); tri(:,[1 2 4])],xy,Nmetric,options); qual2Ds = reshape(qual2Ds,size(tri,1),4);
	%quality = area./mean(L2,2).^(3/2).*mean(qual2Ds,2);
end;


function quality = steiner_ell_metric(tri,xy,Nmetric,options)
if size(tri,2) == 3
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2] = tri_in_metric_space(tri,xy,Nmetric,options);
	XI = elem_metric(tri,xy,xy1I,xy2I,xy3I);

else
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2, xy4I, v3] = tri_in_metric_space(tri,xy,Nmetric,options);
	XI = elem_metric([],[],xy1I,xy2I,xy3I,xy4I);
end;
[eigL,eigR] = analyt_eig(XI); 
if any(eigL(:)<0)
    warning('neg eig'); eigL = abs(eigL);
end;
eigL = 1./sqrt(eigL);
if options.qualM == 4
	quality = min(myfa(eigL),[],2).*sign(area);
elseif options.qualM == 7
	quality = min(relE(eigL),[],2).*sign(area);
else
	quality = area./mean(eigL,2).^(size(tri,2)-1).*myf(mean(eigL,2));
end;


function quality = worst_HL(tri,xy,Nmetric,options)
if size(tri,2) == 3
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2] = tri_in_metric_space(tri,xy,Nmetric,options);
	L = sqrt([sum(v1.^2,2) sum(v2.^2,2) sum((v1-v2).^2,2)]);	
	H = repmat(abs(area),1,3)./L;
	%2D unit element has height sqrt(3)/2
	vals = [L(:,1) H(:,1) L(:,2) H(:,2) L(:,3) H(:,3)];
else
	[Ibad, area, xy1I, xy2I, xy3I, v1, v2, xy4I, v3] = tri_in_metric_space(tri,xy,Nmetric,options);
	L = sqrt([sum(v1.^2,2) sum(v2.^2,2) sum(v3.^2,2) sum((v1-v2).^2,2) sum((v1-v3).^2,2) sum((v2-v3).^2,2)]); %(12,13,14,23,24,34)
	[areaA1,H1] = gen_crosspro(xy1I,xy2I,xy3I,xy4I);
	[areaA2,H2] = gen_crosspro(xy4I,xy2I,xy3I,xy1I);
	[areaA3,H3] = gen_crosspro(xy1I,xy4I,xy3I,xy2I);
        [areaA4,H4] = gen_crosspro(xy1I,xy2I,xy4I,xy3I);
	Wa = [areaA1 areaA1 areaA3 areaA1 areaA2 areaA2]./L;
	Wb = [areaA4 areaA3 areaA4 areaA2 areaA4 areaA3]./L;
	Ha = [H1 H1 H3 H1 H2 H2]; %Ha = repmat(area,1,6)./(L.*Wa);
	Hb = [H4 H3 H4 H2 H4 H3]; %Hb = repmat(area,1,6)./(L.*Wb);
	vals = [L Wa Wb Ha Hb];
end;
if options.qualM == 2 || options.qualM == 6
	nzero = area~=0;
	if options.qualM == 2 
	worst_ = min(myfa(vals(nzero,:)),[],2);
	else
	worst_ = min(relE(vals(nzero,:)),[],2);
	end;
	worst = zeros(size(L,1),1); worst(nzero) = worst_;
	quality = worst.*sign(area);
else
	myf2 = @(x_) area./mean(x_,2).^(size(tri,2)-1).*myf(mean(x_,2));
	if size(tri,2) == 3
		qual = [myf2([L(:,1) H(:,1)]) myf2([L(:,2) H(:,2)]) myf2([L(:,3) H(:,3)])];
	else
		qual = [myf2([L(:,1) Wa(:,1) Ha(:,1)]) myf2([L(:,1) Wb(:,1) Hb(:,1)]) ...
         		        myf2([L(:,2) Wa(:,2) Ha(:,2)]) myf2([L(:,2) Wb(:,2) Hb(:,2)]) ...
         		        myf2([L(:,3) Wa(:,3) Ha(:,3)]) myf2([L(:,3) Wb(:,3) Hb(:,3)]) ...
         		        myf2([L(:,4) Wa(:,4) Ha(:,4)]) myf2([L(:,4) Wb(:,4) Hb(:,4)]) ...         		        
         		        myf2([L(:,5) Wa(:,5) Ha(:,5)]) myf2([L(:,5) Wb(:,5) Hb(:,5)]) ...
         		        myf2([L(:,6) Wa(:,6) Ha(:,6)]) myf2([L(:,6) Wb(:,6) Hb(:,6)])];
	end;
	quality = min(qual,[],2);%.*sign(area);
end;

function [mm,xy1,xy2,xy3,xy4] = calc_tri_metric(tri,xy,Nmetric,options);
xy1 = xy(tri(:,1),:); xy2 = xy(tri(:,2),:); xy3 = xy(tri(:,3),:); 
if size(tri,2) == 4
	xy4 = xy(tri(:,4),:);
end;
mm = metric_avg(tri,Nmetric,options);

function out = myf(x)
out = (min(x,1./x).*(2-min(x,1./x))).^3;

function out = myfa(myx)
out= min(myx,1./myx);

function out = relE(myx)
out = (myx.*(myx<1)+(2-myx).*(myx>=1));

function [crosspro,H] = gen_crosspro(x1,x2,x3,x4)
v1 = x2-x1;
v2 = x3-x1;
v3 = x4 - (x1+x2+x3)/3.;
H = sqrt(sum(v3.^2,2))/sqrt(2/3); %3D unit element has height sqrt(2/3)
crosspro = cross(v1,v2,2); %[v1(:,2).*v2(:,3)-v1(:,3).*v2(:,2) ...
%		   v1(:,3).*v2(:,1)-v1(:,1).*v2(:,3) ...
%		   v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)];
crosspro = sqrt(sum(crosspro.^2,2))*0.5/(sqrt(3)/4);

function [edgnowl,c,wght1] = gen_rel_edgL(edg,xy,Nmetric,options)
mmI = metric_avg(edg,Nmetric,options);
x = xy(:,1); y = xy(:,2);
x_ = reshape(x(edg),size(edg)); y_ = reshape(y(edg),size(edg));
	
if size(xy,2) == 3 %3D
	z = xy(:,3); z_ = reshape(z(edg),size(edg));
	xy1_ = [x_(:,1) y_(:,1) z_(:,1)]; xy2_ = [x_(:,2) y_(:,2) z_(:,2)];
else
	xy1_ = [x_(:,1) y_(:,1)]; xy2_ = [x_(:,2) y_(:,2)];
end;

%xyI1 = analyt_prod(xy1_,mmI); 
%xyI2 = analyt_prod(xy2_,mmI);
%t = xyI1-xyI2;
%edgnowl = sqrt(sum(t.^2,2)); %edge length
t = xy1_-xy2_;
edgnowl = sqrt(sum(analyt_prod(t,mmI).^2,2)); %edge length

if nargout > 1
	if options.spltc == 1
		c = (xy1_ + xy2_)/2; wght1 = repmat(0.5,size(c,1),1);
	else
	%t = xy1_-xy2_; %t = t./repmat(sqrt(sum(t.^2,2)),1,size(t,2)); %compute invLs
	mmI1 = Nmetric(edg(:,1),:); L1 = sqrt(sum(analyt_prod(t,mmI1).^2,2));
	mmI2 = Nmetric(edg(:,2),:); L2 = sqrt(sum(analyt_prod(t,mmI2).^2,2));
	wght1 = L1./(L1+L2); wght1_ = repmat(wght1,1,size(xy1_,2));
	c = xy1_.*wght1_+xy2_.*(1.-wght1_);
	end;
end;

function [Ibad, area, xy1I, xy2I, xy3I, v1, v2, xy4I, v3] = tri_in_metric_space(tri,xy,Nmetric,options)
if size(tri,2) == 3
	[mm,xy1,xy2,xy3] = calc_tri_metric(tri,xy,Nmetric,options);
else
	[mm,xy1,xy2,xy3,xy4] = calc_tri_metric(tri,xy,Nmetric,options);
end;
xy1I = analyt_prod(xy1,mm);
xy2I = analyt_prod(xy2,mm);
xy3I = analyt_prod(xy3,mm);
v1 = xy2I - xy1I; 
v2 = xy3I - xy1I;
if size(tri,2) == 4
	xy4I = analyt_prod(xy4,mm);
	v3 = xy4I - xy1I;
	[Ibad,area] = elem_inv([],[],xy1I,xy2I,xy3I,xy4I);
	area = (1/6)/(sqrt(2)/12)*area; %volume
else
	[Ibad,area] = elem_inv([],[],xy1I,xy2I,xy3I);
	area = 0.5/(sqrt(3)/4)*area;
end;
	

function [badtri,z,zN] = elem_inv(tri,xy,xy1,xy2,xy3,xy4)
if (nargin == 2 && numel(tri) == 0) || (nargin ~= 2 && numel(xy1) == 0) 
    badtri = []; z = []; 
    return;
end;
if nargin == 3 && size(tri,2) == 3 && size(xy,2) == 3
	nvec = xy1;
end;
% orient for check of inversion later on:
if nargin <= 3
	xy1 = xy(tri(:,1),:);
	xy2 = xy(tri(:,2),:);
	xy3 = xy(tri(:,3),:);
	if size(tri,2) == 4
		xy4 = xy(tri(:,4),:);
	end;
end;
v1 = xy2-xy1;
v2 = xy3-xy1;
if size(xy1,2) == 2
	z = v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1);
else
	crosspro = cross(v1,v2,2); %[v1(:,2).*v2(:,3)-v1(:,3).*v2(:,2) ...
%			   v1(:,3).*v2(:,1)-v1(:,1).*v2(:,3) ...
%			   v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)];
	if (size(tri,2) == 3 || nargin == 5) && size(xy1,2) == 3 
		if nargin == 3
		z = sum(nvec.*crosspro,2); %exterior face taking normal into account
		else
		z = sqrt(sum(crosspro.^2,2)); %interior face only requires area calculation
		end;
	else
		v3 = xy4 - (xy1+xy2+xy3)/3;
	         z = sum(v3.*crosspro,2);
        end;
end;
badtri = z<0; %counter-clock-wise
if nargout == 3
 zN = full(sparse(tri(:),ones(numel(tri),1),repmat(z,size(tri,2),1)))/(6+18*(size(xy,2)==3));
end;

function gammaf = fem_filter(tri,xy,bndmesh,Lmin,gamma,bcs,options,X)
% solves the PDE (Lmin**2*inner(grad(v),grad(u))+u*v)*dx == gamma*v*dx (FEniCS notation), 
% where v and u are linear test and trial functions (FunctionSpace(mesh,'CG',1))
%%% INPUT %%%
%tri        : element list, N x 3 or N x 4 (2D or 3D), where N is the number of elements
%xy         : node coordinates, M x 2 or M x 3 (2D or 3D), where M is the number of nodes
%bndmesh.edg: boundary edge list (2D), P x 2, where P is the number of boundary edges 
%bndmesh.fac: boundary face list (3D), P x 3, where P is the number of boundary faces
%bndmesh.IDs: boundary ID list. 
%            The boundary mesh is used to impose Dirchlet boundary conditions, so it is optional.
%Lmin       : filter length
%gamma      : fields to be filtered N x Q, where Q is the number of fields
%bcs        : optional variable with the ID of the boundary on which to impose Dirichlet boundary conditions
%options    : struct containing options for the linear solve (\ or cgs)
%X          : output of fem_getX, which can be used to calculate derivatives of on linear elements

%%% OUTPUT %%%
%gammaf     : nodal filtered fields, M x Q
dim = size(xy,2);
LL = Lmin^2;
gammaf = zeros(size(xy,1),size(gamma,2));
if nargin == 7
  [X,dem_] = fem_getX(tri,xy);
else
  [I,dem_] = elem_inv(tri,xy);
end; 
facts = [1./3.  1./12. 1./60.; 1./6. 1./24. 1./120.; 0.5 1./6. 1./24.; 1. 0.5 1./6.];
Io = {[2 1], [2 3 1 3 1 2],[2 3 4 1 3 4 1 2 4 1 2 3]};
Rd1 = tri(:); Cd1 = Rd1; Sd1 = repmat(dem_,size(tri,2),1)*facts(1,dim);
Co1 = repmat(tri(:),dim,1); Ro1 = reshape(tri(:,Io{dim}),size(Co1)); So1 = repmat(dem_,(dim+1)*dim,1)*facts(2,dim);
%K = diag(repmat(LL,size(xy,2)));
K = diag(repmat(LL,dim,1)); i1 = -3;
C2 = zeros(size(tri,1),4*dim*dim*sum(K(:)~=0)); R2 = C2; S2 = C2;
for hh=1:dim %grad(g_test,x[hh])
for ii=1:dim %grad(g,x[ii]) 
if K(hh,ii) == 0.
	continue;
end;
for jj=1:dim %innerit g_test
for mm=1:dim %innerit g
i1 = i1+4;
i2 = i1+1; i3 = i1+2; i4 = i1+3;
X_ = X(:,dim*(jj-1)+hh).*X(:,dim*(mm-1)+ii).*dem_*facts(4,dim)*K(hh,ii);
R2(:,i1) = tri(:,1);
C2(:,i1) = tri(:,1);
S2(:,i1) = X_;
R2(:,i2) = tri(:,1);
C2(:,i2) = tri(:,jj+1);
S2(:,i2) = -X_;
R2(:,i3) = tri(:,mm+1);
C2(:,i3) = tri(:,1);
S2(:,i3) = -X_;
R2(:,i4) = tri(:,mm+1);
C2(:,i4) = tri(:,jj+1);
S2(:,i4) = X_;
end;
end;
end;
end;

R = [Rd1; Ro1; R2(:)]; C = [Cd1; Co1; C2(:)]; S = [Sd1; So1; S2(:)];
%set bcs
zerodofs = false(size(xy,1),1);
if numel(bcs) ~= 0
for ii=1:numel(bcs{1})
I = bndmesh.IDs == bcs{1}(ii);
if dim==2
 zerodofs(bndmesh.edg(I,:)) = true;
else
 zerodofs(bndmesh.fac(I,:)) = true;
end;
S(and(R ~= C,or(zerodofs(R),zerodofs(C)))) = 0;
S(and(R == C,zerodofs(R))) = 1.; 
end;%for
end;%if
A = sparse(R(:),C(:),S(:),size(xy,1),size(xy,1));

for ii=1:size(gamma,2)
if size(gamma,1) == size(tri,1)
R = tri(:); C = repmat(1,numel(tri),1); S = repmat(dem_.*gamma(:,ii),dim+1,1)*facts(3,dim);	
b = sparse(R,C,S);
elseif size(gamma,1) == size(tri,1)+size(xy,1)
error('not implemented');
else
Rd = tri(:); Cd = repmat(1,numel(tri),1); Sd = repmat(dem_,dim+1,1).*gamma(tri(:))*facts(1,dim);
Co = repmat(1,size(tri,1)*(dim+1)*dim,1); Ro = repmat(tri(:),dim,1); 
So = repmat(dem_,(dim+1)*dim,1).*gamma(reshape(tri(:,Io{dim}),size(Co)),ii)*facts(2,dim);
%So = repmat(dem_,(dim+1)*dim,1).*gamma(Ro,ii)*facts(2,dim); Ro = reshape(tri(:,Io{dim}),size(Co));
b = sparse([Rd; Ro], [Cd; Co], [Sd; So]);
end;
b(zerodofs) = 0;
gammaf(:,ii) = fem_solve(A,b,options); 
end;


function CGf = fem_tri2xy(tri,xy,DGf,options)
% projects element-wise constant fields to nodal fields, u*v*dx == DGf*v*dx (FEniCS notation), 
% where v and u are linear test and trial functions (FunctionSpace(mesh,'CG',1))
%%% INPUT %%%
%tri        : element list, N x 3 or N x 4 (2D or 3D), where N is the number of elements
%xy         : node coordinates, M x 2 or M x 3 (2D or 3D), where M is the number of nodes
%DGf        : fields to be projected N x Q, where Q is the number of fields
%options    : struct containing options for the linear solve (\ or cgs)

%%% OUTPUT %%%
%CGf        : nodal fields, M x Q
CGf = zeros(size(xy,1),size(DGf,2));
dim = size(xy,2);
%nd2tri = inv_table(tri); nd2tri_ = nd2tri(nd2tri~=0);
[I,dem_] = elem_inv(tri,xy);
facts = [1./3.  1./12. 1./60.; 1./6. 1./24. 1./120.; 0.5 1./6. 1./24.];
Io = {[2 1], [2 3 1 3 1 2],[2 3 4 1 3 4 1 2 4 1 2 3]};
Rd = tri(:); Cd = Rd; Sd = repmat(dem_,dim+1,1)*facts(1,dim);
Co = repmat(tri(:),dim,1); Ro = reshape(tri(:,Io{size(xy,2)}),size(Co)); So = repmat(dem_,(dim+1)*dim,1)*facts(2,dim);
A = sparse([Rd;Ro],[Cd;Co],[Sd;So]);
%demb = zeros(size(nd2tri)); demb(nd2tri~=0) = dem_(nd2tri_)*facts(3,dim);
%DGf_ = zeros(size(nd2tri)); 
for ii=1:size(DGf,2)
R = tri(:); C = repmat(1,numel(tri),1); S = repmat(dem_.*DGf,dim+1,1)*facts(3,dim);
b = sparse(R,C,S);
%DGf_(nd2tri~=0) = DGf(nd2tri_,ii);
%b = sum(DGf_.*demb,2);
CGf(:,ii) = fem_solve(A,b,options); %CGf(:,ii) = A\b; %linsolve(A,b);
end;


function [X,dem] = fem_getX(tri,xy)
dim = size(xy,2);
tris = size(tri,1);
Isumdim = tris*dim;
b = reshape(permute(reshape(repmat([1; repmat([zeros(dim,1); 1],dim-1,1)],1,tris),[dim dim tris]),[1 3 2]),tris*dim,dim);
if nargout == 2
[I,dem] = elem_inv(tri,xy);
end;
digs = -dim+1:dim-1;
All = zeros(Isumdim,2*dim-1);
for ii=1:dim
for jj=1:dim
 All(ii:dim:end,dim-jj+ii) = xy(tri(:,1+ii),jj)-xy(tri(:,1),jj);
end;
end;
A = spdiags(All,digs,Isumdim,Isumdim);
X = A\b;
X = reshape(X',dim^2,tris)';


function A = fem_sq_grad(tri,xy,b,X)
A = zeros(size(tri,1),size(b,2));
dim = size(xy,2);
if nargin == 3
  X = fem_getX(tri,xy);
end;
for i=1:size(b,2);
triz = reshape(b(tri(:,2:end),i),size(tri)-[0 1]) - repmat(b(tri(:,1),i),1,dim);
for j=1:dim %grad(b,x[j])
for m=1:dim
for n=1:dim
A(:,i) = A(:,i) + triz(:,m).*X(:,dim*(m-1)+j).*triz(:,n).*X(:,dim*(n-1)+j);
end;
end;
end;
end;

function T = fem_heat(tri,xy,bndmesh,gamma,Emin,simpP,bcs,options,X)
%solves kappa*inner(grad(T),grad(u))*dx == u*dx (FEniCS notation), 
% where u and T are test and trial functions (FunctionSpace(mesh,'CG',1)) and
% kappa = Emin + (1-Emin)*gamma^simpP
%%% INPUT %%%
%tri        : element list, N x 3 or N x 4 (2D or 3D), where N is the number of elements
%xy         : node coordinates, M x 2 or M x 3 (2D or 3D), where M is the number of nodes
%bndmesh.edg: boundary edge list (2D), P x 2, where P is the number of boundary edges 
%bndmesh.fac: boundary face list (3D), P x 3, where P is the number of boundary faces
%bndmesh.IDs: boundary ID list. 
%gamma      : design variables, N x 1
%Emin       : minimum thermal conductivity
%simpP      : penalization exponent
%bcs        : the ID of the boundary on which to impose the Dirichlet boundary condition
%options    : struct containing options for the linear solve (\ or cgs)
%X          : output of fem_getX, which can be used to calculate derivatives of on linear elements

%%% OUTPUT %%%
%T          : nodal temperature field, M x 1
dim = size(xy,2);
T = zeros(size(xy,1),1);
if nargin == 8
  [X,dem_] = fem_getX(tri,xy);
else
  [I,dem_] = elem_inv(tri,xy);
end;
facts = [1./3.  1./12. 1./60.; 1./6. 1./24. 1./120.; 0.5 1./6. 1./24.; 1. 0.5 1./6.];
Io = {[2 1], [2 3 1 3 1 2],[2 3 4 1 3 4 1 2 4 1 2 3]};
i1 = -3;
C = zeros(size(tri,1),4*dim*dim*dim); R = C; S = C;
if size(tri,1) == numel(gamma)
 dem2 = (Emin + (1.-Emin)*gamma.^simpP).*dem_*facts(4,dim);
else
 dem2 = (Emin + (1.-Emin)*mean(gamma(tri).^simpP,2)).*dem_*facts(4,dim);
end;
for hh=1:dim %grad(g_test,x[hh]) grad(g,x[hh]) 
for jj=1:dim %innerit g_test
for mm=1:dim %innerit g
i1 = i1+4;
i2 = i1+1; i3 = i1+2; i4 = i1+3;
X_ = X(:,dim*(jj-1)+hh).*X(:,dim*(mm-1)+hh).*dem2;
R(:,i1) = tri(:,1);
C(:,i1) = tri(:,1);
S(:,i1) = X_;
R(:,i2) = tri(:,1);
C(:,i2) = tri(:,jj+1);
S(:,i2) = -X_;
R(:,i3) = tri(:,mm+1);
C(:,i3) = tri(:,1);
S(:,i3) = -X_;
R(:,i4) = tri(:,mm+1);
C(:,i4) = tri(:,jj+1);
S(:,i4) = X_;
end;
end;
end;
%A = sparse(R(:), C(:), S(:));

R2 = tri(:); C2 = repmat(1,numel(tri),1); S2 = repmat(dem_,dim+1,1)*facts(3,dim);
b = sparse(R2,C2,S2);

%set bcs
zerodofs = false(size(b));
I = bndmesh.IDs == bcs{1};
if dim==2
 zerodofs(bndmesh.edg(I,:)) = true;
else
 zerodofs(bndmesh.fac(I,:)) = true;
end;
S(and(R ~= C,or(zerodofs(R),zerodofs(C)))) = 0;
S(and(R == C,zerodofs(R))) = 1.;  A = sparse(R(:),C(:),S(:),size(xy,1),size(xy,1)); b(zerodofs) = 0.;

T(:) = A\b; %T(:) = fem_solve(A,b,options);


function X = fem_solve(A,b,options)
if options.fem_solve_it
 [X,flg] = cgs(A,b,options.fem_solve_tol,options.fem_solve_it);
 if flg
  if options.verbose ~= -1
   warning('Iterative solver did not converge; resorting to direct solver');
  end;
  X = A\b;
 end;
else
 X = A\b;
end
function [newtri,qualityN,badnds,newtriN,badbadnds,newIDs,newedg,newedgN,newfac,newfacN] = fill_circles(circles,Nmetric,xy,badnds,badIDs,triQtb,options,edg_)
NN = sum(circles~=0,2);
if size(edg_,2) == 3
circles = fix_circles(circles);
end;
if options.prag_crs && size(edg_,2) ~= 2
[newtri,qualityN,badnds,newtriN,badbadnds,newIDs,newedg,newedgN] = prag_coarse(circles,Nmetric,xy,badnds,badIDs,triQtb,NN,options,edg_); return; end;
%if options.rmnonconvec~=2
%[circles,badnds,edg_,badbadnds] = del_dblNcnvx(circles,xy,badnds,edg_,options);
%else
badbadnds = zeros(0,1);
%end;
newtri = zeros(0,3); newtriN = []; qualityN = []; badnds_start = badnds;
newedg = zeros(0,2); newedgN = []; newIDs = [];
%[newtri,circles,badnds,badbadnds,newtriN,qualityN,R2R,NN,edg_] = rm_trinds(circles,badnds,badbadnds,xy,Nmetric,triQtb,options,edg_);
%if size(circles,1) > 0
%NN = sum(circles~=0,2);
qual_ = zeros(size(circles,2),size(circles,1)); qual_(circles'~=0) = NaN;
if size(edg_,2) == 2
newtri = zeros(0,4); 
newfac = zeros(0,3);
newfacN = [];
qual3Da = qual_;
qual3Db = qual_;
elseif not(options.qualRM) 
	qual2 = qual_;
end;
allconvex = true;
while true 
    % FIND ALLOWED NEW ELEMENTS
    [Cc,R] = find(circles'); nbad = size(circles,1);
    thmp1 = [0:size(circles,2)-1]'; thmp2 = [2:size(circles,2)+1]';
    Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
    Cr = thmp2(Cc); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
    I1 = circles(R+(Cr-1)*nbad); I2 = circles(R+(Cc-1)*nbad); I3 = circles(R+(Cl-1)*nbad);
    Ipdt = isnan(qual_(Cc+(R-1)*size(qual_,1))); I123 = [I1(:) I2(:) I3(:)];
    if size(xy,2) == 2
    qual_(isnan(qual_(:))) = elem_qual(I123(Ipdt,:),xy,Nmetric,options,1);
    elseif size(edg_,2) == 2 %swapping
    %qual3Da(isnan(qual_(:))) = elem_qual([edg_(R(Ipdt),1) I123(Ipdt,[3 2 1])], xy,Nmetric,options,1);
    %qual3Db(isnan(qual_(:))) = elem_qual([edg_(R(Ipdt),2) I123(Ipdt,[2 3 1])], xy,Nmetric,options,1);
    qual3Da(isnan(qual_(:))) = elem_qual([I123(Ipdt,[1 2 3]) edg_(R(Ipdt),1)], xy,Nmetric,options,1);
    qual3Db(isnan(qual_(:))) = elem_qual([I123(Ipdt,[3 2 1]) edg_(R(Ipdt),2)], xy,Nmetric,options,1);
    qual3Dmin = min(qual3Da,qual3Db);
    qual_(isnan(qual_(:))) = qual3Dmin(isnan(qual_(:)));
    else %3D
    qual_(isnan(qual_(:))) = elem_qual(I123(Ipdt,:),xy,Nmetric,options,edg_(R(Ipdt),:));
    end;
    
    if options.qualRM || size(edg_,2) == 2
	qual_fix = qual_;
    else 
    	qual2(isnan(qual2(:))) = 1./elem_qual(I123(Ipdt,[1 3]),xy,Nmetric,options);
    	qual_fix = qual2;
    	qual_fix(qual_ < 0) = 0;
    end;
    if not(allconvex) %only triangulate next to concave 
    Ncnvx = any(qual_<0);
    if any(Ncnvx)
    qual_fix_ = qual_fix;
    qual_fix(:,Ncnvx) = NaN;
    Cc = Cc(Ncnvx(R)); Cr = Cr(Ncnvx(R)); Cl = Cl(Ncnvx(R)); R = R(Ncnvx(R));
    Ifix = qual_(Cc+(R-1)*size(qual_,1))<0; 
    qual_fix(Cl(Ifix)+(R(Ifix)-1)*size(qual_,1)) = qual_fix_(Cl(Ifix)+(R(Ifix)-1)*size(qual_,1));
    qual_fix(Cr(Ifix)+(R(Ifix)-1)*size(qual_,1)) = qual_fix_(Cr(Ifix)+(R(Ifix)-1)*size(qual_,1));
    else
    	allconvex = true;
    end;
    end;
    
    %best candidate element found
    [tmp,CI] = max(qual_fix); CI = CI';
    RI = 1:nbad;
    %give up on removing node, if best candidate element is 
    if (options.consRM || size(edg_,2) == 2) && size(edg_,2) ~= 3 %worse than previously
    Ibetter = qual_(CI+(RI'-1)*size(qual_,1)) > triQtb(badnds); 
    else %inverted (we do not store surface quality in 3D)
    Ibetter = qual_(CI+(RI'-1)*size(qual_,1)) > options.minqual; %0 is not safe
    end;
    
    if not(all(Ibetter))
      circles = circles(Ibetter,:); nbad = size(circles,1); 
      qual_   = qual_(:,Ibetter);
      RI = 1:nbad; CI = CI(Ibetter); NN = NN(Ibetter);
      badbadnds = [badbadnds; badnds(not(Ibetter))];
      badnds = badnds(Ibetter);
      badIDs = badIDs(Ibetter);
      if nbad == 0
        break;
      end;
    end;
    Cl = thmp1(CI); I = find(Cl==0); Cl(I) = NN(RI(I));
    Cr = thmp2(CI); I = find(Cr==thmp2(NN(RI))); Cr(I) = 1;
    if size(edg_,2) ~= 0
        edg_ = edg_(Ibetter,:);
    end;
    if size(xy,2) == 2 || size(edg_,2) == 3
    newtri_ = [circles(RI'+(Cr-1)*nbad) circles(RI'+(CI-1)*nbad) circles(RI'+(Cl-1)*nbad)];
    newtriN = [newtriN; badnds];
    newIDs = [newIDs; badIDs];
    qualityN = [qualityN; qual_(CI+(RI'-1)*size(qual_,1))];
    newedg = [newedg; newtri_(NN~=3,[1 3])];
    else %swapping
    qual3Da = qual3Da(:,Ibetter);
    qual3Db = qual3Db(:,Ibetter);
    %newtri_ = [[edg_(:,1) circles(RI'+(Cl-1)*nbad) circles(RI'+(CI-1)*nbad) circles(RI'+(Cr-1)*nbad)]; [edg_(:,2) circles(RI'+(CI-1)*nbad) circles(RI'+(Cl-1)*nbad) circles(RI'+(Cr-1)*nbad)]];
    newtri_ = [[circles(RI'+(Cr-1)*nbad) circles(RI'+(CI-1)*nbad)  circles(RI'+(Cl-1)*nbad) edg_(:,1)]; [circles(RI'+(Cl-1)*nbad) circles(RI'+(CI-1)*nbad)  circles(RI'+(Cr-1)*nbad) edg_(:,2)]];
    newtriN = [newtriN; repmat(badnds,2,1)]; 
    newIDs  = [newIDs;  repmat(badIDs,2,1)];
    qualityN = [qualityN; qual3Da(CI+(RI'-1)*size(qual_,1)); qual3Db(CI+(RI'-1)*size(qual_,1))];
    uppe = [false(nbad,1); true(nbad,1)]; lowe = not(uppe);
    NN1n3 = [NN~=3; false(size(NN))];
    uppe = [false(nbad,1); true(nbad,1)];
    newedg = [newedg; newtri_(NN1n3,[1 3])];
    newfac = [newfac; newtri_(repmat(NN~=3,2,1),[1 3 4]); newtri_(uppe,1:3)];
    newfacN = [newfacN; repmat(badnds(NN~=3),2,1); badnds];
    end;
    newedgN = [newedgN; badnds(NN~=3)];
    newtri = [newtri; newtri_];
    % UPDATE CIRCLES
    circles(RI'+(CI-1)*nbad) = 0;
    qual_(CI+(RI'-1)*size(qual_,1)) = 0;
    qual_(Cr+(RI'-1)*size(qual_,1)) = NaN;
    qual_(Cl+(RI'-1)*size(qual_,1)) = NaN;
    NN = NN - 1;
    I = NN>2; NN = NN(I);
    if not(any(NN>2))
    	break;
    end;
    circles = circles(I,:); qual_ = qual_(:,I); badnds = badnds(I); badIDs = badIDs(I);
    
    if size(edg_,2) == 2
    	qual3Da = qual3Da(:,I);
    	qual3Db = qual3Db(:,I);
    end;
    if size(edg_,2) ~= 0
    	edg_ = edg_(I,:);
    end;
    if size(edg_,2) == 2
    	[circles,qual_,qual3Da,qual3Db] = fix_circles(circles,qual_,qual3Da,qual3Db);
    elseif options.qualRM 
	[circles,qual_] = fix_circles(circles,qual_);
    else
    	qual2 = qual2(:,Ibetter); %if not(all(Ibetter))
    	qual2(CI+(RI'-1)*size(qual2,1)) = 0;
      	qual2(Cr+(RI'-1)*size(qual2,1)) = NaN;
      	qual2(Cl+(RI'-1)*size(qual2,1)) = NaN;
    	qual2 = qual2(:,I);
    	[circles,qual_,qual2] = fix_circles(circles,qual_,qual2);
    end;
end; %while
%end; %ifany
Ikeep = true(max(badnds_start),1); Ikeep(badbadnds) = false;
badnds = badnds_start(Ikeep(badnds_start));
newtri   = newtri(  Ikeep(newtriN),:);
qualityN = qualityN(Ikeep(newtriN));
newIDs = newIDs(Ikeep(newtriN));
newtriN  = newtriN( Ikeep(newtriN));
newedg   = newedg( Ikeep(newedgN),:);
newedgN  = newedgN( Ikeep(newedgN));
if size(edg_,2) == 2
newfac   = newfac( Ikeep(newfacN),:);
newfacN  = newfacN(Ikeep(newfacN));
end;

function [circles,oqual_,oqual2,oqual3] = fix_circles(circles,qual_,qual2,qual3)
ocircles = circles';
[Cg,Rg] = find(ocircles~=0);
[C ,R ] = find(ocircles==0);
move = zeros(size(circles)); move(R+(C-1)*size(circles,1)) = -1; move = cumsum(move')';
Cg = Cg + reshape(move(Rg+(Cg-1)*size(circles,1)),size(Cg));
circles = zeros(size(circles));
circles(Rg+(Cg-1)*size(circles,1)) = ocircles(ocircles~=0);
if nargout >= 2
oqual_ = zeros(size(circles));
oqual_(Rg+(Cg-1)*size(circles,1)) = qual_(find(ocircles)); oqual_ = oqual_';
end;
if nargout >= 3
oqual2 = zeros(size(circles));
oqual2(Rg+(Cg-1)*size(circles,1)) = qual2(find(ocircles)); oqual2 = oqual2';
end;
if nargout >= 4
oqual3 = zeros(size(circles));
oqual3(Rg+(Cg-1)*size(circles,1)) = qual3(find(ocircles)); oqual3 = oqual3';
end;



function [newtri,circles,badnds,badbadnds,newtriN,qualityN,R2R,NN,edg_] = rm_trinds(circles,badnds,badbadnds,xy,Nmetric,triQtb,options,edg_)
NN = sum(circles>0,2);
N3 = NN==3;
R2R = 1:numel(badnds);	
if not(any(N3))
  newtri = []; qualityN = []; newtriN = [];
  return;
end;
if size(xy,2) == 2 || size(edg_,2) == 3
  newtri = circles(N3,[3 2 1]); 
  newtriN = badnds(N3)';
  qualityN = elem_qual(newtri,xy,Nmetric,options); %2D element uninvertible
else %swapping
  newtri = [[edg_(N3,1) circles(N3,1:3)]; ...
            [edg_(N3,2) circles(N3,[2 1 3])]]; 
  qualityN = elem_qual(newtri,xy,Nmetric,options,1);
  newtriN = repmat(badnds(N3)',2,1);
end;
if (options.consRM || size(edg_,2) == 2) && size(edg_,2) ~= 3
  if size(edg_,2) ~= 2
   Ibetter = triQtb(badnds(N3)) < qualityN;
   Ibetter2 = Ibetter;
  elseif size(edg_,2) ~= 3 %swapping
   Ibetter = triQtb(badnds(N3)) < min(reshape(qualityN,[sum(N3) 2]),[],2);
   Ibetter2 = repmat(Ibetter,2,1);
  else %we do not store surface quality in 3D
   Ibetter = 0 < qualityN;
   Ibetter2 = Ibetter;
  end;
  newtri = newtri(Ibetter2,:);
  newtriN = newtriN(Ibetter2);
  qualityN = qualityN(Ibetter2);
  Igood_ = not(N3); Igood_(N3) = Ibetter;
  badbadnds = [badbadnds; badnds(not(Igood_))'];
  badnds = badnds(Igood_); circles = circles(Igood_,:); N3 = N3(Igood_); NN = NN(Igood_);
  if size(xy,2) == 3
  	edg_ = edg_(Igood_,:); edg_ = edg_(not(N3),:); 
  end;
  R2R = (1:numel(badnds))';
end;
circles = circles(not(N3),:); R2R = R2R(not(N3)); NN = NN(not(N3));


function [circles,badnds,edg,badbadnds] = del_dblNcnvx(circles,xy,badnds,edg,options)
NN = sum(circles~=0,2);
[Cc,R] = find(circles'); nbad = size(circles,1);
thmp1 = [0:size(circles,2)-1]'; thmp2 = [2:size(circles,2)+1]';
Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
Cr = thmp2(Cc); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
I1 = circles(R+(Cr-1)*nbad); I2 = circles(R+(Cc-1)*nbad); I3 = circles(R+(Cl-1)*nbad);
I123 = [I1(:) I2(:) I3(:)];
if size(xy,2) == 2
[Ibad,area] = elem_inv(I123,xy);
elseif size(edg,2) == 2 %swapping
[Ibad,area] = elem_inv([[edg(R,1) I123(:,[3 2 1])]; [edg(R,2) I123(:,[2 3 1])] ],xy); 
Ibad = any(reshape(Ibad,size(circles,1),2),2);
else
[Ibad,area] = elem_inv(I123,xy,edg(R,:)); 
end;%Ibad = area < options.minA;
bcircles = false(size(circles));
bcircles(R+(Cc-1)*nbad) = Ibad;
if options.rmnonconvec == 0
	ndblNcnvx = sum(bcircles,2) == 0;
else
	ndblNcnvx = sum(bcircles,2) < 2;
end;
circles = circles(ndblNcnvx,:);
badbadnds = badnds(not(ndblNcnvx));
badnds = badnds(ndblNcnvx);
if size(xy,2) == 3
	edg    = edg(ndblNcnvx,:);
end;

function [newtri,qualityN,badnds,newtriN,badbadnds,newIDs,newedg,newedgN] = prag_coarse(circles,Nmetric,xy,badnds,badIDs,triQtb,NN,options,nvec)
%save for_debug3.mat; error('just stop');
badbadnds = zeros(0,1); badnds_start = badnds;
newtri = zeros(0,3); newtriN = []; qualityN = []; newIDs = []; newedg = zeros(0,2); newedgN = [];

[Cc,R] = find(circles'~=0); nbad = size(circles,1);
edg = [badnds(R) reshape(circles(R+(Cc-1)*nbad),size(R))];
edgL = inf(size(circles));
edgL(R+(Cc-1)*nbad) = elem_qual(edg,xy,Nmetric,options);
[edgL_,Iedg] = sort(edgL,2);
Iedg(isinf(edgL)) = 0;
theC = 1;

while true
[Cc,R] = find(circles'~=0); nbad = size(circles,1);
thmp1 = [0:size(circles,2)-1]'; Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
I2 = circles(R+(Cc-1)*nbad); I3 = circles(R+(Cl-1)*nbad);
C_ = Iedg(R+(theC-1)*nbad);
Ishrt = circles(R+(C_(:)-1)*nbad);
newtri_ = [Ishrt(:) I2(:) I3(:)];

%[newtri__,tmp] = sort(newtri_,2);
%Igood = and(newtri__(:,1) ~= newtri__(:,2),newtri__(:,2) ~= newtri__(:,3));

thmp2 = [2:size(circles,2)+1]';
Cl_ = thmp1(C_);  I = find(Cl_==0); Cl_(I) = NN(R(I));
Cr_ = thmp2(C_);  I = find(Cr_==thmp2(NN(R))); Cr_(I) = 1;
Igood = and(newtri_(:,2) ~= reshape(circles(R+(Cr_-1)*nbad),numel(R),1), ...
            newtri_(:,3) ~= reshape(circles(R+(Cl_-1)*nbad),numel(R),1));
Cr_ = Cr_(Igood); Cl_ = Cl_(Igood); C_ = C_(Igood);

newtri_ = newtri_(Igood,:); R = R(Igood); Cc = Cc(Igood);
if size(xy,2) == 2
   newqual = elem_qual(newtri_,xy,Nmetric,options,1);
else
   newqual = elem_qual(newtri_,xy,Nmetric,options,nvec(R,:));
end;
if options.consRM && size(nvec,2) ~= 3
	Igood = newqual > triQtb(badnds(R));
else
	Igood = newqual > options.minqual;
end;
allgood = true(size(circles)); allgood(R+(Cc-1)*nbad) = Igood;
Igood = all(allgood,2); Ibad = not(Igood); IgoodR = Igood(R);
newtri = [newtri; newtri_(IgoodR,:)];
qualityN = [qualityN; newqual(IgoodR)];
newtriN  = [newtriN ; badnds(R(IgoodR))];
newIDs   = [newIDs  ; badIDs(R(IgoodR))];
if nargout > 6
%thmp2 = [2:size(circles,2)+1]'; C_ = Iedg(R+(theC-1)*nbad); 
%Cl = thmp1(C_); I = find(Cl==0); Cl(I) = NN(R(I));
%Cr = thmp2(C_); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
newedg_ = newtri_(IgoodR,1:2); newedgN_ = badnds(R(IgoodR));
Ie = and(newedg_(:,2) ~= reshape(circles(R(IgoodR)+(Cl_(IgoodR)-1)*nbad),sum(IgoodR),1), ...
         newedg_(:,2) ~= reshape(circles(R(IgoodR)+(Cr_(IgoodR)-1)*nbad),sum(IgoodR),1));
newedg = [newedg; newedg_(Ie,:)]; newedgN = [newedgN; newedgN_(Ie)];
end;
%update circles
circles = circles(Ibad,:); badnds = badnds(Ibad); NN = NN(Ibad); Iedg = Iedg(Ibad,:); badIDs = badIDs(Ibad);
if size(xy,2) == 3
	nvec = nvec(Ibad,:);
end; 
if size(circles,1) == 0
	break;
end;
theC = theC+1;
Ibad = NN < theC; Igood = not(Ibad);
badbadnds = [badbadnds; badnds(Ibad)];
circles = circles(Igood,:); badnds = badnds(Igood); NN = NN(Igood); Iedg = Iedg(Igood,:); badIDs = badIDs(Igood);
if size(xy,2) == 3
	nvec = nvec(Igood,:);
end;
if size(circles,1) == 0
	break;
end; %size(circles)
end; %badbadnds
Ikeep = true(max(badnds_start),1); Ikeep(badbadnds) = false;
badnds = badnds_start(Ikeep(badnds_start));
%next four lines are irrelevant outside the context of region IDs
newtri   = newtri(  Ikeep(newtriN),:);
qualityN = qualityN(Ikeep(newtriN));
newIDs = newIDs(Ikeep(newtriN));
newtriN  = newtriN(Ikeep(newtriN));
if nargout > 6
newedg = newedg(Ikeep(newedgN),:);
newedgN = newedgN(Ikeep(newedgN));
end;


%hold on; trimesh(newtri,xy(:,1),xy(:,2)); hold off;
function [newtri,qualityN,badnds,newtriN,badbadnds,newIDs,newfac,newedg] = fill_spheres(spheres,Nmetric,xy,badnds,badIDs,triQtb,options)
spheres = fix_spheres(spheres);
spheres = spheres(:,1:max(sum(squeeze(spheres(1,:,:))~=0)),:);
if options.prag_crs
	[newtri,qualityN,badnds,newtriN,badbadnds,newIDs] = prag_crs(spheres,Nmetric,xy,badnds,badIDs,triQtb,options); return; end;
NN = sum(squeeze(spheres(1,:,:))~=0)'; 
badnds_start = badnds;
newtri = []; newtriN = []; qualityN = []; badbadnds = []; newIDs = [];
newedg = zeros(0,2); newedgN = []; newfac = zeros(0,3); newfacN = [];
if options.area == 0. %curved geometry
	[spheres,NN,badnds,badIDs,badbadnds] = kill_split_spheres(spheres,NN,badnds,badIDs,badbadnds);
	%edg2tri_ = test_spheres(spheres,xy,badnds);
end;
while true
while true
[newtri,newfac_,spheres,badbadnds,newtriN,newIDs,newfacN_,qualityN,badnds,badIDs,NN,ndone] = rm_frnds(spheres,xy,NN,badnds,badIDs,badbadnds,Nmetric,triQtb,newtriN,newIDs,qualityN,newtri,options); %save for_debug3D.mat;
newfac = [newfac; newfac_];
newfacN = [newfacN; newfacN_];
%if ndone && size(spheres,3) ~= 0
%	spheres = fix_spheres(spheres);
%else 
if not(ndone) || size(spheres,3) == 0
	break;
end;
end; %while
if size(spheres,3) == 0 %nnz(spheres(1,:,:))
	break;
end;

%save for_debug3Db.mat;
%make an edge (and an element)
I = squeeze(spheres(1,:,:))~=0;
R_ = repmat(1:size(spheres,3),size(spheres,2),1);
tris = reshape(repmat(1:size(spheres,2),1,size(spheres,3)),[size(spheres,2) size(spheres,3)]); tris = repmat(tris(I),3,1);
edg = [[spheres(1,I)'; spheres(2,I)'; spheres(3,I)'] [spheres(2,I)'; spheres(3,I)'; spheres(1,I)']];
[edg,Ieven] = sort(edg,2); Ieven = Ieven(:,1) == 1;
edg = [repmat(R_(I),3,1) edg];
[edg,I] = sortrows(edg); tris = tris(I); Ieven = Ieven(I);
d = [true; or(edg(1:end-1,1) ~= edg(2:end,1),or(edg(1:end-1,2) ~= edg(2:end,2),edg(1:end-1,3) ~= edg(2:end,3)))]; ad = not(d);
edg2tri = [tris(d) tris(ad)];
edg = edg(d,:);

Iflp = Ieven(d); edg2tri(Iflp,:) = edg2tri(Iflp,[2 1]);
edgalt = [spheres(:,edg2tri(:,1)+(edg(:,1)-1)*size(spheres,2)); spheres(:,edg2tri(:,2)+(edg(:,1)-1)*size(spheres,2))];
edgalt(edgalt==repmat(edg(:,2)',6,1)) = 0;
edgalt(edgalt==repmat(edg(:,3)',6,1)) = 0;
edgalt = reshape(edgalt(edgalt~=0),2,size(edg,1))';
newtri_ = [edgalt edg(:,2:3)];
qual_ = elem_qual(newtri_,xy,Nmetric,options,[]);
if any(options.RMnd3D == [1 2])
	qual_(gen_badedgalt(edgalt,edg)) = 0; %we prevent splitting of spheres
end;


nN_ = [0; find([edg(1:end-1,1)~=edg(2:end,1); true])];
nN = diff(nN_); nMax = max(nN);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN',nMax,1); I = C <= CnN;
if options.qualRM
qual = zeros(numel(nN),nMax); qual(edg(:,1)+(C(I)-1)*numel(nN)) = qual_;
[maxqual,CI] = max(qual,[],2);
else
qual2_ = 1./elem_qual(edgalt,xy,Nmetric,options);
qual2_(qual_<=0) = 0; %no inversion
qual2 = zeros(numel(nN),nMax); qual2(edg(:,1)+(C(I)-1)*numel(nN)) = qual2_;
[maxqual,CI] = max(qual2,[],2);
end;
qualityN_ = zeros(size(CI)); Igood = maxqual>0;
if any(Igood)
qualityN_(Igood) = qual_(CI(Igood)+nN_([Igood; false]));
end;
if options.consRM
	Ibetter = qualityN_*options.minchg > triQtb(badnds); %badnds(R2R)
else
	Ibetter = qualityN_ > options.minqual;
end;
if not(all(Ibetter))
badbadnds = [badbadnds; badnds(not(Ibetter))]; 
end;
if any(Ibetter)
RC = CI(Ibetter)+nN_([Ibetter; false]);
else
RC = [];
end;
newtri = [newtri; newtri_(RC,:)];
qualityN = [qualityN; qualityN_(Ibetter)];
newtriN = [newtriN; badnds(Ibetter)]; 
newIDs = [newIDs; badIDs(Ibetter)];
NN = NN(Ibetter);
badnds = badnds(Ibetter); 
badIDs = badIDs(Ibetter);
newedg = [newedg; newtri_(RC,[1 2])];
newedgN = [newedgN; badnds];
newfac = [newfac; newtri_(RC,[1 2 3]); newtri_(RC,[1 2 4])];
newfacN = [newfacN; repmat(badnds,2,1)];
spheres(:,edg2tri(RC,1)+(edg(RC,1)-1)*size(spheres,2)) = [edgalt(RC,:) edg(RC,2)]';
spheres(:,edg2tri(RC,2)+(edg(RC,1)-1)*size(spheres,2)) = [edgalt(RC,[2 1]) edg(RC,3)]';
spheres = spheres(:,:,Ibetter);
if size(spheres,3) == 0 
	break;
end;
if not(any(options.RMnd3D == [1 2]))
[spheres,NN,badnds,badIDs,badbadnds] = kill_split_spheres(spheres,NN,badnds,badIDs,badbadnds);
end;
if options.debug == 2
edg2tri = test_spheres(spheres,xy,badnds);
end;
end; %while
Ikeep = true(max(badnds_start),1); 

Ikeep(badbadnds) = false;
badnds = badnds_start(Ikeep(badnds_start));
newtri   = newtri(  Ikeep(newtriN),:);
qualityN = qualityN(Ikeep(newtriN));
newIDs  = newIDs( Ikeep(newtriN));
newtriN  = newtriN( Ikeep(newtriN));
newfac = newfac(Ikeep(newfacN),:);
newfacN = newfacN(Ikeep(newfacN));
newedg = newedg(Ikeep(newedgN),:);
newedgN = newedgN(Ikeep(newedgN));

function [newtri,newfac,spheres,badbadnds,newtriN,newIDs,newfacN,qualityN,badnds,badIDs,NN,ndone] = rm_frnds(spheres,xy,NN,badnds,badIDs,badbadnds,Nmetric,triQtb,newtriN,newIDs,qualityN,newtri,options)
ndone = 0; newfac = zeros(0,3); newfacN = [];
R_ = reshape(repmat(1:size(spheres,3),3*size(spheres,2),1),size(spheres));
tris = reshape(repmat(1:size(spheres,2),3,size(spheres,3)),size(spheres));
Io = spheres~=0; tris = tris(Io);
nd2tri = [R_(Io) spheres(Io)];

[nd2tri,I] = sortrows(nd2tri); tris = tris(I);
I3 = and(and(nd2tri(2:end-1,1)==nd2tri(1:end-2,1), nd2tri(2:end-1,2)==nd2tri(1:end-2,2)), and(nd2tri(2:end-1,1)==nd2tri(3:end,1), nd2tri(2:end-1,2)==nd2tri(3:end,2)));
I3_ = [or(nd2tri(2,1)~=nd2tri(4,1),nd2tri(2,2)~=nd2tri(4,2)); ...
and(or(nd2tri(3:end-2,1)~=nd2tri(1:end-4,1), nd2tri(3:end-2,2)~=nd2tri(1:end-4,2)), or(nd2tri(3:end-2,1)~=nd2tri(5:end,1),nd2tri(3:end-2,2)~=nd2tri(5:end,2))); ...
or(nd2tri(end-1,1)~=nd2tri(end-3,1),nd2tri(end-1,2)~=nd2tri(end-3,2))];
d1 = and(I3,I3_);
if not(any(d1))
	return
end; 

%fix in case of tetrahedrons in spheres
I_ = NN==4;
if any(I_)
	[tmp,fC] = max(spheres(1,:,:)~=0,[],2); %no fix spheres
	first = reshape( repmat(spheres(1,squeeze(fC)+(0:size(spheres,3)-1)'*size(spheres,2)),3*size(spheres,2),1),size(spheres));
	%first = reshape( repmat(squeeze(spheres(1,1,:))',3*size(spheres,2),1),size(spheres));
	first = first(Io); first = first(I);
	I_ = and(I_(nd2tri(:,1)),nd2tri(:,2)~=first);
	d1(I_(2:end-1)) = false;
end;

d2 = [false; false; d1];
d3 = [d1; false; false];
d1 = [false; d1; false];
t1 = tris(d1)+(nd2tri(d1,1)-1)*size(spheres,2)*(size(spheres,3)~=1);
t2 = tris(d2)+(nd2tri(d2,1)-1)*size(spheres,2)*(size(spheres,3)~=1);
t3 = tris(d3)+(nd2tri(d3,1)-1)*size(spheres,2)*(size(spheres,3)~=1);
newtri1 = spheres(:,t1);
newtri2 = spheres(:,t2);
newtri3 = spheres(:,t3);
I1 = newtri1 == repmat(nd2tri(d1,2)',3,1); [C1,R] = find(I1);
I2 = newtri2 == repmat(nd2tri(d2,2)',3,1); [C2,R] = find(I2);
I3 = newtri3 == repmat(nd2tri(d3,2)',3,1); [C3,R] = find(I3);
thmp = [3 1 2]'; C1a = thmp(C1);
thmp = [2 3 1]'; C1 = thmp(C1); C2 = thmp(C2); C3 = thmp(C3);
R = (R-1)*3; newtri_ = [newtri1(C1+R) newtri2(C2+R) newtri3(C3+R) nd2tri(d1,2)];
Ieven = newtri1(C1a+R) ~= newtri2(C2+R);
newtri_(Ieven,1:2) = newtri_(Ieven,[2 1]); 
qualityN_ = elem_qual(newtri_,xy,Nmetric,options,[]);

Rs = nd2tri(d1,1); 
if numel(Rs) ~= 1
d = [Rs(1:end-1)~=Rs(2:end); true];
else
d = true;
end; %save for_debug3D.mat;
nN_ = [0; find(d)]; nN = diff(nN_); nMax = max(nN);
Rs2R = zeros(max(Rs),1); Rs2R(Rs(d)) = 1:numel(nN);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN',nMax,1); I = C <= CnN; 
CI = reshape(C(I),numel(Rs),1); %reshape relevant for nMax==1
qual = inf(numel(nN),nMax); qual(Rs2R(Rs)+(CI-1)*numel(nN)) = qualityN_;
Rsd = Rs(d);

if options.consRM || nargin == 7
Ibetter = min(qual,[],2)*options.minchg > triQtb(badnds(Rsd));
else
Ibetter = min(qual,[],2) > options.minqual;
end;

if any(options.RMnd3D == [2 3]) %we will not kill spheres
IR = Ibetter(Rs2R(Rs)); nN = nN(Ibetter);
Rs = Rs(IR); d = d(IR);
t1 = t1(IR); t2 = t2(IR); t3 = t3(IR);
newtri_ = newtri_(IR,:);
qualityN_ = qualityN_(IR);
end;

%the following three lines are not dangerous, when you account for tetrahedrons in spheres
if any(d)
ndone = 1;
spheres(:,t1) = newtri_(:,1:3)';
spheres(:,t2) = 0;
spheres(:,t3) = 0;
newtriN_ = badnds(Rs);
newIDs_ = badIDs(Rs);
newtriN = [newtriN; newtriN_];
newIDs = [newIDs; newIDs_];
newtri = [newtri; newtri_];
qualityN = [qualityN; qualityN_];
I1 = or(NN(Rs)>6,NN(Rs)==4); I2 = NN(Rs)==6;
newfac = newtri_(I1,1:3); %incomplete
newfacN = newtriN_(I1); %incomplete
NN(Rs(d)) = NN(Rs(d))-2*nN;
end;
if any(not(Ibetter)) && not(any(options.RMnd3D == [2 3]))
Ibttr = true(size(spheres,3),1); Ibttr(Rsd(not(Ibetter))) = false;
badbadnds = [badbadnds; badnds(not(Ibttr))];
spheres = spheres(:,:,Ibttr);
badnds = badnds(Ibttr); 
badIDs = badIDs(Ibttr);
NN = NN(Ibttr);
end;
%finishes spheres
I = NN>3; spheres = spheres(:,:,I); NN = NN(I); badnds = badnds(I); badIDs = badIDs(I);
if size(spheres,3)~=0 && options.debug == 2
	edg2tri = test_spheres(spheres);
end;

function badedgalt = gen_badedgalt(edgalt,edg) %no semi intersecting spheres
nmax = max(max(edg(:,2:3)));
rmax = max(edg(:,1));
edg(:,1) = edg(:,1) + nmax;
nd2edg = inv_table(edg(:,1));
edgalt = sort(edgalt,2); 
inds = nd2edg(edg(:,1),:); I = inds ~= 0;
edg1 = zeros(size(inds)); edg1(I) = edg(inds(I),2);
edg2 = zeros(size(inds)); edg2(I) = edg(inds(I),3);
badedgalt = any(and(edg1==repmat(edgalt(:,1),1,size(nd2edg,2)),edg2==repmat(edgalt(:,2),1,size(nd2edg,2))),2);

%A1 = sparse(edg(:,2),edg(:,3),edg(:,1),nmax,nmax);
%A3 = sparse(edgalt(:,1),edgalt(:,2),1:size(edg,1),nmax,nmax);
%%edg and edgalt might be dubplicated across different sets!
%A2 = sparse(edgalt(:,1),edgalt(:,2),edg(:,1),nmax,nmax);
%I = and(A1==A2,A1~=0);
%badedgalt = full(A3(I));
%

function ospheres = fix_spheres(spheres)
spheres = permute(spheres,[2 1 3]);
[Cg,Rg] = find(spheres~=0);
[C ,R ] = find(spheres==0);
move = zeros(size(spheres)); move(C+(R-1)*size(spheres,1)) = -1; move = cumsum(move);
Cg = Cg + reshape(move(Cg+(Rg-1)*size(spheres,1)),size(Cg));
ospheres = zeros(size(spheres));
ospheres(Cg+(Rg-1)*size(ospheres,1)) = spheres(spheres~=0);
ospheres = permute(ospheres,[2 1 3]);

function [spheres,NN,badnds,badIDs,badbadnds] = kill_split_spheres(spheres,NN,badnds,badIDs,badbadnds)
I = squeeze(spheres(1,:,:))~=0;
R_ = repmat(1:size(spheres,3),size(spheres,2),1);
tris = reshape(repmat(1:size(spheres,2),1,size(spheres,3)),[size(spheres,2) size(spheres,3)]); tris = repmat(tris(I),3,1);
edg = [[spheres(1,I)'; spheres(2,I)'; spheres(3,I)'] [spheres(2,I)'; spheres(3,I)'; spheres(1,I)']];
edg = sort(edg,2);
edg = [repmat(R_(I),3,1) edg];
[edg,I] = sortrows(edg); 
d = [false; false; and(and(edg(3:end-1,1)==edg(4:end,1),edg(3:end-1,1)==edg(1:end-3,1)), ...
and(and(edg(3:end-1,2)==edg(4:end,2),edg(3:end-1,2)==edg(1:end-3,2)), ...
    and(edg(3:end-1,3)==edg(4:end,3),edg(3:end-1,3)==edg(1:end-3,3)))); false];
if any(d)
	Rbad = edg(d,1);
	badbadnds = [badbadnds; badnds(Rbad)]; 
	Ikeep = true(size(spheres,3),1); Ikeep(Rbad) = false;
	spheres = spheres(:,:,Ikeep); badnds = badnds(Ikeep); badIDs = badIDs(Ikeep);
	if numel(NN) ~= 0
	NN = NN(Ikeep);  
	end;
end;

function edg2tri = test_spheres(spheres,xy,badnds)
I = squeeze(spheres(1,:,:))~=0;
R_ = repmat(1:size(spheres,3),size(spheres,2),1);
fac = sortrows([R_(I) sort([spheres(1,I)' spheres(2,I)' spheres(3,I)'],2)]);
d = all(fac(1:end-1,:) == fac(2:end,:),2);
if any(d)
	save for_debug3D.mat; 
	error('double face in spheres');
end;
tris = reshape(repmat(1:size(spheres,2),1,size(spheres,3)),[size(spheres,2) size(spheres,3)]); tris = repmat(tris(I),3,1);
edg = [[spheres(1,I)'; spheres(2,I)'; spheres(3,I)'] [spheres(2,I)'; spheres(3,I)'; spheres(1,I)']];
[edg,Ieven] = sort(edg,2); Ieven = Ieven(:,1) == 1;
edg = [repmat(R_(I),3,1) edg];
[edg,I] = sortrows(edg); tris = tris(I); Ieven = Ieven(I);
d = [true; or(edg(1:end-1,1)~=edg(2:end,1),or(edg(1:end-1,2) ~= edg(2:end,2),edg(1:end-1,3) ~= edg(2:end,3)))]; ad = not(d);
edg = edg(d,:);
edg2tri = [tris(d) tris(ad)];

function [newtri,qualityN,badnds,newtriN,badbadnds,newIDs] = prag_crs(spheres,Nmetric,xy,badnds,badIDs,triQtb,options)
badbadnds = zeros(0,1); 
if options.area == 0. %curved geometry
	[spheres,NN,badnds,badIDs,badbadnds] = kill_split_spheres(spheres,[],badnds,badIDs,badbadnds);
	%edg2tri_ = test_spheres(spheres,xy,badnds);
end;
newtri = zeros(0,4); newtriN = []; qualityN = []; badnds_start = badnds; newIDs = [];

nbad = size(spheres,3); nbad_ = size(spheres,2); nbadE = 3*nbad_;
[Cc,R] = find(reshape(spheres,[nbadE,size(spheres,3)])~=0); 
edg = [badnds(R) reshape(spheres(Cc+(R-1)*nbadE),size(R))];
[edg_,I] = sortrows(edg); [tmp,Iback] = sort(I);
Igood = [true; or(edg_(1:end-1,1) ~= edg_(2:end,1),edg_(1:end-1,2) ~= edg_(2:end,2))];
edg = edg(Igood(Iback),:); R = R(Igood(Iback)); Cc = Cc(Igood(Iback));
edgL = inf(  nbad,nbadE);
edgS = zeros(nbad,nbadE);
edgL(R+(Cc-1)*nbad) = elem_qual(edg,xy,Nmetric,options);
edgS(R+(Cc-1)*nbad) = edg(:,2);
[edgL_,Iedg] = sort(edgL,2);
edgS = reshape(edgS(reshape(repmat((1:nbad)',1,nbadE),nbad*nbadE,1)+(Iedg(:)-1)*nbad),size(edgS));
%Iedg(isinf(edgL)) = 0; 
theC = 1;
NN = sum(edgS~=0,2);
while true
[Cc,R] = find(squeeze(spheres(1,:,:))~=0); nbad = size(spheres,3);
if nbad==1
	R_ = R; R = reshape(Cc,numel(Cc),1); Cc = reshape(R_,size(R));
end;
%C_ = Iedg(R+(theC-1)*nbad);
thshrt = reshape(edgS(R+(theC-1)*nbad),numel(R),1);
newtri_ = [thshrt spheres(1,Cc+(R-1)*nbad_)' spheres(2,Cc+(R-1)*nbad_)' spheres(3,Cc+(R-1)*nbad_)'];
[newtri__,tmp] = sort(newtri_,2);
Igood = and(and(newtri__(:,1)~=newtri__(:,2),newtri__(:,2)~=newtri__(:,3)),newtri__(:,3)~=newtri__(:,4));
newtri_ = newtri_(Igood,:); Cc = Cc(Igood); R = R(Igood);
newqual = elem_qual(newtri_,xy,Nmetric,options,1);
if options.consRM
	Igood = newqual > triQtb(badnds(R));
else
	Igood = newqual > options.minqual;
end;
allgood = true(nbad,nbad_); allgood(R+(Cc-1)*nbad) = Igood;
Igood = all(allgood,2); Ibad = not(Igood);
newtri = [newtri; newtri_(Igood(R),:)];
qualityN = [qualityN; newqual(Igood(R))];
newtriN  = [newtriN ; badnds(R(Igood(R)))];
newIDs   = [newIDs  ; badIDs(R(Igood(R)))];
%update circles
spheres = spheres(:,:,Ibad); badnds = badnds(Ibad); NN = NN(Ibad); edgS = edgS(Ibad,:); badIDs = badIDs(Ibad);
if size(spheres,3) == 0
	break;
end;
theC = theC+1;
Ibad = NN < theC; Igood = not(Ibad);
badbadnds = [badbadnds; badnds(Ibad)];
spheres = spheres(:,:,Igood); badnds = badnds(Igood); NN = NN(Igood); edgS = edgS(Igood,:); badIDs = badIDs(Igood);
if size(spheres,3) == 0
	break;
end; %[size(spheres,3) size(edgS,1) theC]
end;
Ikeep = true(max(badnds_start),1); Ikeep(badbadnds) = false;
badnds = badnds_start(Ikeep(badnds_start));
%next four lines are irrelevant outside the context of region IDs
newtri   = newtri(  Ikeep(newtriN),:);
qualityN = qualityN(Ikeep(newtriN));
newIDs = newIDs(Ikeep(newtriN));
newtriN  = newtriN( Ikeep(newtriN));
function edgalt = bks_edgalt(edg,tri,edg2tri,bndedg_)
edgs = size(edg,1);
dim = size(tri,2);
bndedg = false(edgs,1); bndedg(bndedg_) = true;
Iedg = not(bndedg); edgs = edgs-numel(bndedg_);
edgalt = zeros(size(edg));
edgnds = zeros(dim*2,edgs);
edgnds(1:dim,:) = reshape(tri(edg2tri(Iedg,1),:)',dim,edgs);
edgnds(1+dim:2*dim,:) = reshape(tri(edg2tri(Iedg,2),:)',dim,edgs);
edgnds(edgnds==repmat(edg(Iedg,1)',2*dim,1)) = 0;
edgnds(edgnds==repmat(edg(Iedg,2)',2*dim,1)) = 0;
if dim==4
edgnds(edgnds==repmat(edg(Iedg,3)',2*dim,1)) = 0;
end;
edgnds = reshape(edgnds(edgnds~=0),2,edgs)'; edgalt = zeros(numel(Iedg),2); edgalt(Iedg,:) = edgnds;
%clc; clear all; close all;
%options = gen_options(); [xy,tri,bndmesh,options,geomfunc] = mesh_rect(20,options,0);
%metric = @(x_) metric_uniform(x_, [0.1 0.1]); Nmetric = metric(xy);
%[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
%delete('test.pvd'); delete('test0.vtu'); export_vtk('test.pvd',tri,xy,sin(xy),1,{'x1','y2'})

function [] = export_vtk(flname,tri,xy,sclrvr,ii,vrnm)
%%% INPUT %%%
%flname   : string specifying the output pvd filename
%tri      : element list, N x 3 or N x 4 (2D or 3D), where N is the number of elements
%xy       : node coordinates, M x 2 or M x 3 (2D or 3D), where M is the number of elements
%sclrvr   : scalar variables, M x Q or N x Q (node or element wise values), where Q is the number of variables
%ii       : the output number in the series
%vrnm     : the names of the scalar variables
tmp = pwd;
if tmp(1) == '/'
 mslsh = '/';
else
 mslsh = '\';
end;
fid = fopen(flname,'w');
output1 = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n  <Collection>\n';
output2 = [];
flstrt = max([0 find(flname==mslsh)])+1;
for i=1:ii
output2 = [output2 sprintf('    <DataSet timestep="%d" part="0" file="%s%d.vtu" />\n', i-1, flname(flstrt:end-4), i-1)];
end;
output3 = '  </Collection>\n</VTKFile>';
output = [output1 output2 output3];
fwrite(fid,sprintf(output),'char');
fclose(fid);

tabchar = '  ';
fid = fopen(sprintf('%s%d.vtu',flname(1:end-4),ii-1),'w');
output1 = sprintf('<?xml version="1.0"?>\\n<VTKFile type="UnstructuredGrid"  version="0.1"  >\\n<UnstructuredGrid>\\n<Piece  NumberOfPoints="%d" NumberOfCells="%d">\\n<Points>\\n<DataArray  type="Float64"  NumberOfComponents="%d"  format="ascii">',size(xy,1),size(tri,1),3);
if size(xy,2) == 2
 output2 = sprintf(['%.12e %.12e %.12e' tabchar],[xy zeros(size(xy,1),1)]');
 output4 = sprintf(['%d %d %d' tabchar],tri'-1); 
else %==3
 output2 = sprintf(['%.12e %.12e %.12e' tabchar],xy');
 output4 = sprintf(['%d %d %d %d' tabchar],tri'-1);
end;
output6 = sprintf('%d ',size(tri,2):size(tri,2):numel(tri));
output8 = sprintf('%d ',(5+5*(size(xy,2)==3))*ones(size(tri,1),1));
output3 = '</DataArray>\n</Points>\n<Cells>\n<DataArray  type="UInt32"  Name="connectivity"  format="ascii">';
if size(sclrvr,1) == size(xy,1)
output9 = '</DataArray>\n</Cells>\n<PointData>\n';
output10 = [];
for jj=1:size(sclrvr,2)
 output10 = [output10 sprintf('<DataArray  type="Float64"  Name="%s"  format="ascii">',vrnm{jj}) sprintf(['%.12e' tabchar], sclrvr(:,jj)) '</DataArray>\n'];
end;
output11 = '</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>';
else
output9 = '</DataArray>\n</Cells>\n<CellData>\n';
output10 = [];
for jj=1:size(sclrvr,2)
 output10 = [output10 sprintf('<DataArray  type="Float64"  Name="%s"  format="ascii">',vrnm{jj}) sprintf(['%.12e' tabchar], sclrvr(:,jj)) '</DataArray>\n'];
end;
output11 = '</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>';	
end;
output5 = '</DataArray>\n<DataArray  type="UInt32"  Name="offsets"  format="ascii">';
output7 = '</DataArray>\n<DataArray  type="UInt8"  Name="types"  format="ascii">';
output = [output1 output2 output3 output4 output5 output6 output7 output8 output9 output10 output11];
fwrite(fid,sprintf(output),'char');
fclose(fid);
%clc; clear all; close all;
%options = gen_options(); [xy,tri,bndmesh,options,geomfunc] = mesh_rect(20,options,0);
%metric = @(x_) metric_uniform(x_, [0.1 0.1]); Nmetric = metric(xy);
%[tri,xy,Nmetric,bndmesh,triQ] = adapt_mesh(tri,xy,Nmetric,bndmesh,geomfunc,options);
%delete('test.pvd'); delete('test0.vtu'); export_vtk_binary('test.pvd',tri,xy,sin(xy),1,{'x1','y2'})

function [] = export_vtk_binary(flname,tri,xy,sclrvr,ii,vrnm)
%%% INPUT %%%
%flname   : string specifying the output pvd filename
%tri      : element list, N x 3 or N x 4 (2D or 3D), where N is the number of elements
%xy       : node coordinates, M x 2 or M x 3 (2D or 3D), where M is the number of elements
%sclrvr   : scalar variables, M x Q or N x Q (node or element wise values), where Q is the number of variables
%ii       : the output number in the series
%vrnm     : the names of the scalar variables
ndn = 'l'; %LittleEndian ('b' for BigEndian)
tmp = pwd; 
if tmp(1) == '/'
 mslsh = '/';
else
 mslsh = '\';
end;
fid = fopen(flname,'w');
output1 = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n  <Collection>\n';
output2 = [];
flstrt = max([0 find(flname==mslsh)])+1;
for i=1:ii
output2 = [output2 sprintf('    <DataSet timestep="%d" part="0" file="%s%d.vtu" />\n', i-1, flname(flstrt:end-4), i-1)];
end;
output3 = '  </Collection>\n</VTKFile>';
output = [output1 output2 output3];
fwrite(fid,sprintf(output),'char',ndn);
fclose(fid);

fid = fopen(sprintf('%s%d.vtu',flname(1:end-4),ii-1),'w');
if ndn ~= 'l'
fwrite(fid,sprintf('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid"  version="0.1"  byte_order="BigEndian">\n<UnstructuredGrid>\n<Piece  NumberOfPoints="%d" NumberOfCells="%d">\n',size(xy,1),size(tri,1)),'char',ndn);
else
fwrite(fid,sprintf('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid"  version="0.1"  byte_order="LittleEndian">\n<UnstructuredGrid>\n<Piece  NumberOfPoints="%d" NumberOfCells="%d">\n',size(xy,1),size(tri,1)),'char',ndn);
end;
offset = 0;
if numel(sclrvr)~=0
if size(sclrvr,1) == size(xy,1)
fwrite(fid,sprintf('<PointData>\n'),'char',ndn);
output10 = [];
for jj=1:size(sclrvr,2)
 fwrite(fid,sprintf('<DataArray  type="Float32"  Name="%s"  format="appended" offset="%d"/>\n',vrnm{jj},offset),'char',ndn);
 offset = offset + size(xy,1)*4+4;
end;
fwrite(fid,sprintf('</PointData>\n'),'char',ndn);
else
fwrite(fid,sprintf('%s','<CellData>\n'),'char',ndn);
for jj=1:size(sclrvr,2)
 fwrite(fid,sprintf('<DataArray  type="Float32"  Name="%s"  format="appended" offset="%d"/>\n',vrnm{jj},offset),'char',ndn);
 offset = offset + size(tri,1)*4+4;
end;
fwrite(fid,sprintf('</CellData>\n'),'char',ndn);
end;
end;
fwrite(fid,sprintf('<Points>\n<DataArray  type="Float32"  NumberOfComponents="%d"  format="appended" offset="%d"/>',3,offset),'char',ndn);
offset = offset + size(xy,1)*3*4+4;
fwrite(fid,sprintf('\n</Points>\n<Cells>\n<DataArray  type="UInt32"  Name="connectivity"  format="appended" offset="%d"/>',offset),'char',ndn);
offset = offset + numel(tri)*4+4;
fwrite(fid,sprintf('\n<DataArray  type="UInt32"  Name="offsets"  format="appended" offset="%d"/>',offset),'char',ndn);
offset = offset + size(tri,1)*4+4;
fwrite(fid,sprintf('\n<DataArray  type="UInt8"  Name="types"  format="appended" offset="%d"/>',offset),'char',ndn);
offset = offset + size(tri,1)+4;
fwrite(fid,sprintf('\n</Cells>\n'),'char',ndn);

fwrite(fid,sprintf('</Piece>\n</UnstructuredGrid>\n<AppendedData encoding="raw">\n   _'),'char',ndn);
tmpdt = 'int32';
for jj=1:size(sclrvr,2)
 fwrite(fid,size(sclrvr,1)*4,tmpdt,ndn);
 fwrite(fid,sclrvr(:,jj),'float32',ndn);
end;
fwrite(fid,size(xy,1)*3*4,tmpdt,ndn);
if size(xy,2) == 2
 fwrite(fid,reshape([xy zeros(size(xy,1),1)]',numel(xy)*3/2,1),'float32',ndn);
else %==3
 fwrite(fid,reshape(xy',numel(xy),1),'float32',ndn);
end;
fwrite(fid,numel(tri)*4,tmpdt,ndn);
fwrite(fid,reshape(tri'-1,numel(tri),1),'uint32',ndn);
fwrite(fid,size(tri,1)*4,tmpdt,ndn);
fwrite(fid,size(tri,2):size(tri,2):numel(tri),'uint32',ndn);
fwrite(fid,size(tri,1),tmpdt,ndn);
fwrite(fid,(5+5*(size(xy,2)==3))*ones(size(tri,1),1),'uint8',ndn);
fwrite(fid,sprintf('\n  </AppendedData>\n</VTKFile>'),'char',ndn);	
fclose(fid);
function [xy,tri,bndmesh,options,geomfunc] = mesh_rect(N,options,ustruct)
if nargin ~= 2
x = rand(1,N*N);
y = rand(1,N*N);
x = [x 0 0 1 1 zeros(1,N) rand(1,N) ones(1,N) rand(1,N)]';
y = [y 0 1 1 0 rand(1,N) ones(1,N) rand(1,N) zeros(1,N)]';
tri = delaunay(x,y);
tri = elem_fix(tri,[x y]);
else
if numel(N) == 1
N = [N N];
end;
[x,y] = meshgrid(0:1/N(1):1,0:1/N(2):1); x=x(:); y=y(:);
[xn,yn] = meshgrid(1:N(1),1:N(2)); xn=xn(:); yn=yn(:);
trin = yn+(N(2)+1)*(xn-1); 
trin = [trin trin+1 trin+N(2)+1 trin+N(2)+2];
tri = [trin(:,[2 1 3]); trin(:,2:4)];
end;
[C,I] = sort(rand(size(tri,1),1)); tri = tri(I,:); %shuffle elements
[C,I] = sort(rand(size(x))); x = x(I); y = y(I); %shuffle nodes
[C,I2] = sort(I); tri = I2(tri);
xy = [x y];
bndmesh = [];
options.area = 1.;
if nargout == 5
	geomfunc = [];
	%geomfunc = @(x_) geom_rect(x_,1,1); 
end;
function [xy,tri,bndmesh,options,geomfunc] = mesh_box(N,options,ustruct)
if numel(N) == 1
 N = repmat(N,1,3);
end;
if ustruct ~= 0
x = rand(1,N^3); y = rand(1,N^3); z = rand(1,N^3);
%add edges
x = [x 0 0 0 0 1 1 1 1];
y = [y 0 0 1 1 0 0 1 1];
z = [z 0 1 0 1 0 1 0 1];
x = [x rand(1,(N-1)*4)];
y = [y ones(1,N-1) ones(1,N-1)  zeros(1,N-1) zeros(1,N-1)];
z = [z ones(1,N-1) zeros(1,N-1) zeros(1,N-1) ones(1,N-1)];
x = [x ones(1,N-1) ones(1,N-1)  zeros(1,N-1) zeros(1,N-1)];
y = [y rand(1,(N-1)*4)];
z = [z ones(1,N-1) zeros(1,N-1) zeros(1,N-1) ones(1,N-1)];
x = [x ones(1,N-1) ones(1,N-1)  zeros(1,N-1) zeros(1,N-1)]';
y = [y ones(1,N-1) zeros(1,N-1) zeros(1,N-1) ones(1,N-1)]';
z = [z rand(1,(N-1)*4)]';
if ustruct == 2
%%add faces
%[XYZ,YZX] = meshgrid(1/N:1/N:1-1/N); XYZ = repmat(XYZ(:),2,1); YZX = repmat(YZX(:),2,1); XYZ0 = [zeros((N-1)^2,1); ones((N-1)^2,1)];
XYZ = rand((N-1)^2,1); YZX = rand((N-1)^2,1); XYZ = repmat(XYZ(:),2,1); YZX = repmat(YZX(:),2,1); XYZ0 = [zeros((N-1)^2,1); ones((N-1)^2,1)];
x = [x; XYZ0; YZX; XYZ]; 
y = [y; XYZ; XYZ0; YZX]; 
z = [z; YZX; XYZ; XYZ0];
end;
tri = delaunay3(x,y,z);
tri = elem_fix(tri,[x y z]);
else
[x,y,z] = meshgrid(0:1/N(1):1,0:1/N(2):1,0:1/N(3):1); x=x(:); y=y(:); z=z(:);
[xn,yn,zn] = meshgrid(1:N(1),1:N(2),1:N(3)); xn=xn(:); yn=yn(:); zn=zn(:);
trin = yn+(N(2)+1)*(xn-1)+(N(1)+1)*(N(2)+1)*(zn-1);
trin = [trin trin+1 trin+N(2)+1 trin+N(2)+2];
trin = [trin trin+(N(1)+1)*(N(2)+1)];
%tri = [trin(:,[7 3 1 4]); trin(:,[2 8 6 5]); trin(:,[8 1 4 7]); trin(:,[8 1 7 5]); trin(:,[8 1 2 4]); trin(:,[8 1 5 2])];
tri = [trin(:,[7 3 1 4]); trin(:,[2 8 6 5]); trin(:,[4 5 7 8]); trin(:,[4 5 1 7]); trin(:,[4 5 2 1]); trin(:,[4 5 8 2])];
end;
[C,I] = sort(rand(size(tri,1),1)); tri = tri(I,:); %shuffle elements
[C,I] = sort(rand(size(x))); x = x(I); y = y(I); z = z(I); %shuffle nodes
[C,I2] = sort(I); tri = I2(tri);

xy = [x y z];
bndmesh = [];
options.area = 1;
if nargout == 5
	geomfunc = [];
end;
function metric = metric_uniform(x_, dxy);
n = size(x_,1);
metric = [1./dxy(1) 0 1./dxy(2)];
if size(x_,2) == 3
	metric = [metric 0 0 1./dxy(3)];
end;
metric = repmat(metric,n,1);

function Nmetric = metric_pnorm(tri,xy,z,eta,pnorm,options,nd2tri)
%p-norm scaling to the metric, as in Chen, Sun and Xu, Mathematics of
%  Computation, Volume 76, Number 257, January 2007, pp. 179-204.
%%% INPUT %%%
%tri        : element list, N x 3 or N x 4 (2D or 3D), where N is the number of elements
%xy         : node coordinates, M x 2 or M x 3 (2D or 3D), where M is the number of nodes
%z          : scalar nodal variable, M x 1
%eta        : scaling factor
%pnorm      : norm in which to minimize interpolation error
%options    : for controlling the solver uses for galerkin projection, see fem_tri2xy
%nd2tri     : inverse element list facilitates derivative recovery by means of taking the area/volume 
%             weighted average in the elements involving the node.
%             Alternatively the 7th argument can be the output of fem_getX (doproj == false)
%%% OUTPUT %%%
%Nmetric    : scalar and nodal symmetric positive definite tensor, M x 3 or M x 6 (2D or 3D)

doproj = true; %use fem_tri2xy to project DG to CG
if nargin < 7 && not(doproj)
 nd2tri = inv_table(tri);
 nd2tri_ = nd2tri~=0;
end;
dim = size(xy,2);
xys = size(xy,1);
tris = size(tri,1);
triz = z(tri(:,2:end)) - repmat(z(tri(:,1)),1,dim);
if doproj && nargin==7
X = nd2tri;
[I_,dem] = elem_inv(tri,xy);
else
[X,dem] = fem_getX(tri,xy);
end;
if not(doproj) %dem strictly also unnecessary
nd2trim = zeros(size(nd2tri)); nd2trim(nd2tri_) = dem(nd2tri(nd2tri_));
nd2trim = nd2trim./repmat(sum(nd2trim,2),1,size(nd2tri,2));
end;
gradf = zeros(size(triz));
hes = zeros(tris,dim*dim);
hesn = zeros(xys,dim*dim);
for i=1:dim
for j=1:dim
	gradf(:,i) = gradf(:,i) + triz(:,j).*X(:,dim*(j-1)+i);
end;
if doproj
gradfn_ = fem_tri2xy(tri,xy,gradf(:,i),options);
else
gradfn_ = zeros(size(nd2tri)); gradfn_(nd2tri_) = gradf(nd2tri(nd2tri_),i);
gradfn_ = sum(gradfn_.*nd2trim,2);
end;
triz_ = gradfn_(tri(:,2:end)) - repmat(gradfn_(tri(:,1)),1,dim);
for j=1:dim
for k=1:dim
	hes(:,(i-1)*dim+j) = hes(:,(i-1)*dim+j) + triz_(:,k).*X(:,dim*(k-1)+j);
end;
if doproj
hesn(:,(i-1)*dim+j) = fem_tri2xy(tri,xy,hes(:,(i-1)*dim+j),options);
else
hesn_ = zeros(size(nd2tri)); hesn_(nd2tri_) = hes(nd2tri(nd2tri_),(i-1)*dim+j);
hesn(:,(i-1)*dim+j) = sum(hesn_.*nd2trim,2);
end;
end;
end;

if size(xy,2) == 2
	hesn(:,2) = mean(hesn(:,2:3),2);
	hesn = hesn(:,[1 2 4]);
    hesn(:,[1 3]) = hesn(:,[1 3]) + 1e-12;
else
	hesn(:,2) = mean(hesn(:,[2 4]),2);
	hesn(:,3) = mean(hesn(:,[3 7]),2);
	hesn(:,6) = mean(hesn(:,[6 8]),2);
	hesn = hesn(:,[1 2 5 3 6 9]);
    hesn(:,[1 3 6]) = hesn(:,[1 3 6]) + 1e-12;
end;
[eigL,eigR] = analyt_eig(hesn);
eigL = reshape(max(abs(eigL),1e-12),size(eigL));
detn = prod(eigL,2);
hesn = analyt_prod(analyt_fulleig(eigL.^0.5),eigR);
exponent = -0.5/(2.*pnorm + size(xy,2));
Nmetric = 1/sqrt(eta)*repmat(detn.^exponent,1,size(hesn,2)).*hesn;
function mmI = metric_avg(edg,Nmetric,options);
if options.log_m_add
 if size(edg,2) == 2
   mmI = add_log_metric(Nmetric(edg(:,1),:),Nmetric(edg(:,2),:));
 elseif size(edg,2) == 3
   mmI = add_log_metric(Nmetric(edg(:,1),:),Nmetric(edg(:,2),:),Nmetric(edg(:,3),:));
 else
   mmI = add_log_metric(Nmetric(edg(:,1),:),Nmetric(edg(:,2),:),Nmetric(edg(:,3),:),Nmetric(edg(:,4),:));
 end;
else
 if size(edg,2) == 2
  mmI = (Nmetric(edg(:,1),:)+Nmetric(edg(:,2),:))/2;
 elseif size(edg,2) == 3
  mmI = (Nmetric(edg(:,1),:)+Nmetric(edg(:,2),:)+Nmetric(edg(:,3),:))/3;
 else
  mmI = (Nmetric(edg(:,1),:)+Nmetric(edg(:,2),:)+Nmetric(edg(:,3),:)+Nmetric(edg(:,4),:))/4;
 end;
end;


function out_sum = add_log_metric(metric1,metric2,metric3,metric4)
%[lambda1 ,lambda2 ,v1xn ,v1yn ] = analyt_eig(metric1);
%mm_log1 = [v1xn.^2.*log(lambda1)+v1yn.^2.*log(lambda2) ...
       %v1yn.*v1xn.*(log(lambda1)-log(lambda2)) ...
       %v1yn.^2.*log(lambda1)+v1xn.^2.*log(lambda2)];
%[lambda1,lambda2,v1xn,v1yn] = analyt_eig(metric2);
%mm_log2 = [v1xn.^2.*log(lambda1)+v1yn.^2.*log(lambda2) ...
       %v1yn.*v1xn.*(log(lambda1)-log(lambda2)) ...
       %v1yn.^2.*log(lambda1)+v1xn.^2.*log(lambda2)];
%logsum = mm_log1+mm_log2;
%if nargin == 3
	%[lambda1,lambda2,v1xn,v1yn] = analyt_eig(metric3);
	%logsum = logsum+[v1xn.^2.*log(lambda1)+v1yn.^2.*log(lambda2) ...
	       %v1yn.*v1xn.*(log(lambda1)-log(lambda2)) ...
	       %v1yn.^2.*log(lambda1)+v1xn.^2.*log(lambda2)];
%end;
%logsum = logsum/nargin;
%[lambda1,lambda2,v1xn,v1yn] = analyt_eig(logsum);
%out_sum = [v1xn.^2.*exp(lambda1)+v1yn.^2.*exp(lambda2) ...
       %v1yn.*v1xn.*(exp(lambda1)-exp(lambda2)) ...
       %v1yn.^2.*exp(lambda1)+v1xn.^2.*exp(lambda2)];
if size(metric1,2) == 1 %1D
 outsum = exp((log(metric1)+log(metric2))/2.);
 return;
end;
[eigL, eigR] = analyt_eig(metric1);
mm_log1 = analyt_prod(analyt_fulleig(log(eigL)),eigR);
[eigL, eigR] = analyt_eig(metric2);
mm_log2 = analyt_prod(analyt_fulleig(log(eigL)),eigR);
logsum = mm_log1+mm_log2;
if nargin == 3
	[eigL, eigR] = analyt_eig(metric3);
	logsum = logsum+analyt_prod(analyt_fulleig(log(eigL)),eigR);
end;
if nargin == 4
	[eigL, eigR] = analyt_eig(metric4);
	logsum = logsum+analyt_prod(analyt_fulleig(log(eigL)),eigR);
end;
logsum = logsum/nargin;
[eigL, eigR] = analyt_eig(logsum);
out_sum = analyt_prod(analyt_fulleig(exp(eigL)),eigR);
function dist = geom_circ(x_,r,c) %or sphere
if nargin == 1
	r=1;
end;
if nargin < 3
	c=repmat(0,1,size(x_,2));
end;
x_ = x_-repmat(c,size(x_,1),1);
rr = sqrt(sum(x_.^2,2)+eps);
dir_r = -x_./repmat(rr,1,size(x_,2)); 
dir_r(rr < eps,:) = 0;
dist = [(r-rr) ...
          dir_r];

function [crnds,edg,edg2,edg2tri,edg2ID] = geom_crnds(bndmesh,addcrs)
if not(isfield(bndmesh,'fac'))
nds = bndmesh.edg(:);
IDs = repmat(bndmesh.IDs,2,1);
[nds,I] = sort(nds);
IDs = IDs(I);
I2 = [false; and(nds(1:end-1)==nds(2:end),IDs(1:end-1)~=IDs(2:end))];
crnds = nds(I2);
%I3_ = find(I3);
%snglnds = nds(I3_(diff([0; I3_]) == 1));
snglnds = nds([nds(1) ~= nds(2); and(nds(1:end-2)~=nds(2:end-1),nds(2:end-1)~=nds(3:end)); nds(end-1) ~= nds(end)]);
crnds = [crnds; snglnds]; %an edge can stop in 2D space, but we do not check for isolated nodes
%bndmesh.crnds = [bndmesh.crnds; crnds];
else %3D
nds = bndmesh.fac(:);
IDs = repmat(bndmesh.IDs,3,1);
%get corner nodes
[edg2,edg2ID,edg2tri] = books(bndmesh.fac,bndmesh.IDs);
I = edg2ID(:,1)~=edg2ID(:,2);
if size(edg2ID,2) > 2
I = or(I,sum(edg2ID~=0,2)>2);
end;
edg = edg2(I,:);
if nargin == 2
%edg2ID = edg2ID(I,:);
%edg2tri = edg2tri(I,:);
edgs = [1:size(edg,1) 1:size(edg,1)]';
[nds,I] = sort(edg(:));
%IDs = IDs(I);
edgs = edgs(I);
I2 = [false; and(and(nds(2:end-1)==nds(3:end),edgs(2:end-1)~=edgs(3:end)),and(nds(2:end-1)==nds(1:end-2),edgs(2:end-1)~=edgs(1:end-2))); false];
crnds = nds(I2);
if numel(crnds) ~=0
crnds = crnds([crnds(1:end-1)~=crnds(2:end); true]);
else
crnds = [];
end;
else
crnds = [];
end; %nargin==2 
end; %3D
 %a face can stop in 3D space, but we do not check for isolated nodes or lines going through 3D space

function [edg,edg2ID,edg2tri] = books(tria,IDs)
edga = [tria(:,2) tria(:,1); tria(:,3) tria(:,2); tria(:,1) tria(:,3)];
edga2tri = [1:size(tria,1) 1:size(tria,1) 1:size(tria,1)]';
IDs = repmat(IDs,3,1);
edg = sort(edga,2); %biggest in last column
[edg,Is] = sortrows(edg);
tris = edga2tri(Is);
IDs = IDs(Is);
d =[true; or(edg(1:end-1,1)~=edg(2:end,1), ...
            edg(1:end-1,2)~=edg(2:end,2))]; Nd = sum(d);
edg = edg(d,:);
edgs = cumsum(d);
edg2tri = rpval2M(edgs,tris); edg2tri = rpval2M_clean(edg2tri);
edg2ID = rpval2M(edgs,IDs); edg2tri = rpval2M_clean(edg2tri);
function out = geom_diff(in1,in2)
inside   = and(in1(:,1)>=0,in2(:,1)<=0);
outsideA = and(in1(:,1)< 0,in2(:,1) <0);
outsideB = in2(:,1)>0;
L1 = abs(in1(:,1))<abs(in2(:,1)); L2 = not(L1);
out = zeros(size(in1));
insideL1 = and(inside,L1);
insideL2 = and(inside,L2);
out(insideL1 ,:) =  in1(insideL1,:);
out(insideL2 ,:) = -in2(insideL2,:);
out(outsideA,:) =  in1(outsideA,:);
out(outsideB,:) = -in2(outsideB,:);
%insideA1 = and(insideA,L1);
%insideA2 = and(insideA,L2);
%out(insideA1,:) = in2(insideA1,:);
%out(insideA2,:) = in1(insideA2,:);
function dist = geom_rect(x_,dxy,poscr)
if nargin == 1
	dxy = repmat(1,1,size(x_,2));
end;
if nargin < 3
	poscr = repmat(0,1,size(x_,2));
end;
x_ = x_ - repmat(dxy/2+poscr,size(x_,1),1);
LH = dxy(1)/2-abs(x_(:,1));
LV = dxy(2)/2-abs(x_(:,2));
chooseH = LH < LV;
chooseV = not(chooseH);
if size(x_,2) == 3
LD = dxy(3)/2-abs(x_(:,3));
chooseH = and(chooseH,LH <= LD);
chooseV = and(LV <= LH, LV < LD);
chooseD = and(not(chooseH),not(chooseV));
dist = [chooseH.*LH+ ...
        chooseV.*LV+ ...
        chooseD.*LD  ...
        (-sign(x_(:,1)).*chooseH) ... 
        (-sign(x_(:,2)).*chooseV) ...
        (-sign(x_(:,3)).*chooseD)]; %dist vecx vecy vecz
else
dist = [chooseH.*LH+...
        chooseV.*LV ...
        (-sign(x_(:,1)).*chooseH) ... 
        (-sign(x_(:,2)).*chooseV)]; %dist vecx vecy
end;
function dist = geom_rotate(x_,func,angle,origin)
if nargin == 3
origin = zeros(1,size(x_,2));
end;
x_ = x_ - repmat(origin,size(x_,1),1);
if size(x_,2) == 2
x_ = [cos(angle)*x_(:,1) + sin(angle)*x_(:,2), -sin(angle)*x_(:,1) + cos(angle)*x_(:,2)];
else
x_ = [cos(angle)*x_(:,1) + sin(angle)*x_(:,2), -sin(angle)*x_(:,1) + cos(angle)*x_(:,2) x_(:,3)];
%warning('experimental 3D rotation');
end;
dist = func(x_+repmat(origin,size(x_,1),1));
dist(:,2:3) = [cos(angle)*dist(:,2) - sin(angle)*dist(:,3), sin(angle)*dist(:,2) + cos(angle)*dist(:,3)];
function invtbl = inv_table(table_)
[nds,I] = sort(table_(:));
edgs = repmat(1:size(table_,1),1,size(table_,2))'; edgs = edgs(I);
invtbl = rpval2M(nds,edgs);
function out = rpval2M(R,vals)
%C    = cell2mat(accumarray(R,1   ,[],@(x){cumsum(x)})); nMax = max(C);
%vals = cell2mat(accumarray(R,vals,[],@(x){x}));
%R    = cell2mat(accumarray(R,R   ,[],@(x){x}));
%out = zeros(Rmax,nMax); out(R+(C-1)*Rmax) = vals;

%faster alternative using sort()
if not(issorted(R,'rows'))
	[R,I] = sort(R); %[size(vals) size(I)]
	vals = vals(I);
end;
nN = diff([0 find([R(1:end-1)'~=R(2:end)' true])]); nMax = max(nN); Rmax = max(R);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN,nMax,1); I = C <= CnN; 
out = zeros(Rmax,nMax); out(R+(reshape(C(I),size(R))-1)*Rmax) = vals;
function out = rpval2M_clean(vin)
% takes out doubles
vin = sort(vin,2);
I = and([true(size(vin,1),1) vin(:,1:end-1)~=vin(:,2:end)],vin~=0);
vin(not(I)) = 0;
vin = sort(vin,2);
nN = sum(vin~=0,2);
out = vin(:,size(vin,2):-1:size(vin,2)-max(nN)+1);
