function [F] = createFaultGridPoints3D(faultTri, rho)

F.f.pts = [];
F.c.CC = [];
F.c.R = [];
F.f.pri = [];
F.f.Gs = [];
F.l.fPos = 1;
F.l.f =[];
F.l.cPos = 1;

for i = 1:numel(faultTri)
  dF = faultTri{i};
  R = rho(dF.Points);

  % cI1 = dF.ConnectivityList(:,1);
  % cI2 = dF.ConnectivityList(:,2:3);
  % R1  = R([cI1;cI2(:,1)]);
  % R2  = R(cI2')';
  % CC1 = dF.Points([cI1;cI2(:,1)],:);
  % CC2 = reshape(dF.Points(cI2',:)',6,[])';

  % cI = dF.ConnectivityList(:);
  % cI = repmat(cI,1,3);
  % cI = [cI(:,1),circshift(cI(:,2),size(cI,1)/3,1),circshift(cI(:,3),-size(cI,1)/3,1)];
  % 
  % 
  % 
  % n = cross(dF.Points(cI2(:,1),:) - dF.Points(cI1,:),dF.Points(cI2(:,2),:) - dF.Points(cI1,:));
  % n = bsxfun(@rdivide,n,sqrt(sum(n.^2,2)));

  CC = reshape(dF.Points(dF.ConnectivityList',:)',9,[])';
  R  = reshape(R(dF.ConnectivityList',:),3,[])';
  [fPts,Gs] = circCircCircInt(CC(:,1:3),R(:,1),CC(:,4:6),R(:,2),CC(:,7:9),R(:,3));
  fPts = uniquetol(fPts,'byrows',true);
  CC = reshape(CC',3,[])';
  F.l.fPos = [F.l.fPos; size(F.f.pts,1)+1+size(fPts,1)];
  F.l.cPos = [F.l.cPos; size(F.c.CC,1)+1+size(CC,1)];
  F.f.pts = [F.f.pts;fPts];
  F.c.CC =[F.c.CC; CC];
  F.c.R  = [F.c.R; reshape(R',[],1)];
  F.f.pri = [F.f.pri;i*ones(size(fPts,1),1)];
  F.f.Gs = [F.f.Gs;Gs];

end
 F.l.f = (1:F.l.fPos(end)-1)';
 F.l.c = (1:F.l.cPos(end)-1)';
end





