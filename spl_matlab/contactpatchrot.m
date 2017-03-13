function [phi] = contactpatchrot(pe,pc,varargin)
% [phi] = contactpatchrot(pf,pe) computes the rotational angles of valid contact
% between the foot and the environment patches.
%
%   INPUTS
%
%   The parameters pf and pe are (scalar) Matlab structs describing the foot and
%   environment patches respectively.  See the documentation for the function
%   patchchk() for details.
%
%   OUTPUTS
%
%   The output phi is the set of the ccw angles (in rad) around the z local axis
%   of the foot patch pf, that are valid for contact.  The 0 rad angle is
%   the one of the patch's local axis (x-up, y-left, z-out).  The phi Na-by-2
%   matrix contains the starting (first column) and the ending (second column)
%   angles in radians.
%
%   OPTIONS
%
%   A list of (name,value) pairs may be given to set options.  Unrecognized 
%   names cause warnings.  If a name is given more than once the last-given 
%   (rightmost) value takes precedence.
%
%   Option 'dbg' (default 0): may have any nonnegative value.  dbg=0 disables
%   debug.  Larger values enable successively more verbose debug.
%
%
%   TBD
%   1. Refactor the code to take advantage of 
%
% Copyright (C) 2016 Dimitrios Kanoulas

tstart = tic(); % timing

% process (name,value) options
nva = length(varargin);
if (mod(nva,2)~=0)
  warning('expect even number of varargs (name,value pairs)');
  nva = nva-1;
end
nopt = nva/2;

% option defaults
dbg = 0;

for iopt=1:nopt
  n = varargin{2*iopt-1}; v = varargin{2*iopt};
  if (ischar(n))
    switch (n)
      case 'dbg'; dbg = v;
      otherwise; warning('unexpected optional arg %s',n);
    end
  else
    warning('non-string, expected name');
  end
end

% Read or infer surface and boundary type
if (isfield(pc,'b')); btf = pc.b; end;
if (isfield(pe,'b')); bte = pe.b; end;

% If the foot patch is a planar rect
if strcmp('r',btf) 
  if (dbg); fprintf('Foot patch is a planar rect.\n'); end

  rx = pc.d(1); ry = pc.d(2); % dx and dy radii size 
  rf = sqrt(rx^2+ry^2); % the hypotenuse of the radii (for the bounded circle)
  
  % vertices of the foot rect in 3D in the local frame
  ve1 = [ rx -ry true]; ve2 = [ rx  ry true]; % upper right/left
  ve3 = [-rx -ry true]; ve4 = [-rx  ry true]; % lower left/right
  ve = [ve1 0; ve2 0; ve3 0; ve4 0];
  
  ipts=[]; isec=0; %if boundaries intersect
  switch (bte) % environment patch boundary type
    case 'c' % circle
      phi = [0 2*pi]; isec = 0;
    case 'e' % ellipse
      if (dbg); disp ('Env patch is an ellipse'); end
      
      % intersection between ellipse and circle
      % see: http://mathworld.wolfram.com/Circle-EllipseIntersection.html
      re = min(pe.d(1),pe.d(2));
      if (rf<=re)
        phi = [0 2*pi]; isec = 0;
      else
        iptx = pe.d(1) * sqrt(((rf^2) - (pe.d(2)^2))/((pe.d(1)^2)-(pe.d(2)^2)));
        ipty = pe.d(2) * sqrt(((pe.d(1)^2) - (rf^2))/((pe.d(1)^2)-(pe.d(2)^2)));
        
        % intersection points
        % ipts: x,y,vertex(1)/intersection(0), angle from angle 0
        radius = sqrt(rx^2+ry^2);
        ipts_start = [radius 0 true]; %angle starting point
        ipts_ul = [iptx ipty false]; ipts_ll = [-iptx ipty false]; %upper/lower left
        ipts_lr = [-iptx -ipty false]; ipts_ur = [iptx -ipty false]; %lower/upper right
        ipts = [ipts_start; ipts_ul; ipts_ll; ipts_lr; ipts_ur]; %ccw order
        ipts = [ipts; ve1; ve2; ve3; ve4]; %include all the points
        ipts(:,4) = mod(atan2(ipts(:,2),ipts(:,1)),2*pi); %compute all the angles
        ipts = sortrows(ipts,4); %sorted wrt the angle from C_0
        
        % mark the valid angle for the first vertex
        % First vertex
        idx_v = find (ipts(2:end,3),1,'first'); % find the first vertex instance and shift
        ipts = circshift(ipts,-idx_v);
        ipts = comp_valid_angles (ipts, 5);
        phi_origin(:,1) = mod((ipts(:,4)-ipts(1,4)),2*pi);
        phi_origin(:,2) = ipts(:,5);
        
        % Second vertex
        idx_v = find (ipts(2:end,3),1,'first'); % find the first vertex instance and shift
        ipts = circshift(ipts,-idx_v);
        ipts = comp_valid_angles (ipts,6);
        phi_origin(:,3) = mod((ipts(:,4)-ipts(1,4)),2*pi);
        phi_origin(:,4) = ipts(:,6);
        
        % Third vertex
        idx_v = find (ipts(2:end,3),1,'first'); % find the first vertex instance and shift
        ipts = circshift(ipts,-idx_v);
        ipts = comp_valid_angles (ipts,7);
        phi_origin(:,5) = mod((ipts(:,4)-ipts(1,4)),2*pi);
        phi_origin(:,6) = ipts(:,7);
        
        % Fourth vertex
        idx_v = find (ipts(2:end,3),1,'first'); % find the first vertex instance and shift
        ipts = circshift(ipts,-idx_v);
        ipts = comp_valid_angles (ipts,8);
        phi_origin(:,7) = mod((ipts(:,4)-ipts(1,4)),2*pi);
        phi_origin(:,8) = ipts(:,8);
                
        % all the angles splits
        phi_all = uniquetol([phi_origin(:,1);
                             phi_origin(:,3);
                             phi_origin(:,5);
                             phi_origin(:,7)]);
        phi_all_rows = size(phi_all,1);
        epsilon = 0.00001;
        for i=1:phi_all_rows
          [r,~] = find(phi_origin(:,1)<=phi_all(i,1)+epsilon,1,'last');
          phi_all (i,2) = phi_origin(r,2);
          
          [r,~] = find(phi_origin(:,3)<=phi_all(i,1)+epsilon,1,'last');
          phi_all (i,3) = phi_origin(r,4);
          
          [r,~] = find(phi_origin(:,5)<=phi_all(i,1)+epsilon,1,'last');
          phi_all (i,4) = phi_origin(r,6);
          
          [r,~] = find(phi_origin(:,7)<=phi_all(i,1)+epsilon,1,'last');
          phi_all (i,5) = phi_origin(r,8);
        end
        phi_all(:,6) = all (phi_all(:,2:5),2);
        
        phi=[];
        for i=1:size(phi_all,1)-1
          if phi_all(i,6)==true
            phi = [phi; phi_all(i,1) phi_all(i+1,1)];
          end
        end
        if phi_all(end,6)==true
          phi = [phi; phi_all(end,1) 2*pi];
        end
        
        % Compute the z-value of the intersection points
        iptsz = pe.zl(ipts(:,1),ipts(:,2)); isec = 1;
      end
      
    case {'r','q'} % axis aligned rectangle
      %TBD
      if (dbg); disp ('axis aligned rectangle'); end;
      % vertices of the env rect in 3D in the local frame
      if (strcmp(bte,'r'))
        if (dbg);  disp ('bte: r'); end;
        vpe = [ pe.d(1)  pe.d(2) pe.zl( pe.d(1),  pe.d(2));
               -pe.d(1)  pe.d(2) pe.zl( pe.d(1), -pe.d(2));
               -pe.d(1) -pe.d(2) pe.zl(-pe.d(1), -pe.d(2));
                pe.d(1) -pe.d(2) pe.zl(-pe.d(1),  pe.d(2))];          
      end
      
      if (strcmp(bte,'q'))
        if (dbg); disp ('bte: q'); end;
        g = pe.d(5); p = [g; pi-g; pi+g; -g];
        r = [pe.d(1); pe.d(2); pe.d(3); pe.d(4)];
        v = [cos(p).*r, sin(p).*r];
        vpe = [v(1,1) v(1,2) pe.zl(v(1,1),v(1,2));
               v(2,1) v(2,2) pe.zl(v(2,1),v(2,2));
               v(3,1) v(3,2) pe.zl(v(3,1),v(3,2));
               v(4,1) v(4,2) pe.zl(v(4,1),v(4,2))];
      end
      
      ipts = [];
      new_ipts = line_circle_int (vpe(1,1:2),vpe(2,1:2),0,0,rf);
      ipts = [ipts; new_ipts];
      new_ipts = line_circle_int (vpe(2,1:2),vpe(3,1:2),0,0,rf);
      ipts = [ipts; new_ipts];
      new_ipts = line_circle_int (vpe(3,1:2),vpe(4,1:2),0,0,rf);
      ipts = [ipts; new_ipts];
      new_ipts = line_circle_int (vpe(4,1:2),vpe(1,1:2),0,0,rf);
      ipts = [ipts; new_ipts];

      radius = sqrt(rx^2+ry^2);
      ipts_start = [radius 0 true]; %angle starting point
      ipts = [ipts; ipts_start; ve1; ve2; ve3; ve4]; %include all the points
      ipts(:,4) = mod(atan2(ipts(:,2),ipts(:,1)),2*pi); %compute all the angles
      ipts = sortrows(ipts,4); %sorted wrt the angle from C_0
      
      % mark the valid angle for the first vertex
      % First vertex
      idx_v = find (ipts(2:end,3),1,'first'); % find the first vertex instance and shift
      ipts = circshift(ipts,-idx_v);
      ipts = comp_valid_angles (ipts, 5);
      phi_origin(:,1) = mod((ipts(:,4)-ipts(1,4)),2*pi);
      phi_origin(:,2) = ipts(:,5);
      
      % Second vertex
      idx_v = find (ipts(2:end,3),1,'first'); % find the first vertex instance and shift
      ipts = circshift(ipts,-idx_v);
      ipts = comp_valid_angles (ipts,6);
      phi_origin(:,3) = mod((ipts(:,4)-ipts(1,4)),2*pi);
      phi_origin(:,4) = ipts(:,6);
      
      % Third vertex
      idx_v = find (ipts(2:end,3),1,'first'); % find the first vertex instance and shift
      ipts = circshift(ipts,-idx_v);
      ipts = comp_valid_angles (ipts,7);
      phi_origin(:,5) = mod((ipts(:,4)-ipts(1,4)),2*pi);
      phi_origin(:,6) = ipts(:,7);
      
      % Fourth vertex
      idx_v = find (ipts(2:end,3),1,'first'); % find the first vertex instance and shift
      ipts = circshift(ipts,-idx_v);
      ipts = comp_valid_angles (ipts,8);
      phi_origin(:,7) = mod((ipts(:,4)-ipts(1,4)),2*pi);
      phi_origin(:,8) = ipts(:,8);
      
      % all the angles splits
      phi_all = uniquetol([phi_origin(:,1);
                           phi_origin(:,3);
                           phi_origin(:,5);
                           phi_origin(:,7)]);
      phi_all_rows = size(phi_all,1);
      epsilon = 0.00001;
      for i=1:phi_all_rows
        [r,~] = find(phi_origin(:,1)<=phi_all(i,1)+epsilon,1,'last');
        phi_all (i,2) = phi_origin(r,2);
        
        [r,~] = find(phi_origin(:,3)<=phi_all(i,1)+epsilon,1,'last');
        phi_all (i,3) = phi_origin(r,4);
        
        [r,~] = find(phi_origin(:,5)<=phi_all(i,1)+epsilon,1,'last');
        phi_all (i,4) = phi_origin(r,6);
        
        [r,~] = find(phi_origin(:,7)<=phi_all(i,1)+epsilon,1,'last');
        phi_all (i,5) = phi_origin(r,8);
      end
      phi_all(:,6) = all (phi_all(:,2:5),2);
      
      phi=[];
      for i=1:size(phi_all,1)-1
        if phi_all(i,6)==true
          phi = [phi; phi_all(i,1) phi_all(i+1,1)];
        end
      end
      if phi_all(end,6)==true
        phi = [phi; phi_all(end,1) 2*pi];
      end
      
      % Compute the z-value of the intersection points
      iptsz = pe.zl(ipts(:,1),ipts(:,2)); isec = 1;
    otherwise
      iee('unknown bounds type');
  end
end


% If the foot patch is a sphere
if strcmp('c',btf) 
  if (dbg); fprintf('Foot patch is a sphere.\n'); end
  phi = [0 2*pi];
end

t = toc(tstart);
if (dbg); fprintf('Time for contact search=%g\n',t); end
if (nargout>1); sec = t; end
  
  % function that computes the angles between points
  function phi = phicomp(pts)
  % phi = phicomp(ipts, firstin) computes the angles for which a foot patch 
  % (expressed as its bounding circle) is inside the env patch.
    
  %TBD: refactor the code
  end

  function ipts = comp_valid_angles (ipts, index)
    ipts_size = size(ipts);
    ipts(1,index) = 1;
    for ipts_i=2:ipts_size(1)
      if ipts(ipts_i,3)==0
        ipts(ipts_i,index) = ~ipts(ipts_i-1,index);
      else
        ipts(ipts_i,index) = ipts(ipts_i-1,index);
      end
    end
  end

  function ipts = line_circle_int (lpt1,lpt2,ccx,ccy,cr)
  % ipts = line_circle_int (lpt1,lpt2,ccx,ccy,cr) returns the intersection
  % points between a line segment (lpt1, lpt2) and a circle centered at (ccx,
  % ccy), with radius cr.
  %
  % lpt1,lpt2 are 2D points that characterize the line segment
  % ccx, ccy are the xy coordinates of the circle center
  % cr is the radius of the circle
  slope = (lpt1(2)-lpt2(2))/(lpt1(1)-lpt2(1));
  if isinf(slope) %x-intercept
    intercpt = lpt1(1); %x-value
    if(dbg); plot(intercpt,0,'*'); end
  else %y-intercept
    intercpt = lpt1(2)-(slope*lpt1(1)); %y-value
    if(dbg); plot(0,intercpt,'*'); end
  end
  [xout,yout] = linecirc(slope,intercpt,ccx,ccy,cr);
  if(dbg); refline([slope intercpt]); end
  
  % check if the intersection points are in between the line segment
  ipts = [];
  for iint=1:size(xout,1)+1
    if (~all(isnan(xout)))
      sqdist = norm (lpt1-lpt2);
      sqdist_1 = norm (lpt1-[xout(iint) yout(iint)]);
      sqdist_2 = norm ([xout(iint) yout(iint)]-lpt2);
      if (abs(sqdist-sqdist_1-sqdist_2) <= eps)
        ipts = [ipts; xout(iint) yout(iint) false];
      else
        ipts = [];
      end
    else
      ipts = [];
    end
  end
  end %line_circle_int

end